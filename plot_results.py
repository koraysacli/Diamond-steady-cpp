#!/usr/bin/env python3
"""Make README figures from a saved iteration: pressure, Mach, velocity,
numerical schlieren (shock lines), surface Cp and convergence history.

Usage: python3 plot_results.py [--iteration N]   (default: latest saved)
Output: docs/images/*.png and docs/results.json
"""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/diamond-steady-mpl')
import argparse, csv, json, math
from pathlib import Path
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.colors import LogNorm

vtk.vtkLogger.SetStderrVerbosity(vtk.vtkLogger.VERBOSITY_WARNING)
root = Path(__file__).resolve().parent
out = root / 'docs' / 'images'
out.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- inputs
cfg = {}
for line in (root / 'case-summary.txt').read_text().splitlines():
    if '=' in line:
        k, v = line.split('=', 1)
        try: cfg[k] = float(v)
        except ValueError: pass
pinf, rhoinf, Uinf = cfg['p_Pa'], cfg['rho'], cfg['U']
q = 0.5 * rhoinf * Uinf**2
g, R = 1.4, 287.05287
chord, half_t = 0.2, 0.01
theta = math.atan(half_t / (chord / 2))

def saved_times(base):
    ts = []
    for d in base.iterdir():
        try: t = float(d.name)
        except ValueError: continue
        if t > 0 and (d / 'p').exists(): ts.append(t)
    return set(ts)

ap = argparse.ArgumentParser(); ap.add_argument('--iteration', type=float); a = ap.parse_args()
recon = saved_times(root)
procs = sorted(root.glob('processor[0-9]*'))
decomp = set.intersection(*[saved_times(p) for p in procs]) if procs else set()
it = a.iteration if a.iteration is not None else max(recon | decomp)
reconstructed = it in recon
if not reconstructed and it not in decomp:
    raise SystemExit(f'Iteration {it} not found')

# ---------------------------------------------------------------- theory
def beta_weak(M, th):
    f = lambda b: math.tan(th) - 2/math.tan(b)*(M*M*math.sin(b)**2-1)/(M*M*(g+math.cos(2*b))+2)
    lo, hi = math.asin(1/M) + 1e-9, math.radians(64)
    for _ in range(200):
        mid = (lo+hi)/2
        if f(lo)*f(mid) <= 0: hi = mid
        else: lo = mid
    return (lo+hi)/2
def after_shock(M, b):
    Mn = M*math.sin(b)
    p = 1 + 2*g/(g+1)*(Mn*Mn-1)
    Mn2 = math.sqrt((1+(g-1)/2*Mn*Mn)/(g*Mn*Mn-(g-1)/2))
    return p, Mn2/math.sin(b-theta)
def nu(M):
    k = math.sqrt((g+1)/(g-1))
    return k*math.atan(math.sqrt((M*M-1)/k**2)) - math.atan(math.sqrt(M*M-1))
def nu_inv(v):
    lo, hi = 1.0001, 20
    for _ in range(200):
        mid = (lo+hi)/2
        if nu(mid) < v: lo = mid
        else: hi = mid
    return (lo+hi)/2
M0 = 3.0
b1 = beta_weak(M0, theta); p1, M1 = after_shock(M0, b1)
M2 = nu_inv(nu(M1) + 2*theta)
p2 = p1 * ((1+(g-1)/2*M1*M1)/(1+(g-1)/2*M2*M2))**(g/(g-1))
b3 = beta_weak(M2, theta)          # trailing-edge shock, relative to flow at -theta
mu1, mu2 = math.asin(1/M1), math.asin(1/M2)
cd_th = 2*(p1-p2)*pinf*half_t / (q*chord)

# ---------------------------------------------------------------- read case
reader = vtk.vtkPOpenFOAMReader()
reader.SetFileName(str(root / 'diamond.foam'))
reader.SetCaseType(1 if reconstructed else 0)
reader.CreateCellToPointOff(); reader.UpdateInformation()
reader.DisableAllPatchArrays()
for n in ('internalMesh', 'patch/diamond'): reader.SetPatchArrayStatus(n, 1)
reader.DisableAllCellArrays()
for n in ('p', 'U', 'T', 'rho'): reader.SetCellArrayStatus(n, 1)
reader.UpdateTimeStep(it)
def leaves(b):
    for i in range(b.GetNumberOfBlocks()):
        c = b.GetBlock(i); name = b.GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME())
        if isinstance(c, vtk.vtkMultiBlockDataSet): yield from leaves(c)
        elif c is not None: yield str(name), c
blocks = dict(leaves(reader.GetOutput()))
mesh, wall = blocks['internalMesh'], blocks['diamond']
cd_ = mesh.GetCellData()
U = vtk_to_numpy(cd_.GetArray('U')); T = vtk_to_numpy(cd_.GetArray('T'))
Umag = np.linalg.norm(U, axis=1)
for name, v in (('Umag', Umag), ('Mach', Umag/np.sqrt(g*R*T))):
    arr = numpy_to_vtk(v, deep=True); arr.SetName(name); cd_.AddArray(arr)
c2p = vtk.vtkCellDataToPointData(); c2p.SetInputData(mesh); c2p.Update()

def sample(x0, x1, y0, y1, nx, ny):
    img = vtk.vtkImageData(); img.SetDimensions(nx, ny, 1)
    img.SetOrigin(x0, y0, 0.0005); img.SetSpacing((x1-x0)/(nx-1), (y1-y0)/(ny-1), 1)
    pr = vtk.vtkProbeFilter(); pr.SetInputData(img); pr.SetSourceData(c2p.GetOutput()); pr.Update()
    o = pr.GetOutput(); pd = o.GetPointData()
    mask = vtk_to_numpy(pd.GetArray('vtkValidPointMask')).reshape(ny, nx) == 0
    return {n: np.ma.array(vtk_to_numpy(pd.GetArray(n)).reshape(ny, nx), mask=mask)
            for n in ('p', 'Mach', 'Umag', 'rho')}

X0, X1, Y0, Y1 = -0.05, 0.55, -0.22, 0.22
NX, NY = 1500, 1100
f = sample(X0, X1, Y0, Y1, NX, NY)
ext = (X0, X1, Y0, Y1)
poly = np.array([[0, 0], [chord/2, half_t], [chord, 0], [chord/2, -half_t]])

def theory_lines(ax, color='white'):
    kw = dict(color=color, lw=0.9, ls='--', alpha=0.85)
    L = 1.0
    for s in (1, -1):
        ax.plot([0, L*math.cos(b1)], [0, s*L*math.sin(b1)], **kw)                      # LE shock
        for ang in (theta+mu1, -theta+mu2):                                             # expansion fan
            ax.plot([chord/2, chord/2+L*math.cos(ang)], [s*half_t, s*(half_t+L*math.sin(ang))], **kw)
        ang = b3 - theta                                                                # TE shock
        ax.plot([chord, chord+L*math.cos(ang)], [0, s*L*math.sin(ang)], **kw)

def field_plot(data, fname, title, label, cmap, vmin=None, vmax=None, norm=None, lines=True, lc='white'):
    fig, ax = plt.subplots(figsize=(10, 7.2), layout='constrained')
    im = ax.imshow(data, origin='lower', extent=ext, cmap=cmap, vmin=vmin, vmax=vmax,
                   norm=norm, interpolation='bilinear')
    ax.add_patch(Polygon(poly, color='#222222', zorder=10))
    if lines: theory_lines(ax, lc)
    ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1); ax.set_aspect('equal')
    ax.set_xlabel('x [m]'); ax.set_ylabel('y [m]')
    ax.set_title(title + (' — dashed: inviscid theory' if lines else ''))
    fig.colorbar(im, ax=ax, label=label, shrink=0.8)
    fig.savefig(out / fname, dpi=150); plt.close(fig)

field_plot(f['p']/pinf, 'pressure.png', 'Static pressure', r'$p/p_\infty$', 'turbo', 0.55, 1.6)
field_plot(f['Mach'], 'mach.png', 'Mach number', 'Mach', 'turbo', 2.6, 3.4)
uv = f['Umag'].compressed()
field_plot(f['Umag'], 'velocity.png', 'Velocity magnitude', '|U| [m/s]', 'viridis',
           float(np.percentile(uv, 0.5)), float(np.percentile(uv, 99.5)))

# numerical schlieren: |grad rho| on the sampled grid
rho = f['rho'].filled(np.nan)
dy, dx = (Y1-Y0)/(NY-1), (X1-X0)/(NX-1)
gy, gx = np.gradient(rho, dy, dx)
gmag = np.ma.masked_invalid(np.hypot(gx, gy))
gmax = float(np.nanpercentile(gmag.filled(np.nan), 99.9))
schl = np.exp(-12*np.clip(gmag/gmax, 0, 1))       # classic numerical-schlieren mapping
field_plot(schl, 'schlieren.png', 'Numerical schlieren (shock and expansion lines)',
           r'$\exp(-12\,|\nabla\rho|/|\nabla\rho|_{max})$', 'gray', 0, 1, lc='tab:red')

# ---------------------------------------------------------------- surface Cp + forces
xyz = vtk_to_numpy(wall.GetPoints().GetData())
conn = vtk_to_numpy(wall.GetPolys().GetConnectivityArray()).reshape(-1, 4)
v = xyz[conn]; c = v.mean(axis=1)
sf = (np.cross(v[:, 1]-v[:, 0], v[:, 2]-v[:, 0]) + np.cross(v[:, 2]-v[:, 0], v[:, 3]-v[:, 0]))/2
sf[np.sum(sf*(np.array([chord/2, 0, 0.0005])-c), axis=1) < 0] *= -1
pw = vtk_to_numpy(wall.GetCellData().GetArray('p'))
F = np.sum((pw-pinf)[:, None]*sf, axis=0)
span = float(np.ptp(xyz[:, 2]))
Cd, Cl = F[:2] / (q*chord*span)
cp = (pw-pinf)/q
fig, ax = plt.subplots(figsize=(9, 4.5), layout='constrained')
for lab, side, ls in (('Upper surface', c[:, 1] > 0, '-'), ('Lower surface', c[:, 1] < 0, '--')):
    i = np.flatnonzero(side); i = i[np.argsort(c[i, 0])]
    ax.plot(c[i, 0]/chord, cp[i], ls, lw=1.6, label=lab)
cp1, cp2 = (p1-1)*pinf/q, (p2-1)*pinf/q
ax.hlines([cp1, cp2], [0, 0.5], [0.5, 1], colors='k', linestyles=':', lw=1.4, label='Theory')
ax.set(xlabel='x / c', ylabel=r'$C_p$', title=f'Surface pressure coefficient (iteration {int(it)})')
ax.grid(alpha=0.3); ax.legend()
fig.savefig(out / 'surface-cp.png', dpi=150); plt.close(fig)

# ---------------------------------------------------------------- convergence
rows = list(csv.DictReader((root/'postProcessing/steadyDiagnostics/coefficients.csv').open()))
rows = [r for r in rows if float(r['iteration']) > 0]
itn = np.array([float(r['iteration']) for r in rows])
cdh = np.array([float(r['Cd_pressure']) for r in rows])
clh = np.array([float(r['Cl_pressure']) for r in rows])
fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True, layout='constrained')
a1.plot(itn, cdh, lw=1.4, label='CFD'); a1.axhline(cd_th, color='k', ls=':', label=f'Theory {cd_th:.6f}')
a1.set_ylabel(r'$C_d$'); a1.set_ylim(cd_th*0.9, cd_th*1.1); a1.legend(); a1.grid(alpha=0.3)
a2.semilogy(itn, np.abs(clh)+1e-14, lw=1.2, color='tab:orange'); a2.set_ylabel(r'$|C_l|$')
a2.set_xlabel('Iteration (local time stepping)'); a2.grid(alpha=0.3, which='both')
a1.set_title('Convergence history')
fig.savefig(out / 'convergence.png', dpi=150); plt.close(fig)

# ---------------------------------------------------------------- shock angle from schlieren
xs = np.linspace(0.03, 0.09, 40); yy = np.linspace(Y0, Y1, NY)
angles = {}
for s, name in ((1, 'upper'), (-1, 'lower')):
    pts = []
    for x in xs:
        j = int(round((x-X0)/dx)); col = gmag[:, j].filled(0)
        win = np.flatnonzero((s*yy > 0.3*x) & (s*yy < 0.6*x))
        if len(win): pts.append((x, abs(yy[win[np.argmax(col[win])]])))
    pts = np.array(pts); k = np.polyfit(pts[:, 0], pts[:, 1], 1)[0]
    angles[name] = math.degrees(math.atan(k))
res = {
    'iteration': int(it),
    'Cd_pressure': float(Cd), 'Cl_pressure': float(Cl), 'Cd_theory': cd_th,
    'Cd_error_percent': float(100*(Cd-cd_th)/cd_th),
    'shock_angle_deg_cfd': angles, 'shock_angle_deg_theory': math.degrees(b1),
    'p_front_over_pinf_theory': p1, 'p_rear_over_pinf_theory': p2,
    'p_front_over_pinf_cfd': float(rows[-1]['p_front_over_pinf']) if rows else None,
    'p_rear_over_pinf_cfd': float(rows[-1]['p_rear_over_pinf']) if rows else None,
}
(root/'docs'/'results.json').write_text(json.dumps(res, indent=2)+'\n')
print(json.dumps(res, indent=2))
