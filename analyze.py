#!/usr/bin/env python3
"""Read a complete snapshot. Plot domain, Mach, surface Cp and estimate shock angle."""
import os
os.environ.setdefault('VTK_SMP_MAX_THREADS','1')
os.environ.setdefault('MPLCONFIGDIR','/tmp/diamond-steady-mpl')
from pathlib import Path
import argparse,json,csv,time,subprocess,sys
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy as array,numpy_to_vtk
vtk.vtkLogger.SetStderrVerbosity(vtk.vtkLogger.VERBOSITY_WARNING)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon,Rectangle
root=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--iteration',type=float);args=ap.parse_args()
roots=sorted(root.glob('processor[0-9]*')) or [root]
valid=[]
for base in roots:
    times=set()
    for d in base.iterdir():
        try:t=float(d.name)
        except ValueError:continue
        if t<=0 or not d.is_dir():continue
        complete=True
        for name in ('p','U','T','rho'):
            f=d/name
            if not f.exists() or time.time()-f.stat().st_mtime<5:complete=False;break
            with f.open('rb') as stream:
                stream.seek(max(0,f.stat().st_size-120));tail=stream.read()
            if b'// *****' not in tail:complete=False;break
        if complete:times.add(t)
    valid.append(times)
times=set.intersection(*valid)
if not times:sys.exit('Tamamlanmis kayit bekleniyor.')
t=args.iteration if args.iteration is not None else max(times)
if t not in times:sys.exit('Secilen kayit tum sureclerde tamamlanmis degil.')
params={}
for line in (root/'case-summary.txt').read_text().splitlines():
    if '=' in line:
        k,v=line.split('=',1)
        try:params[k]=float(v)
        except ValueError:pass
pinf=params['p_Pa'];q=0.5*params['rho']*params['U']**2
reader=vtk.vtkPOpenFOAMReader();reader.SetFileName(str(root/'diamond.foam'))
reader.SetCaseType(0 if len(roots)>1 else 1);reader.CreateCellToPointOff();reader.UpdateInformation()
reader.DisableAllPatchArrays()
for name in ('internalMesh','patch/diamond'):reader.SetPatchArrayStatus(name,1)
reader.DisableAllCellArrays()
for name in ('p','U','T','rho'):reader.SetCellArrayStatus(name,1)
reader.UpdateTimeStep(t)
def leaves(b):
    if isinstance(b,vtk.vtkMultiBlockDataSet):
        for i in range(b.GetNumberOfBlocks()):
            name=b.GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME())
            child=b.GetBlock(i)
            if isinstance(child,vtk.vtkMultiBlockDataSet):yield from leaves(child)
            elif child is not None:yield str(name),child
blocks=dict(leaves(reader.GetOutput()));mesh=blocks['internalMesh'];wall=blocks['diamond']
assert mesh.GetNumberOfCells()==1200000,mesh.GetNumberOfCells()
vals={n:array(mesh.GetCellData().GetArray(n)) for n in ('p','T','rho','U')}
assert all(np.isfinite(v).all() for v in vals.values())
assert all((vals[n]>0).all() for n in ('p','T','rho'))
# Cell-based physical derived fields, before interpolation for visualization.
M=np.linalg.norm(vals['U'],axis=1)/np.sqrt(1.4*287.05287*vals['T'])
for name,v in [('Mach',M),('Cp',(vals['p']-pinf)/q),('p0',vals['p']*(1+0.2*M*M)**3.5)]:
    a=numpy_to_vtk(v,deep=True);a.SetName(name);mesh.GetCellData().AddArray(a)
c2p=vtk.vtkCellDataToPointData();c2p.SetInputData(mesh);c2p.Update()
def sample(bounds,nx,ny):
    x0,x1,y0,y1=bounds
    image=vtk.vtkImageData();image.SetDimensions(nx,ny,1);image.SetOrigin(x0,y0,.0005);image.SetSpacing((x1-x0)/(nx-1),(y1-y0)/(ny-1),1)
    probe=vtk.vtkProbeFilter();probe.SetInputData(image);probe.SetSourceData(c2p.GetOutput());probe.Update()
    out=vtk.vtkImageData();out.DeepCopy(probe.GetOutput())
    data={n:array(out.GetPointData().GetArray(n)).reshape(ny,nx) for n in ('p','Mach','Cp','rho','p0')}
    mask=array(out.GetPointData().GetArray('vtkValidPointMask')).reshape(ny,nx)==0
    data={k:np.ma.array(v,mask=mask) for k,v in data.items()}
    return out,data
outdir=root/'visualization'/str(int(t));outdir.mkdir(parents=True,exist_ok=True)
full,fd=sample((-1,4,-2.5,2.5),900,900)
zoom,zd=sample((-.03,.5,-.18,.18),900,600)
for name,data in [('full-domain.vti',full),('near-profile.vti',zoom)]:
    writer=vtk.vtkXMLImageDataWriter();writer.SetFileName(str(outdir/name));writer.SetInputData(data);writer.Write()
poly=np.array([[0,0],[.1,.01],[.2,0],[.1,-.01]])
fig,axs=plt.subplots(1,2,figsize=(14,6),layout='constrained')
im=axs[0].imshow(fd['p'],origin='lower',extent=(-1,4,-2.5,2.5),cmap='turbo',vmin=1600,vmax=4500)
fig.colorbar(im,ax=axs[0],label='Static pressure [Pa]',shrink=.75)
axs[0].add_patch(Rectangle((-.1,-.2),.6,.4,fill=False,edgecolor='black',ls='--',lw=1,label='Original domain'))
axs[0].legend(loc='upper left');axs[0].set_title('Expanded 5 m × 5 m flow domain')
im=axs[1].imshow(zd['Mach'],origin='lower',extent=(-.03,.5,-.18,.18),cmap='turbo',vmin=2.5,vmax=3.5)
fig.colorbar(im,ax=axs[1],label='Mach number',shrink=.75);axs[1].set_title('Near profile: shocks and expansions')
for ax in axs:
    ax.add_patch(Polygon(poly,color='#333333',zorder=10));ax.set_xlabel('x [m]');ax.set_ylabel('y [m]');ax.set_aspect('equal')
fig.suptitle(f'Mach 3 · α = 0° · 1,200,000 cells · pseudo-iteration {int(t)}\nSnapshot; convergence must be checked separately')
fig.savefig(outdir/'flow.png',dpi=170);plt.close(fig)
# Independent wall-force integration from saved surface fields.
xyz=array(wall.GetPoints().GetData()); polys=wall.GetPolys(); offsets=array(polys.GetOffsetsArray());conn=array(polys.GetConnectivityArray())
assert np.all(np.diff(offsets)==4)
v=xyz[conn.reshape(-1,4)];centres=v.mean(axis=1)
sf=(np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0])+np.cross(v[:,2]-v[:,0],v[:,3]-v[:,0]))/2
sf[np.sum(sf*(np.array([.1,0,.0005])-centres),axis=1)<0]*=-1
pwall=array(wall.GetCellData().GetArray('p'));force=np.sum((pwall-pinf)[:,None]*sf,axis=0)
cd,cl=force[:2]/(q*.0002)
cp=(pwall-pinf)/q
fig,ax=plt.subplots(figsize=(9,4),layout='constrained')
for label,side in [('Upper',centres[:,1]>0),('Lower',centres[:,1]<0)]:
    inds=np.flatnonzero(side);inds=inds[np.argsort(centres[inds,0])]
    ax.plot(centres[inds,0]/.2,cp[inds],label=label,ls='-' if label=='Upper' else '--')
ax.set(xlabel='x / chord',ylabel='Cp',title=f'Surface pressure — pseudo-iteration {int(t)}');ax.grid(alpha=.3);ax.legend()
fig.savefig(outdir/'surface-cp.png',dpi=160);plt.close(fig)
np.savetxt(outdir/'surface-cp.csv',np.column_stack((centres[:,:2],cp)),delimiter=',',header='x,y,Cp',comments='')
_,shock=sample((.04,.18,-.12,.12),150,900)
xx=np.linspace(.04,.18,150);yy=np.linspace(-.12,.12,900)
grad=np.gradient(shock['p'].filled(np.nan),yy,axis=0)
angles={}
for sign,name in [(1,'upper'),(-1,'lower')]:
    ys=[]
    for j,x in enumerate(xx):
        window=np.flatnonzero((sign*yy>.25*x)&(sign*yy<.6*x))
        indices=window[np.isfinite(grad[window,j])]
        if len(indices):ys.append((x,yy[indices[np.argmax(-sign*grad[indices,j])]]))
    data=np.array(ys);fit=np.polyfit(data[:,0],data[:,1],1)
    angles[name]={'beta_deg':float(np.degrees(np.arctan(abs(fit[0])))),'intercept_m':float(fit[1]),'fit_rms_m':float(np.sqrt(np.mean((data[:,1]-np.polyval(fit,data[:,0]))**2)))}
status=json.loads(subprocess.check_output([sys.executable,str(root/'monitor.py'),'--json']))
rows=list(csv.DictReader((root/'postProcessing/steadyDiagnostics/coefficients.csv').open()))
row=next((r for r in reversed(rows) if float(r['iteration'])==t),None)
error=None if row is None else max(abs(cd-float(row['Cd_pressure'])),abs(cl-float(row['Cl_pressure'])))
assert error is None or error<1e-6,('CSV and saved wall forces disagree',error)
summary={'iteration':t,'cell_count':mesh.GetNumberOfCells(),'snapshot_Cd_pressure':float(cd),'snapshot_Cl_pressure':float(cl),'independent_force_error':error,'shock_estimate':angles,'latest_monitor_iteration':status['last']['iteration'],'latest_steady_candidate':status['steady_candidate'],'caution':'Gradient-fit shock estimate; sample/window sensitivity not tested. Visualization is interpolated. No grid independence claim.'}
(outdir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2));print(outdir/'flow.png')
