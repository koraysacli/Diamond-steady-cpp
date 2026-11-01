#!/usr/bin/env python3
"""Read-only steady monitoring; no automatic changes to a running solver."""
from pathlib import Path
import csv, json, sys
root=Path(__file__).resolve().parent
path=root/'postProcessing/steadyDiagnostics/coefficients.csv'
if not path.exists():
    print('Katsayi kaydi henuz yok.'); sys.exit(0)
rows=[]
for r in csv.DictReader(path.open()):
    try: v={k:float(x) for k,x in r.items()}
    except (ValueError,TypeError): continue
    if rows and v['iteration']<=rows[-1]['iteration']: rows=[]
    rows.append(v)
if not rows: sys.exit('Tamamlanmis kayit yok.')
last=rows[-1]; window=[r for r in rows if r['iteration']>=last['iteration']-2000]
def spread(key): return max(r[key] for r in window)-min(r[key] for r in window)
enough=last['iteration']>=5000 and len(window)>=100
criteria={
 'enough_history':enough,
 'Cd_stable':spread('Cd_pressure')<1e-5,
 'Cl_symmetric':max(abs(r['Cl_pressure']) for r in window)<1e-5,
 'mass_flux_balanced':max(r['mass_flux_imbalance'] for r in window)<1e-4,
 'pressure_field_stable':all(0<=r['deltaCp_volume_rms']<1e-5 for r in window),
 'front_pressure_stable':spread('p_front_over_pinf')<1e-3,
 'rear_pressure_stable':spread('p_rear_over_pinf')<1e-3,
 'positive_p_T':all(r['p_min']>0 and r['T_min']>0 for r in window),
}
status={'last':last,'window_samples':len(window),'criteria':criteria,'steady_candidate':all(criteria.values()),'note':'Engineering thresholds, not proof of grid independence or experimental validation. Time is pseudo-iteration.'}
if '--json' in sys.argv: print(json.dumps(status,indent=2))
else:
 print(f"Iterasyon: {last['iteration']:.0f} (fiziksel saniye degil)")
 print(f"Cd (basinc): {last['Cd_pressure']:.8f}    Cl: {last['Cl_pressure']:.8f}")
 print(f"Kutle akisi dengesizligi: {last['mass_flux_imbalance']:.3e}")
 print(f"20 adim arasi hacim agirlikli deltaCp RMS: {last['deltaCp_volume_rms']:.3e}")
 print(f"p_min / p_max: {last['p_min']:.1f} / {last['p_max']:.1f} Pa; T_min: {last['T_min']:.2f} K")
 print('Son 2000 iterasyon: '+('STEADY ADAYI; sok/goruntu/ag kontrolu gerekli.' if all(criteria.values()) else 'Yakinsama henuz dogrulanmadi.'))
 for k,v in criteria.items(): print(f"  {'OK' if v else '--'} {k}")
