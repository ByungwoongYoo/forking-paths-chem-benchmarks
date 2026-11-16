import pandas as pd,numpy as np,json
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator as R
cfg={'bit_r2_2048':dict(c=False,ch=False),'count_r2_2048':dict(c=True,ch=False),'bit_r2_2048_chiral':dict(c=False,ch=True),'count_r2_2048_chiral':dict(c=True,ch=True)}
out={}
for n in ['ESOL','FreeSolv','Lipophilicity','BACE']:
  d=pd.read_csv(f'parsed_{n}.csv'); ms=[Chem.MolFromSmiles(s) for s in d.can]; out[n]={}
  for k,v in cfg.items():
    g=R.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=v['ch'])
    keys=[(g.GetCountFingerprintAsNumPy(m) if v['c'] else g.GetFingerprintAsNumPy(m)).tobytes() for m in ms]
    d['k']=keys; gg=d.groupby('k'); sz=gg.size(); 
    multi=gg.filter(lambda x: x.can.nunique()>1)
    grp=multi.groupby('k')
    spread=grp.y.agg(lambda v:v.max()-v.min())
    out[n][k]=dict(groups_distinct_structures=int(len(spread)),rows=int(len(multi)),spread_median=float(spread.median()) if len(spread) else None,spread_max=float(spread.max()) if len(spread) else None, groups_spread_ge1=int((spread>=1).sum()))
print(json.dumps(out,indent=1)); json.dump(out,open('fp_fix_stats.json','w'),indent=1)
