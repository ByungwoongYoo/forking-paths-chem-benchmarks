import pandas as pd,numpy as np,json
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator as R
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
res={}
for n in ['ESOL','BACE']:
  d=pd.read_csv(f'parsed_{n}.csv'); ms=[Chem.MolFromSmiles(s) for s in d.can]; y=d.y.values
  feats={}
  for name,c,ch in [('bit',0,0),('count',1,0),('count_chiral',1,1)]:
    g=R.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=bool(ch))
    feats[name]=np.array([(g.GetCountFingerprintAsNumPy(m) if c else g.GetFingerprintAsNumPy(m)) for m in ms]).astype(np.float32)
  rr={k:[] for k in feats}
  for seed in range(10):
    p=np.random.RandomState(seed).permutation(len(d)); tr=p[:int(.8*len(d))]; te=p[int(.9*len(d)):]
    for k,X in feats.items():
      rf=RandomForestRegressor(n_estimators=200,n_jobs=-1,random_state=seed).fit(X[tr],y[tr]); rr[k].append(r2_score(y[te],rf.predict(X[te])))
  res[n]={k:dict(mean=float(np.mean(v)),sd=float(np.std(v,ddof=1)),wins_vs_bit=int(sum(a>b for a,b in zip(v,rr['bit'])))) for k,v in rr.items()}
  print(n,res[n])
json.dump(res,open('fix_r2.json','w'),indent=1)
