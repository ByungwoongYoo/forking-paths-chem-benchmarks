import pandas as pd, numpy as np, json
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.SaltRemover import SaltRemover
specs={'ESOL':('delaney-processed.csv','smiles','measured log solubility in mols per litre'),
'FreeSolv':('SAMPL.csv','smiles','expt'),
'Lipophilicity':('Lipophilicity.csv','smiles','exp'),
'BACE':('bace.csv','mol','pIC50'),
'BBBP':('BBBP.csv','smiles','p_np')}
out={}
for name,(f,sc,yc) in specs.items():
    df=pd.read_csv(f)
    rows=[]
    for i,(s,y) in enumerate(zip(df[sc],df[yc])):
        m=Chem.MolFromSmiles(s)
        if m is None: rows.append((i,None,None,None,y,None)); continue
        ik=Chem.MolToInchiKey(m)
        can=Chem.MolToSmiles(m)
        # largest fragment, no stereo
        frags=Chem.GetMolFrags(m,asMols=True)
        lf=max(frags,key=lambda x:x.GetNumHeavyAtoms())
        flat=Chem.MolToSmiles(lf,isomericSmiles=False)
        rows.append((i,can,ik,flat,y,m))
    d=pd.DataFrame(rows,columns=['i','can','ik','flat','y','m'])
    ok=d[d.m.notna()].copy()
    r={'n_rows':len(df),'parse_fail':int(d.m.isna().sum()),'n_parsed':len(ok)}
    for key in ['can','ik','flat']:
        g=ok.groupby(key)
        sz=g.size(); dup=sz[sz>1]
        r[f'{key}_dup_groups']=int(len(dup)); r[f'{key}_rows_in_dup_groups']=int(dup.sum())
        spread=g.y.agg(lambda v:v.max()-v.min())[sz>1]
        r[f'{key}_dup_label_spread_median']=float(spread.median()) if len(spread) else None
        r[f'{key}_dup_label_conflict_groups']=int((spread>0).sum())
        r[f'{key}_dup_label_spread_max']=float(spread.max()) if len(spread) else None
    out[name]=r
    ok[['i','can','ik','flat','y']].to_csv(f'parsed_{name}.csv',index=False)
print(json.dumps(out,indent=1))
json.dump(out,open('dup_stats.json','w'),indent=1)
