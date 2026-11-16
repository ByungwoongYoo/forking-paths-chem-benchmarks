import pandas as pd,numpy as np,json,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator as R
d=pd.read_csv('parsed_ESOL.csv')
g=R.GetMorganGenerator(radius=2,fpSize=2048)
d['k']=[g.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)).tobytes() for s in d.can]
ref=g.GetFingerprintAsNumPy(Chem.MolFromSmiles('CCCCCCCCO')).tobytes()
a=d[d.k==ref].copy(); a['nC']=[Chem.MolFromSmiles(s).GetNumAtoms()-1 for s in a.can]; a=a.sort_values('nC')
a[['can','nC','y']].to_csv('esol_alcohol_collision.csv',index=False)
fx=json.load(open('fix_r2.json'))['ESOL']
plt.rcParams.update({'font.size':20,'font.family':'DejaVu Sans'})
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(14,9),gridspec_kw={'width_ratios':[1.6,1]})
ax1.bar([f'C{n}' for n in a.nC],a.y,color='#c0392b')
ax1.set_ylabel('measured log solubility (mol/L)')
ax1.set_title('10 ESOL alcohols, ONE fingerprint',fontweight='bold',fontsize=22)
ax1.text(0.03,0.06,'hexanol → octadecanol\nidentical 2048-bit Morgan (r=2)\nlabels span 7.2 log units',transform=ax1.transAxes,fontsize=17,va='bottom')
ax1.tick_params(axis='x',labelsize=16)
ax1.axhline(0,color='k',lw=0.8)
vals={'bit\n(default)':fx['bit']['mean'],'count':fx['count']['mean']}
ax2.bar(list(vals),list(vals.values()),color=['#7f8c8d','#2471a3'])
for i,v in enumerate(vals.values()): ax2.text(i,v+0.01,f'{v:.2f}',ha='center',fontsize=22,fontweight='bold')
ax2.set_ylim(0,1.12); ax2.set_ylabel('test R², random forest')
ax2.set_title('Same data, same RF',fontweight='bold',fontsize=22)
ax2.text(0.5,0.97,'mean of 10 random splits\ncount wins 10 of 10',transform=ax2.transAxes,ha='center',va='top',fontsize=16)
fig.suptitle('Your fingerprint may not see the molecule',fontsize=28,fontweight='bold')
fig.text(0.01,0.01,'MoleculeNet ESOL (delaney-processed.csv, DeepChem S3). RDKit Morgan r=2, 2048 bits; RF 200 trees, no tuning. @woongworks',fontsize=12,color='#555')
fig.tight_layout(rect=[0,0.03,1,0.95])
fig.savefig('fingerprint-blind-spot.png',dpi=100)
fig.savefig('fingerprint-blind-spot.svg')
print(a[['can','nC','y']])
