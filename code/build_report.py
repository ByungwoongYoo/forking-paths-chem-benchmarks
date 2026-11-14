"""Fill the official Apart template (Google Doc export) with the report text. New sprint code."""
import re, copy
import pandas as pd
import docx
from docx.shared import Inches, Pt
from pathlib import Path
ROOT = str(Path(__file__).resolve().parents[1])
TITLE = 'Forking Paths in Molecular Property Benchmarks: How Far a Reported Score Moves Under Defensible Choices, and a Training-Free Audit Signal'
AUTHOR = 'Byungwoong Yoo\nBigWonder (independent), Seoul, South Korea'
ABSTRACT = open(f'{ROOT}/report/abstract.txt').read().strip()

TEMPLATE = Path(ROOT) / 'report' / 'apart_template_original.docx'
if not TEMPLATE.exists():
    raise SystemExit('report/apart_template_original.docx is not in this repo (it is the official Apart sprint template, not ours to redistribute). The finished report is already in report/report.pdf and report/report.md. To regenerate the .docx, place the official template at that path and re-run.')
d = docx.Document(str(TEMPLATE))
hdr = d.tables[0]
def all_cell_paras(t):
    for r in t.rows:
        for c in r.cells:
            for p in c.paragraphs: yield p
            for tt in c.tables: yield from all_cell_paras(tt)
for p in all_cell_paras(hdr):
    for r in p.runs:
        if r.text.strip().startswith('PROJECT TITLE'): r.text = TITLE
        m = re.match(r'Author name (\d)', r.text)
        if m: r.text = AUTHOR if m.group(1) == '1' else ''
    if p.text.startswith('Summarize your project'):
        for r in p.runs[1:]: r.text = ''
        p.runs[0].text = ABSTRACT; p.runs[0].italic = False

from docx.oxml.ns import qn
def widen(t):
    for r in t.rows:
        for c in r.cells:
            for tt in c.tables:
                if len(tt.columns) == 3:
                    for gc in tt._tbl.tblGrid.findall(qn('w:gridCol')): gc.set(qn('w:w'), '3100')
                    for tc in tt._tbl.iter(qn('w:tc')):
                        tcPr = tc.get_or_add_tcPr(); w = tcPr.find(qn('w:tcW'))
                        if w is None: w = tcPr.makeelement(qn('w:tcW'), {}); tcPr.append(w)
                        w.set(qn('w:w'), '3100'); w.set(qn('w:type'), 'dxa')
                    tblW = tt._tbl.tblPr.find(qn('w:tblW'))
                    if tblW is not None: tblW.set(qn('w:w'), '9300'); tblW.set(qn('w:type'), 'dxa')
                widen(tt)
widen(hdr)
# drop everything after the header table except sectPr
body = d.element.body
kids = list(body); start = kids.index(hdr._tbl) + 1
for el in kids[start:]:
    if not el.tag.endswith('sectPr'): body.remove(el)

MD=[f'# {TITLE}', '', AUTHOR.replace(chr(10), ', '), '', '**Abstract.** '+ABSTRACT, '']
def H(t, lvl=2):
    MD.append('#'*lvl+' '+t); MD.append(''); d.add_paragraph(t, style=f'Heading {lvl}')
def P(t):
    MD.append(t); MD.append('')
    p = d.add_paragraph(style='normal')
    # **bold** markup
    for i, part in enumerate(re.split(r'\*\*', t)):
        if part: r = p.add_run(part); r.bold = (i % 2 == 1)
    return p
def FIG(path, cap, w=6.3):
    MD.append(f'![]({path.replace(ROOT, "..")})')
    d.add_picture(path, width=Inches(w)); c = P(cap); c.runs[0].italic = True
def TAB(df, cap, fs=8):
    MD.append(df.to_markdown(index=False)); MD.append('')
    c = P(cap); c.runs[0].italic = True
    t = d.add_table(rows=1, cols=len(df.columns)); t.style = d.styles['Table Grid'] if 'Table Grid' in [s.name for s in d.styles] else None
    for j, col in enumerate(df.columns): t.rows[0].cells[j].text = str(col)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for j, v in enumerate(row): cells[j].text = str(v)
    for row in t.rows:
        for cell in row.cells:
            for p in cell.paragraphs:
                for r in p.runs: r.font.size = Pt(fs)
    d.add_paragraph(style='normal')

exec(open(f'{ROOT}/report/body.py').read())
d.save(f'{ROOT}/report/report.docx')
open(f'{ROOT}/report/report.md','w').write(chr(10).join(MD))
print('saved')
