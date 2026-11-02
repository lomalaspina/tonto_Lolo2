"""Synthetic reader regression: precise precedence, labels, skew cells, fallback.

Uses the native import probe; no experimental dataset or density file needed.
"""
from pathlib import Path
import argparse, json, math, os, subprocess, tempfile
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--program', required=True)
parser.add_argument('--basis-directory', required=True)
args = parser.parse_args()
PROGRAM = str(Path(args.program).resolve())
BASIS = str(Path(args.basis_directory).resolve())
workspace = tempfile.TemporaryDirectory(prefix='tonto-precise-adps-')
WORK = Path(workspace.name)
U={'H1':[.021234567891,.032345678912,.043456789123,.001234567891,-.002345678912,.003456789123],
   'O1':[.011111234567,.022222345678,.033333456789,.000300123456,0.,.000500234567]}
TERMS=['11','22','33','12','13','23']
def fixture(skew,precise,rounded,errors=True,bad=False):
    angles=(80,93,103) if skew else (90,90,90)
    text='data_test\n_symmetry_Int_Tables_number 1\n'
    for key,val in zip(['length_a','length_b','length_c','angle_alpha','angle_beta','angle_gamma'],[8,9,10,*angles]):
        text+=f'_cell_{key} {val}\n'
    text+='loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_U_iso_or_equiv\n_atom_site_adp_type\nH1 .2 .3 .4 .032345 Uani\nO1 .4 .5 .6 .022222 Uani\n'
    text+='loop_\n_atom_site_aniso_label\n'+''.join('_atom_site_aniso_U_'+t+'\n' for t in TERMS)
    for label,values in U.items():
        text+=label+' '+' '.join((f'{v:.3f}(1)' if rounded else f'{v:.12f}(100000000)') for v in values)+'\n'
    if precise:
        cols=TERMS[:-1] if bad else TERMS
        text+='loop_\n_atom_site_Cryst_ADP2_U_label\n'
        text+=''.join('_atom_site_Cryst_ADP2_U_'+t+'\n' for t in cols)
        if errors: text+=''.join('_atom_site_Cryst_ADP2_U_'+t+'_esu\n' for t in cols)
        for label in ('O1','H1'): # Deliberately different from standard/atom order.
            text+=label+' '+' '.join(f'{v:.12f}' for v in U[label][:len(cols)])
            if errors: text+=' '+ ' '.join(['.000100']*len(cols))
            text+='\n'
    return text
def run(name,text,export=False):
    work=WORK/name; work.mkdir()
    (work/'input.cif').write_text(text)
    commands='\nput_cif\nwrite_xtal23_xyz_file\n' if export else '\n'
    (work/'stdin').write_text('{\n name= synthetic\n cif= { file_name= input.cif }\n process_cif'+commands+'}\n')
    cmd=[PROGRAM,str(work/'stdin'),BASIS]
    p=subprocess.run(cmd,cwd=work,text=True,capture_output=True,timeout=30,
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
    (work/'console.log').write_text(p.stdout+p.stderr)
    rows={}
    for line in p.stdout.splitlines():
        fields=line.split()
        if len(fields)==13 and fields[0] in U:
            rows[fields[0]]=list(map(float,fields[1:]))
            assert all(math.isfinite(v) for v in rows[fields[0]])
    return p.returncode,rows,p.stdout+p.stderr
results=[]
for skew in (False,True):
    suffix='skew' if skew else 'orthogonal'
    rc,reference,_=run('legacy-'+suffix,fixture(skew,False,False))
    assert rc==0 and len(reference)==2
    for errors in (True,False):
        name=f'precise-{suffix}-errors-{errors}'
        rc,values,_=run(name,fixture(skew,True,True,errors))
        assert rc==0 and values.keys()==reference.keys()
        error=max(abs(a-b) for key in reference for a,b in zip(reference[key][:6],values[key][:6]))
        uncertainty_error=max(abs(a-b) for key in reference
            for a,b in zip(reference[key][6:] if errors else [0.]*6,values[key][6:]))
        results.append(dict(case=name,max_error_au2=error,
            max_uncertainty_error_au2=uncertainty_error,
            passed=error<1e-12 and uncertainty_error<1e-12))
    rc,_,message=run('incomplete-'+suffix,fixture(skew,True,True,bad=True))
    results.append(dict(case='incomplete-'+suffix,returncode=rc,
        passed=rc>0 and 'incomplete precise crystal ADP2 tensor' in message))
    unknown=fixture(skew,True,True).replace('O1 0.011111','X1 0.011111')
    rc,_,message=run('unknown-label-'+suffix,unknown)
    results.append(dict(case='unknown-label-'+suffix,returncode=rc,
        passed=rc>0 and 'unmatching atom_site_Cryst_ADP2_U_label' in message))
    def mixed(text):
        text=text.replace('O1 .4 .5 .6 .022222 Uani','O1 .4 .5 .6 .022222 Uiso')
        return '\n'.join(line for line in text.splitlines()
                         if not line.startswith('O1 0.'))+'\n'
    rc,reference,_=run('mixed-reference-'+suffix,mixed(fixture(skew,False,False)))
    assert rc==0 and len(reference)==2
    rc,values,_=run('mixed-precise-'+suffix,mixed(fixture(skew,True,True)))
    assert rc==0 and values.keys()==reference.keys()
    error=max(abs(a-b) for key in reference for a,b in zip(reference[key],values[key]))
    results.append(dict(case='mixed-isotropic-'+suffix,max_error_au2=error,passed=error<1e-12))
    name='roundtrip-'+suffix
    source=fixture(skew,False,False).replace('H1 .2 .3 .4',
        'H1 .201234567891 .312345678912 .423456789123')
    rc,reference,_=run(name,source,export=True)
    assert rc==0 and len(reference)==2
    output=(WORK/name/'synthetic.fractional.cif1').read_text()
    rc,values,_=run(name+'-readback',output)
    assert rc==0 and values.keys()==reference.keys()
    error=max(abs(a-b) for key in reference for a,b in zip(reference[key],values[key]))
    # Twelve decimal Angstrom^2 fields allow sub-1e-11 au^2 round-trip error,
    # including the non-orthogonal axis conversion; six decimals do not.
    results.append(dict(case=name,max_error_au2=error,passed=error<2e-11))
    xyz=(WORK/name/'synthetic.xyz').read_text().splitlines()
    hydrogen=[list(map(float,line.split()[1:])) for line in xyz[2:] if int(line.split()[0])==1]
    assert len(hydrogen)==1
    error=max(abs(a-b) for a,b in zip(hydrogen[0],[.201234567891,.312345678912,.423456789123]))
    results.append(dict(case='crystal-handoff-'+suffix,max_error_fractional=error,passed=error<2e-12))
(WORK/'results.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
assert all(r['passed'] for r in results), results
workspace.cleanup()
