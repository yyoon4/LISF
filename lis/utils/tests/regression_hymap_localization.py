#!/usr/bin/env python3
"""Synthetic-only regression checks for the actual HYMAP2/HYMAP3 routines.

Without a compiler this reports only source structure and selector semantics.
With --compiler it builds each full original/candidate subroutine against tiny
ESMF/LIS stubs, compares legacy-safe cases, and exercises the previous overflow.
"""
from pathlib import Path
import argparse, atexit, hashlib, json, os, random, re, shutil, subprocess, sys, tempfile

HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--baseline-root',dest='source_root',type=Path,required=True,
                    help='Unmodified LIS checkout for legacy comparison (no data required)')
parser.add_argument('--candidate-root',type=Path,
                    default=next((p for p in (Path.cwd(),*HERE.parents)
                                  if (p/'lis/routing/HYMAP3_router').is_dir()),None),
                    help='Patched LIS checkout; defaults to this repository')
parser.add_argument('--compiler',default=shutil.which('gfortran'))
parser.add_argument('--sysroot')
args=parser.parse_args()
if args.candidate_root is None:
    parser.error('Pass --candidate-root with a LIS checkout')
build=Path(tempfile.mkdtemp(prefix='hymap-localization-regression-'))
atexit.register(shutil.rmtree,build)

def structure_check(text):
    """Check DO/IF block nesting after joining free-form continuations."""
    statements=[]; current=''
    for line in text.splitlines():
        clean=line.split('!',1)[0].strip().lower()
        if not clean or clean.startswith('#'): continue
        current += ' '+clean.lstrip('&').rstrip('&')
        if not clean.endswith('&'):
            statements.append(current.strip()); current=''
    assert not current
    stack=[]
    for s in statements:
        if re.match(r'^if\s*\(.*\)\s*then\s*$',s): stack.append('if')
        elif re.match(r'^else\s*if\s*\(',s) or s.startswith('elseif') or s=='else':
            assert stack and stack[-1]=='if',s
        elif s in ['endif','end if']:
            assert stack and stack.pop()=='if',s
        elif re.match(r'^do(?:\s|$)',s): stack.append('do')
        elif s in ['enddo','end do']:
            assert stack and stack.pop()=='do',s
    assert not stack,stack

def selector_regression():
    rng=random.Random(13579)
    cases=0
    # Include a direct-write event to preserve the legacy missing-site behavior.
    # Earlier winner history stays intact after that direct assignment.
    for size in [0,1,2,4,16,101,441]:
        for trial in range(40):
            history=[]; old=0.; new=0.; has=False; best=0.
            for i in range(size):
                if rng.random()<.1:
                    value=rng.choice([-5.,-1.,0.,1.,5.])
                    old=value; new=value
                    continue
                w=rng.choice([-2.,-.5,0.,.25,.5,1.,1.5])
                increment=rng.choice([-5.,-1.,0.,1.,5.])
                history.append(w)
                if history.index(max(history))==len(history)-1: old=increment*w
                if not has:
                    has=True; best=w; new=increment*w
                elif w>best:
                    best=w; new=increment*w
                assert old==new,(size,trial,i,old,new)
            cases+=1
    return cases

result={'scope':'Synthetic only; no operational inputs or event data',
        'selector_cases':selector_regression(),
        'native_compiler_available':bool(args.compiler),'routines':{}}
template=(HERE/'hymap_localization_regression_driver.F90').read_text()
for router in ['HYMAP2','HYMAP3']:
    rel=Path(f'lis/routing/{router}_router/da_SWOT/{router}_updateSWOT.F90')
    baseline=(args.source_root/rel).read_text()
    candidate=(args.candidate_root/rel).read_text()
    structure_check(baseline); structure_check(candidate)
    assert 'checkWeight' not in candidate and 'nsfcelevIncr_tmp' not in candidate
    assert 'localupdDX*' not in candidate
    assert 'elseif (localweight(c,r)>bestWeight(t1)) then' in candidate
    assert max(map(len,candidate.splitlines()))<=132
    legacy=baseline.replace(f'subroutine {router}_updateSWOT',f'subroutine {router}_updateSWOT_legacy')
    driver=template.replace('HYMAP3',router).replace('! INSERT_ACTUAL_ROUTINES_HERE',legacy+'\n'+candidate)
    name=f'{router}_regression.f90'
    (build/name).write_text(driver)
    entry={'baseline_sha256':hashlib.sha256(baseline.encode()).hexdigest(),
           'candidate_sha256':hashlib.sha256(candidate.encode()).hexdigest(),
           'structure_check':'passed','actual_subroutine_driver':name}
    if args.compiler:
        exe=f'{router}_regression'
        sysroot_flags=['-isysroot',args.sysroot] if args.sysroot else []
        compiler_env=dict(os.environ)
        compiler_env['PATH']=os.pathsep.join([str(Path(args.compiler).resolve().parent),
                                              '/usr/bin','/bin',os.environ.get('PATH','')])
        compile_run=subprocess.run([args.compiler]+sysroot_flags+['-O0','-g','-fcheck=all','-fbacktrace',
                                   '-Wall','-Wextra','-ffree-line-length-none',name,'-o',exe],
                                  cwd=build,env=compiler_env,text=True,capture_output=True)
        if compile_run.returncode:
            print(compile_run.stdout+compile_run.stderr)
            raise SystemExit(compile_run.returncode)
        actual=subprocess.run([str((build/exe).resolve())],cwd=build,text=True,capture_output=True)
        assert actual.returncode==0,actual.stdout+actual.stderr
        assert 'PASS synthetic Fortran fixtures: 9' in actual.stdout,actual.stdout
        legacy_dense=subprocess.run([str((build/exe).resolve()),'--legacy-dense'],
                                    cwd=build,text=True,capture_output=True)
        error=legacy_dense.stderr.lower()
        assert (legacy_dense.returncode!=0 and 'checkweight' in error and 'above upper bound' in error), error
        entry['native_fixtures']=9
        entry['native_stdout']=actual.stdout.strip()
        entry['legacy_dense_expected_failure']=legacy_dense.stderr.splitlines()[:3]
        entry['compiler']=Path(args.compiler).name
    else:
        entry['native_status']='not executed: no Fortran compiler supplied'
    result['routines'][router]=entry

print(json.dumps(result,indent=2))
if not args.compiler:
    print('SKIP: Fortran regression not executed; compiler required',file=sys.stderr)
    sys.exit(77)
