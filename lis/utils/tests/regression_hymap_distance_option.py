#!/usr/bin/env python3
"""Compile actual routing init modules with synthetic ESMF/NetCDF API stubs.

No operational maps, MPI/ESMF/NetCDF libraries or LIS jobs are used. This tests
the actual initializer's control flow and typed calls, not the nf90 binding.
Pass a pre-option source tree via --baseline-root for bitwise-default checks.
"""
from pathlib import Path
import argparse,json,os,shutil,subprocess,tempfile

p=argparse.ArgumentParser()
p.add_argument('--source-root',type=Path,default=Path(__file__).resolve().parents[3])
p.add_argument('--baseline-root',type=Path,required=True)
p.add_argument('--compiler',default=shutil.which('gfortran'),required=not bool(shutil.which('gfortran')))
p.add_argument('--sysroot')
p.add_argument('--result-json',type=Path)
a=p.parse_args()
template=Path(__file__).with_name('hymap_distance_init_stubs.F90').read_text()
env=dict(os.environ)
env['PATH']=os.pathsep.join([str(Path(a.compiler).resolve().parent),'/usr/bin','/bin',env.get('PATH','')])
flags=['-isysroot',a.sysroot] if a.sysroot else []
success=[0,1,2,11,12,13,23,24]
failure=[3,4,5,6,7,8,9,10,14,15,16,17,18,19,20,21,22,25,26]
results={'scope':'Actual init modules with typed synthetic API stubs; not actual ESMF/NetCDF/MPI or model execution',
         'routines':{},'compiler':Path(a.compiler).name}
with tempfile.TemporaryDirectory(prefix='hymap_distance_init_') as tmp:
    build=Path(tmp)
    (build/'LIS_misc.h').write_text('! Synthetic compile: runtime header intentionally empty.\n')
    (build/'synthetic_localization_map.nc').write_bytes(b'Synthetic stub fixture only\n')
    for router in ['HYMAP2','HYMAP3']:
        rel=Path(f'lis/routing/{router}_router/da_SWOT/{router}_daSWOT_Mod.F90')
        baseline=(a.baseline_root/rel).read_text().replace(f'module {router}_daSWOT_Mod',
                                                         f'module {router}_daSWOT_Mod_baseline')
        candidate=(a.source_root/rel).read_text()
        assert 'network_distance_km' in candidate
        assert 'allocatable :: distance(:,:)' in candidate
        source=template.replace('HYMAP3',router).replace('! INSERT_MODULES_HERE',baseline+'\n'+candidate)
        source=source.replace('SLICE_X_START','1' if router=='HYMAP2' else '3')
        source=source.replace('SLICE_Y_START','1' if router=='HYMAP2' else '2')
        name=f'{router}_init_regression.F90'
        (build/name).write_text(source)
        entry={'success':[],'expected_errors':[]}
        for netcdf in [True,False]:
            exe=f'{router}_init_'+('netcdf_stub' if netcdf else 'no_netcdf')
            macros=['-DUSE_NETCDF4'] if netcdf else []
            c=subprocess.run([a.compiler]+flags+['-cpp']+macros+
                ['-O0','-g','-fcheck=all','-ffree-line-length-none',name,'-o',exe],
                cwd=build,env=env,text=True,capture_output=True)
            assert c.returncode==0,c.stdout+c.stderr
            cases=success+failure if netcdf else [12,2]
            for case in cases:
                r=subprocess.run([str(build/exe),str(case)],cwd=build,text=True,capture_output=True)
                output=r.stdout+r.stderr
                should_pass=case in success if netcdf else case==12
                if should_pass:
                    assert r.returncode==0 and 'PASS INIT:' in output,(router,case,output)
                    entry['success'].append({'case':case,'netcdf_stub':netcdf,
                                             'result':next(l.strip() for l in r.stdout.splitlines() if 'PASS INIT:' in l)})
                else:
                    assert r.returncode!=0 and ('[ERR]' in output or 'EXPECTED INIT ERROR' in output),(router,case,output)
                    message=next(l.strip() for l in output.splitlines() if '[ERR]' in l or 'EXPECTED INIT ERROR' in l)
                    entry['expected_errors'].append({'case':case,'netcdf_stub':netcdf,'message':message})
        results['routines'][router]=entry
if a.result_json:
    a.result_json.parent.mkdir(parents=True,exist_ok=True)
    a.result_json.write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
