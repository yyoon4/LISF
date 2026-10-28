#!/usr/bin/env python3
"""Compile and test the production SNIP updater on a configured server.
No packages are installed. Only LIS state containers/logging are stubbed.
The actual Noah-MP Combo/Divide implementation is compiled unchanged.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, default=None)
    args = ap.parse_args()
    here = Path(__file__).resolve().parent
    lis = here.parents[1]
    out = (args.output or Path('snip-unit-' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))).resolve()
    out.mkdir(parents=True, exist_ok=False)
    fc = shlex.split(os.environ.get('FC', 'gfortran'))
    result = {'status': 'not_run', 'compiler': fc, 'tests': [], 'sources': {}}

    def finish(status, code):
        result['status'] = status
        (out / 'unit_results.json').write_text(json.dumps(result, indent=2) + '\n')
        print(json.dumps({'status': status, 'results': str(out)}, indent=2))
        return code

    if not fc or not shutil.which(fc[0]):
        result['error'] = 'Compiler not found. Load the server compiler and set FC; nothing was installed.'
        return finish('compiler_unavailable', 2)
    version = subprocess.run(fc + ['--version'], capture_output=True, text=True)
    banner = version.stdout + version.stderr
    (out / 'compiler.txt').write_text(banner)
    if 'FFLAGS' in os.environ:
        flags = shlex.split(os.environ['FFLAGS'])
    elif 'GNU Fortran' in banner:
        flags = ['-O0', '-g', '-Wall', '-Wextra', '-fcheck=all', '-fbacktrace',
                 '-ffpe-trap=invalid,zero,overflow', '-ffree-line-length-none']
    elif 'Intel' in banner or Path(fc[0]).name in ('ifort', 'ifx'):
        flags = ['-O0', '-g', '-check', 'all', '-traceback', '-fpe0']
    else:
        result['error'] = 'Unknown compiler; set FFLAGS explicitly for a debug build.'
        return finish('compiler_flags_required', 2)
    sources = [here / 'stubs.F90', lis / 'surfacemodels/land/noahmp.5.0/NoahMP50_snowphys_updateMod.F90',
               lis / 'surfacemodels/land/noahmp.5.0/da_snip/noahmp50_snip_update.F90',
               here / 'test_snip_update.F90']
    for p in sources:
        result['sources'][str(p.relative_to(lis))] = hashlib.sha256(p.read_bytes()).hexdigest()
    executable = out / 'test_snip_update'
    command = fc + flags + [str(p) for p in sources] + ['-o', str(executable)]
    result['compile_command'] = command
    with (out / 'compile.log').open('w') as log:
        compiled = subprocess.run(command, cwd=out, stdout=log, stderr=subprocess.STDOUT)
    if compiled.returncode:
        return finish('compile_failed', 1)
    jobs = [('', None), ('negative', 'negative analysis target'),
            ('inconsistent', 'inconsistent depth/SWE targets'), ('nan', 'non-finite snow state or increment')]
    for argument, expected in jobs:
        run = subprocess.run([str(executable)] + ([argument] if argument else []),
                             cwd=out, capture_output=True, text=True)
        output = run.stdout + run.stderr
        name = argument or 'valid_targets'
        (out / (name + '.log')).write_text(output)
        passed = (run.returncode == 0 and 'SNIP_UNIT_PASS' in output) if expected is None else (
            run.returncode != 0 and expected in output)
        result['tests'].append({'name': name, 'returncode': run.returncode, 'passed': passed})
    return finish('passed' if all(t['passed'] for t in result['tests']) else 'tests_failed',
                  0 if all(t['passed'] for t in result['tests']) else 1)


if __name__ == '__main__':
    sys.exit(main())
