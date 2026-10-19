#!/usr/bin/env python3
"""Compile SNIP adapters with real NoahMP types and minimal LIS/ESMF doubles."""
import configparser
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
MODEL = ROOT / 'lis/surfacemodels/land/noahmp.5.0'


def run(command, cwd):
    return subprocess.run(command, cwd=cwd, check=True, text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout


def main():
    compiler = shlex.split(os.environ.get('FC', 'gfortran'))
    if not shutil.which(compiler[0]):
        raise SystemExit('GNU Fortran is required; set FC to its executable.')
    submodule = MODEL / 'phys/noahmp'
    if not (submodule / 'utility/Machine.F90').is_file():
        raise SystemExit('First run: git submodule update --init '
                         'lis/surfacemodels/land/noahmp.5.0/phys/noahmp')
    with tempfile.TemporaryDirectory(prefix='noahmp50-snip-') as directory:
        work = Path(directory)
        config = configparser.ConfigParser()
        config.read(ROOT / 'lis/make/default.cfg')
        shutil.copy(ROOT / 'lis/make/default.cfg', work)
        shutil.copy(ROOT / 'lis/make/plugins.py', work)
        plugin = ROOT / 'lis/plugins/LIS_lsmda_pluginMod.F90'
        # Exercise the actual build-path generator and CPP guards.
        for model, snow, da in [('NoahMP.5.0', True, True),
                                ('NoahMP.5.0', False, True),
                                ('NoahMP.4.0.1', True, True),
                                ('NoahMP.5.0', True, False)]:
            enabled = {model}
            if snow:
                enabled.add('DA OBS SNIP')
            if da:
                enabled.add('EnKF')
            (work / 'user.cfg').write_text(''.join(
                f'{name}: {name in enabled}\n' for name in config.sections()))
            run([sys.executable, str(work / 'plugins.py')], work)
            paths = (work / 'Filepath').read_text()
            output = run(compiler + ['-cpp', '-E', '-P', '-I', str(work), str(plugin)], work)
            active = model == 'NoahMP.5.0' and da
            assert ('noahmp.5.0/da_snip' in paths) == active
            # da_usafsi is available even with its observation reader disabled.
            assert ('noahmp.5.0/da_usafsi' in paths) == active
            count = output.count('trim(LIS_snipobsId)//char(0),noahmp50_')
            assert count == (10 if active and snow else 0), (model, snow, da, count)
        print('PASS: build-path and plugin guards in four configurations.')

        # Keep production type declarations verbatim; omit only framework
        # initialization, which requires a configured full LIS executable.
        declarations = (MODEL / 'NoahMP50_lsmMod.F90').read_text().split('\ncontains', 1)[0]
        declarations = re.sub(r'^\s*public\s*::\s*NoahMP50_ini\s*$', '',
                              declarations, flags=re.I | re.M)
        (work / 'NoahMP50_lsmMod.F90').write_text(declarations + '\nend module\n')
        sources = [HERE / 'stubs.F90', submodule / 'utility/Machine.F90',
                   submodule / 'drivers/lis/LisNoahmpParamType.F90',
                   MODEL / 'NoahMP50_module.F90', work / 'NoahMP50_lsmMod.F90',
                   MODEL / 'NoahMP50_snowphys_updateMod.F90',
                   MODEL / 'da_usafsi/noahmp50_usafsi_update.F90']
        sources += sorted((MODEL / 'da_snip').glob('*.F90'))
        sources.append(HERE / 'check_interfaces.F90')
        print(run(compiler + ['-cpp', '-O0', '-g', '-fcheck=all',
                             '-ffpe-trap=invalid,zero,overflow',
                             '-o', str(work / 'check_interfaces')]
                  + [str(source) for source in sources], work), end='')
        print(run([str(work / 'check_interfaces')], work), end='')


if __name__ == '__main__':
    try:
        main()
    except subprocess.CalledProcessError as error:
        print(error.stdout, file=sys.stderr)
        print(error.stderr, file=sys.stderr)
        raise SystemExit(error.returncode)
