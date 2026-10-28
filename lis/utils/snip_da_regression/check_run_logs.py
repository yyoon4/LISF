#!/usr/bin/env python3
"""Summarize target-conservation checks from ALL MPI rank logs of a new LIS run."""
import argparse
import json
import math
from pathlib import Path
import sys


def summarize(paths):
    result = {'status': 'not_checked', 'files': len(paths), 'files_with_checks': 0, 'checks': 0,
              'patch_updates': 0, 'changed': 0, 'reduced': 0, 'layers_changed': 0, 'failed': 0,
              'max_depth_error_m': 0., 'max_SWE_error_kgm2': 0.,
              'max_layer_mass_error_kgm2': 0., 'max_layer_depth_error_m': 0., 'errors': []}
    for p in paths:
        seen = False
        with p.open(errors='replace') as stream:
            for line_number, line in enumerate(stream, 1):
                if '[ERR]' in line and ('SNIP_UPDATE' in line or 'SNIP_TARGET_FAIL' in line):
                    result['errors'].append({'file': str(p), 'line': line_number, 'message': line.strip()})
                if 'SNIP_TARGET_CHECK ' not in line:
                    continue
                try:
                    fields = line.split('SNIP_TARGET_CHECK ', 1)[1].split()
                    if len(fields) != 11:
                        raise ValueError('expected 7 integers and 4 reals')
                    values = list(map(int, fields[:7]))
                    reals = [float(x.replace('D', 'E').replace('d', 'e')) for x in fields[7:]]
                    if any(not math.isfinite(v) or v < 0 for v in reals) or any(v < 0 for v in values):
                        raise ValueError('non-finite or negative diagnostic')
                    result['checks'] += 1
                    seen = True
                    for key, value in zip(['patch_updates', 'changed', 'reduced', 'layers_changed', 'failed'], values[2:]):
                        result[key] += value
                    for key, value in zip(['max_depth_error_m', 'max_SWE_error_kgm2',
                                           'max_layer_mass_error_kgm2', 'max_layer_depth_error_m'], reals):
                        result[key] = max(result[key], value)
                except ValueError as exc:
                    result['errors'].append({'file': str(p), 'line': line_number, 'message': str(exc)})
        result['files_with_checks'] += int(seen)
    if result['errors'] or result['failed']:
        result['status'] = 'failed'
    elif not result['checks'] or not result['changed'] or not result['reduced']:
        result['status'] = 'incomplete_no_reduction_coverage'
    else:
        result['status'] = 'passed_target_mapping_checks'
    result['scope'] = ('Per-call mapping checks only; not proof of completed simulation, '
                       'correct EnKF analysis, or scientific accuracy. Check normal LIS completion separately.')
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--logs', type=Path, required=True, help='New run logs directory, not a mix of old/new runs')
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    paths = sorted(p for p in args.logs.rglob('lislog*') if p.is_file())
    if not paths:
        ap.error('No lislog* files found; retain all MPI rank logs.')
    result = summarize(paths)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
    return 0 if result['status'] == 'passed_target_mapping_checks' else 1


if __name__ == '__main__':
    sys.exit(main())
