#!/usr/bin/env python3
"""Summarize temporary reduction samples; this does NOT certify a completed LIS run."""
import argparse
import json
import math
from pathlib import Path
import sys


def summarize(paths):
    result = {'status': 'no_reduction_samples', 'files': len(paths), 'samples': 0,
              'post_physics_samples': 0, 'max_sample_mapping_depth_error_m': 0.,
              'max_sample_mapping_SWE_error_kgm2': 0.,
              'max_sample_post_physics_depth_difference_m': 0.,
              'max_sample_post_physics_SWE_difference_kgm2': 0., 'warnings': []}
    for path in paths:
        with path.open(errors='replace') as stream:
            for number, line in enumerate(stream, 1):
                marker = next((m for m in ('SNIP_REDUCTION_MAP_WARN', 'SNIP_REDUCTION_SKIPPED',
                                          'SNIP_REDUCTION_MAP', 'SNIP_REDUCTION_POST') if m+' ' in line), None)
                if marker is None:
                    if 'SNIP_TARGET_CHECK ' in line or 'SNIP_TRACE_SNOW_QC ' in line:
                        result['warnings'].append({'file': str(path), 'line': number, 'message': 'obsolete diagnostics: check binary/log provenance'})
                    elif '[ERR]' in line:
                        result['warnings'].append({'file': str(path), 'line': number, 'message': line.strip()})
                    continue
                if marker in ('SNIP_REDUCTION_MAP_WARN', 'SNIP_REDUCTION_SKIPPED'):
                    result['warnings'].append({'file': str(path), 'line': number, 'message': line.strip()})
                    continue
                try:
                    fields = line.split(marker, 1)[1].split()
                    if len(fields) != 7:
                        raise ValueError('expected nest, local call, tile and four real values')
                    tuple(map(int, fields[:3]))
                    values = [float(v.replace('D', 'E').replace('d', 'e')) for v in fields[3:]]
                    if not all(math.isfinite(v) for v in values):
                        raise ValueError('non-finite diagnostic')
                    if marker == 'SNIP_REDUCTION_MAP':
                        result['samples'] += 1
                        h,s = values[2:]
                        if h < 0 or s < 0:
                            raise ValueError('negative absolute mapping error')
                        result['max_sample_mapping_depth_error_m'] = max(result['max_sample_mapping_depth_error_m'], h)
                        result['max_sample_mapping_SWE_error_kgm2'] = max(result['max_sample_mapping_SWE_error_kgm2'], s)
                    else:
                        result['post_physics_samples'] += 1
                        h,s = abs(values[2]-values[0]),abs(values[3]-values[1])
                        result['max_sample_post_physics_depth_difference_m'] = max(result['max_sample_post_physics_depth_difference_m'],h)
                        result['max_sample_post_physics_SWE_difference_kgm2'] = max(result['max_sample_post_physics_SWE_difference_kgm2'],s)
                except ValueError as exc:
                    result['warnings'].append({'file': str(path), 'line': number, 'message': str(exc)})
    if result['samples'] != result['post_physics_samples']:
        result['warnings'].append({'message': 'unpaired mapping/post-physics sample counts'})
    if result['warnings']:
        result['status'] = 'needs_review'
    elif result['samples']:
        result['status'] = 'samples_recorded_no_reported_warnings'
    result['scope'] = ('At most three successful reduction samples per setter call/rank; all applied reductions '
                       'are checked for mapping residuals in the updater. No global reduction count or proof of '
                       'run completion, full rank coverage, or scientific accuracy. Post-physics differences '
                       'are descriptive, not mapping failures. Retain all rank logs and check completion separately.')
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--logs', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    paths = sorted(p for p in args.logs.rglob('lislog*') if p.is_file())
    if not paths:
        ap.error('No lislog* files found; retain all MPI-rank logs from the same run.')
    result = summarize(paths)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
    return 0 if result['status'] == 'samples_recorded_no_reported_warnings' else 1


if __name__ == '__main__':
    sys.exit(main())
