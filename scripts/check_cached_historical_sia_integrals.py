"""Reject changed source/operands/products; restore original bytes and timestamps."""
import argparse
import os
from pathlib import Path

from collinsflow.historical_sia_integrals import verify_integral_job
from collinsflow.provenance import atomic_json, file_hash, read_json

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('request')
parser.add_argument('output')
parser.add_argument('--julia', default='julia')
args = parser.parse_args()
request, output = Path(args.request).resolve(), Path(args.output).resolve()
verify_integral_job(request, output, julia=args.julia)
paths = [request, *[output/name for name in ('predictions.json', 'audit.json', 'invocation.json', 'point-comparison.csv', 'point-comparison.png')]]
times = {p: (p.stat().st_atime_ns, p.stat().st_mtime_ns) for p in paths}
original = {p: p.read_bytes() for p in paths}
controls = ('figure', 'table', 'rehashed_request_geometry', 'rehashed_request_angular', 'rehashed_request_Q2_label',
    'rehashed_window_bounds', 'rehashed_window_measure', 'rehashed_window_lower_provenance', 'rehashed_integrated_as_point',
    'rehashed_P_mean_as_bound', 'rehashed_error', 'rehashed_covariance', 'rehashed_prediction',
    'rehashed_derivative_samples', 'production_scope')
for mutation in controls:
    rejected = False
    try:
        invocation = read_json(output/'invocation.json')
        if mutation in ('figure', 'table'):
            name = 'point-comparison.png' if mutation == 'figure' else 'point-comparison.csv'
            (output/name).write_bytes(b'explicit altered-product control; original restored')
        elif mutation in ('rehashed_prediction', 'rehashed_derivative_samples'):
            report = read_json(output/'predictions.json')
            if mutation == 'rehashed_prediction':
                report['points'][0]['amplitude'] += .001
            else:
                report['representative_derivatives'][0]['finite_differences'][0]['gradient_samples_decimal'][6][0][0] = '1.0'
            atomic_json(output/'predictions.json', report)
            audit = read_json(output/'audit.json')
            audit['report_sha256'] = file_hash(output/'predictions.json')
            atomic_json(output/'audit.json', audit)
            for name in ('predictions.json', 'audit.json'):
                invocation['outputs'][str(output/name)] = file_hash(output/name)
        elif mutation == 'production_scope':
            invocation['scope'] = 'production'
        else:
            data = read_json(request)
            point = data['points'][0]
            if mutation == 'rehashed_request_geometry':
                point['kinematics']['z1'] += .001
            elif mutation == 'rehashed_request_angular':
                point['source_angular']['average_type'] = 'ratio_of_angular_means'
            elif mutation == 'rehashed_request_Q2_label':
                point['legacy_Q2_status'] = 'primary_verified'
            elif mutation == 'rehashed_window_bounds':
                data['momentum_window']['qT_bounds_GeV'][1] = 3.
            elif mutation == 'rehashed_window_measure':
                data['momentum_window']['measure'] = 'dqT; omitted Jacobian'
            elif mutation == 'rehashed_window_lower_provenance':
                data['momentum_window']['lower_limit_source'] = 'experimental_acceptance_verified'
            elif mutation == 'rehashed_integrated_as_point':
                data['model_convention']['kinematic_prescription'] = 'fixed_reported_P_point'
            elif mutation == 'rehashed_P_mean_as_bound':
                point['mapped_P_bounds_GeV'][0] = point['reported_pT_mean_GeV']
            elif mutation == 'rehashed_error':
                point['measured']['statistical_error'] += .001
            else:
                point['measured']['covariance'] = [[1.]]
            atomic_json(request, data)
            invocation['input_hashes'][str(request)] = file_hash(request)
        atomic_json(output/'invocation.json', invocation)
        try:
            verify_integral_job(request, output, julia=args.julia)
        except ValueError:
            rejected = True
    finally:
        for path, raw in original.items():
            path.write_bytes(raw)
            os.utime(path, ns=times[path])
    if not rejected:
        raise ValueError('cached integral diagnostic accepted '+mutation)
    verify_integral_job(request, output, julia=args.julia)
print('Historical SIA integral cache rejected fifteen controls and accepted restored source/numerical artifacts')
