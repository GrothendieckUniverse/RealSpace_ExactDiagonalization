"""Exercise generated job resources, flux protocol, and submission filtering.

Run: python3 test/hpc_workflow.py
Uses a temporary checkout and mocked Slurm commands; never submits real jobs.
"""
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class HPCWorkflow(unittest.TestCase):
    def test_generated_flux_campaign(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hpc = root / 'phase_exploration/hpc'
            hpc.mkdir(parents=True)
            generator = hpc / 'hyak_slurm_gen.sh'
            shutil.copy(ROOT / 'phase_exploration/hpc/hyak_slurm_gen.sh', generator)
            env = dict(os.environ, REPO_DIR=str(root), FLOW_STEPS='17', PUMP_STEPS='17',
                       FLOW_CYCLES='1', PUMP_CYCLES='1', DIAGNOSTIC_OBSERVABLES='flow,pump',
                       FLUX_DIRECTION='1', POLARIZATION_DIRECTION='2')
            subprocess.run(['bash', str(generator)], env=env, check=True, capture_output=True)
            generated = hpc / 'generated'
            jobs = list(generated.glob('diagnostics_*.sbatch'))
            self.assertEqual(len(jobs), 33)
            for job in generated.glob('*.sbatch'):
                subprocess.run(['bash', '-n', str(job)], check=True)
            for job in jobs:
                text = job.read_text()
                self.assertIn('--flow-steps "17"', text)
                self.assertIn('--pump-steps "17"', text)
                self.assertIn('--flow-cycles "1"', text)
                self.assertIn('--pump-cycles "1"', text)
                self.assertIn('--observables "flow,pump"', text)
                self.assertIn('--refresh true', text)
                self.assertIn('diagnostics_v11_', text)
                self.assertIn('#SBATCH --ntasks=1\n', text)
                cpus = re.search(r'--cpus-per-task=(\d+)', text)[1]
                self.assertIn(f'export JULIA_NUM_THREADS={cpus}\n', text)
                required = re.findall(r'^# PHASE_STUDY_REQUIRED_OUTPUT=(.*)$', text, re.M)
                self.assertFalse(any('structure_' in Path(r).name or 'entanglement' in Path(r).name for r in required))
            # Old completion markers cannot satisfy the new protocol.
            for job in jobs:
                required = re.findall(r'^# PHASE_STUDY_REQUIRED_OUTPUT=(.*)$', job.read_text(), re.M)
                for name in required:
                    path = Path(name)
                    path.parent.mkdir(parents=True, exist_ok=True)
                    if 'v11_' not in path.name:
                        path.write_text('legacy data\n')
                    else:
                        path.with_name('characteristic_points_v9.done').write_text('old protocol\n')
            mockbin = root / 'mockbin'
            mockbin.mkdir()
            (mockbin / 'squeue').write_text('#!/bin/sh\nexit 0\n')
            (mockbin / 'sbatch').write_text('''#!/usr/bin/env python3
import json,os,sys
from pathlib import Path
p=Path(os.environ['MOCK_SUBMISSIONS'])
rows=json.loads(p.read_text()) if p.exists() else []
rows.append(sys.argv[1:]);p.write_text(json.dumps(rows))
print(1000+len(rows))
''')
            for file in mockbin.iterdir():
                file.chmod(0o755)
            records = root / 'submissions.json'
            env.update(PATH=str(mockbin) + os.pathsep + os.environ['PATH'],
                       MOCK_SUBMISSIONS=str(records), USER='test')
            subprocess.run(['bash', str(generated / 'submit_all.sh'), '--kind', 'diagnostics'],
                           env=env, check=True, capture_output=True)
            submissions = json.loads(records.read_text())
            data = [r for r in submissions if Path(r[-1]).name.startswith('diagnostics_')]
            self.assertEqual(len(data), 33)
            self.assertEqual(len(submissions), 35)  # setup, 33 data, plot
            self.assertTrue(all('--dependency=afterok:1001' in r for r in data))
            self.assertTrue(Path(submissions[-1][-1]).name == 'plot_results.sbatch')
            self.assertTrue(any(x.startswith('--dependency=afterok:') for x in submissions[-1]))
            # The larger all-observables campaign remains available.
            env['DIAGNOSTIC_OBSERVABLES'] = 'all'
            subprocess.run(['bash', str(generator)], env=env, check=True, capture_output=True)
            for job in generated.glob('diagnostics_*.sbatch'):
                text = job.read_text()
                if '_fci_' in job.name:
                    self.assertIn('--observables "structure,flow,pump,pes"', text)
                else:
                    self.assertIn('--observables "structure,flow,pump"', text)


if __name__ == '__main__':
    unittest.main()
