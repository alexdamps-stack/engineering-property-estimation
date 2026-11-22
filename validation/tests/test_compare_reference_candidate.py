"""Check rejection of nonfinite numeric values in reference reproduction."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


COMPARATOR = Path(__file__).resolve().parents[1] / "compare-reference-candidate.py"


class CompareReferenceCandidate(unittest.TestCase):
    def test_nonfinite_values_fail_with_field_path(self):
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / "candidate.json"
            committed = Path(directory) / "committed.json"
            for fresh, expected in [("NaN", "1.0"), ("1.0", "NaN"), ("1e400", "1e400")]:
                with self.subTest(candidate=fresh, committed=expected):
                    candidate.write_text(f'{{"values": [{fresh}]}}', encoding="utf-8")
                    committed.write_text(f'{{"values": [{expected}]}}', encoding="utf-8")
                    result = subprocess.run(
                        [sys.executable, str(COMPARATOR), str(candidate), str(committed)],
                        capture_output=True, text=True,
                    )
                    self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                    self.assertIn("$.values[0]: expected finite numbers", result.stderr)
