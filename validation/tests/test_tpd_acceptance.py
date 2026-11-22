"""Exercise reference stability evidence without invoking Rust or EOS backends."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np


VALIDATION = Path(__file__).resolve().parents[1]


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, VALIDATION / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TpdAcceptance(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cubic = load("cubic_reference_test", "generate-flash-references.py")
        cls.gerg = load("gerg_reference_test", "generate-gerg-equilibrium-references.py")

    def search(self, module, results):
        backend = module.a4 if module is self.gerg else module
        model = SimpleNamespace(trial_records=[])
        with patch.object(backend, "min_gibbs", return_value=(1.0, 1.0, np.zeros(2))), \
                patch.object(backend, "wilson", return_value=np.ones(2)), \
                patch.object(backend, "RANDOM_STARTS", 0), \
                patch.object(module, "minimize", side_effect=results):
            return module.tpd_search(model, ["methane", "propane"], 300., 1e6, np.array([.5, .5]))

    @staticmethod
    def result(value, success, x=None):
        return SimpleNamespace(fun=value, success=success, message="synthetic stop",
                               x=np.log([.5, .5]) if x is None else x)

    def test_failed_nonnegative_or_roundoff_minima_are_unresolved(self):
        for module in (self.cubic, self.gerg):
            for value in (1., 0., -1e-12, float("nan"), float("inf")):
                with self.subTest(module=module.__name__, value=value):
                    with self.assertRaisesRegex(AssertionError, "TPD unresolved"):
                        self.search(module, [self.result(value, False)] * 2)

    def test_successful_minimum_supports_stability(self):
        for module in (self.cubic, self.gerg):
            with self.subTest(module=module.__name__):
                result = self.search(module, [self.result(-1e-12, False), self.result(0., True)])
                self.assertEqual(result[0], 0.)

    def test_failed_negative_trial_is_an_instability_witness(self):
        for module in (self.cubic, self.gerg):
            with self.subTest(module=module.__name__):
                result = self.search(module, [self.result(-.01, False), self.result(0., True)])
                self.assertEqual(result[0], -.01)

    def test_nonfinite_coordinates_cannot_support_stability(self):
        for module in (self.cubic, self.gerg):
            with self.subTest(module=module.__name__):
                with self.assertRaisesRegex(AssertionError, "TPD unresolved"):
                    self.search(module, [self.result(0., True, np.array([np.nan, 0.]))] * 2)


if __name__ == "__main__":
    unittest.main()
