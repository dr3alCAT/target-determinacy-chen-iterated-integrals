"""Regression tests for the public API and command-line interface."""

import contextlib
import io
import unittest

import sympy as sp
import mini


class MiniTests(unittest.TestCase):
    def test_weight_one_matrix_keeps_all_relation_rows(self):
        data = mini.build(1)
        self.assertEqual(len(data["rows"]), 4)
        self.assertEqual(data["A"].shape, (4, 0))
        self.assertEqual(data["b"].shape, (4, 1))
        self.assertEqual(data["rank"], data["augmented_rank"])

    def test_build_weight_two_without_prior_build_call(self):
        data = mini.build(2)
        self.assertEqual((len(data["words"]), len(data["cols"]), len(data["rows"])),
                         (16, 12, 32))
        self.assertEqual(data["A"].shape, (32, 12))
        self.assertEqual(data["rank"], 12)
        self.assertEqual(data["augmented_rank"], 12)

    def test_weight_two_exact_value(self):
        mini.build(2)
        expected = (sp.Rational(3, 8) * mini.L**2
                    + sp.I * mini.L * mini.pi / 8 - mini.pi**2 / 96)
        self.assertEqual(sp.simplify(mini.VAL[(1, 2)] - expected), 0)

    def test_weight_three_rank_defect_is_explicit(self):
        data = mini.build(3)
        self.assertEqual((len(data["words"]), len(data["cols"]), len(data["rows"])),
                         (80, 76, 188))
        self.assertEqual(data["A"].shape, (188, 76))
        self.assertEqual(data["rank"], 75)
        self.assertEqual(data["augmented_rank"], 75)
        self.assertEqual(len(data["cols"]) - data["rank"], 1)

    def test_each_build_resets_and_reconstructs_lower_weights(self):
        mini.build(3)
        data = mini.build(2)
        self.assertEqual(data["rank"], 12)
        self.assertIn((1, 2), mini.VAL)

    def test_invalid_build_weights_fail_clearly(self):
        for weight in (0, 4, -1):
            with self.subTest(weight=weight):
                with self.assertRaisesRegex(ValueError, "supported weights"):
                    mini.build(weight)
        with self.assertRaises(TypeError):
            mini.build(2.0)

    def test_matrix_printer_prints_labels_and_augmented_matrix(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            mini.print_matrix(mini.build(2))
        rendered = output.getvalue()
        self.assertIn("Column labels:", rendered)
        self.assertIn("Row labels:", rendered)
        self.assertIn("Exact augmented matrix [A | b]:", rendered)
        self.assertNotIn("\\nColumn labels:", rendered)

    def test_cli_matrix_option_and_range_validation(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertEqual(mini.main(["2", "--matrix"]), 0)
        self.assertIn("Exact augmented matrix [A | b]", output.getvalue())
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as raised:
                mini.main(["4"])
        self.assertEqual(raised.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
