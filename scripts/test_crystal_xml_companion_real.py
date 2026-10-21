#!/usr/bin/env python3
"""Guard exact phase-wrap equivalence without tolerating changed reflections."""
import unittest

from check_crystal_xml_companion_real import equivalent_phase_rows


HEADER = (b'data_test\nloop_\n_refln_index_h\n_refln_index_k\n_refln_index_l\n'
          b'_refln_F_calc\n_refln_phase_calc\n')
LEFT = HEADER + b'0 0 -1 3.591 180.000\n'


class PhaseRepresentationTests(unittest.TestCase):
    def test_exact_whole_turns(self):
        for phase in (b'-180.000', b'540.000', b'180.0'):
            with self.subTest(phase=phase):
                self.assertEqual(equivalent_phase_rows(
                    LEFT, HEADER + b'0 0 -1 3.591 ' + phase + b'\n'), 1)

    def test_zero_and_full_turn(self):
        self.assertEqual(equivalent_phase_rows(
            HEADER + b'1 2 3 4.0 0.000\n', HEADER + b'1 2 3 4.0 360.000\n'), 1)

    def test_not_an_amplitude_or_hkl_tolerance(self):
        for row in (b'0 0 -1 3.592 -180.000\n', b'0 0 -2 3.591 -180.000\n'):
            with self.subTest(row=row), self.assertRaises(ValueError):
                equivalent_phase_rows(LEFT, HEADER + row)

    def test_not_a_phase_tolerance(self):
        for phase in (b'-179.999', b'179.999999999999', b'-179.' + b'9' * 80,
                      b'0', b'NaN', b'Infinity', b'?'):
            with self.subTest(phase=phase), self.assertRaises(ValueError):
                equivalent_phase_rows(LEFT, HEADER + b'0 0 -1 3.591 ' + phase + b'\n')

    def test_changed_metadata_rejected(self):
        with self.assertRaises(ValueError):
            equivalent_phase_rows(LEFT, LEFT.replace(b'data_test', b'data_other'))

    def test_missing_or_extra_row_rejected(self):
        for text in (HEADER, LEFT + b'0 0 -2 3.591 -180.000\n'):
            with self.subTest(text=text), self.assertRaises(ValueError):
                equivalent_phase_rows(LEFT, text)

    def test_untagged_columns_rejected(self):
        with self.assertRaises(ValueError):
            equivalent_phase_rows(b'0 0 -1 3.591 180.000\n', b'0 0 -1 3.591 -180.000\n')

    def test_reordered_phase_column(self):
        header = (b'loop_\n_refln_phase_calc\n_refln_index_h\n_refln_index_k\n'
                  b'_refln_index_l\n_refln_F_calc\n')
        self.assertEqual(equivalent_phase_rows(
            header + b'180.000 0 0 -1 3.591\n', header + b'-180.000 0 0 -1 3.591\n'), 1)

    def test_new_loop_resets_phase_column(self):
        a = LEFT + b'loop_\n_x\n_y\n_z\n_a\n_b\n0 0 -1 3.591 180.000\n'
        b = LEFT + b'loop_\n_x\n_y\n_z\n_a\n_b\n0 0 -1 3.591 -180.000\n'
        with self.assertRaises(ValueError):
            equivalent_phase_rows(a, b)

    def test_duplicate_phase_tags_rejected(self):
        header = HEADER + b'_refln_phase_calc\n'
        with self.assertRaises(ValueError):
            equivalent_phase_rows(header + b'0 0 -1 3.591 180.000 180.000\n',
                                  header + b'0 0 -1 3.591 180.000 -180.000\n')

    def test_data_boundary_ends_reflection_loop(self):
        for boundary in (b'data_other\n', b'DATA_other\n'):
            with self.subTest(boundary=boundary), self.assertRaises(ValueError):
                equivalent_phase_rows(
                    LEFT + boundary + b'1 2 3 4.0 180.000\n',
                    LEFT + boundary + b'1 2 3 4.0 -180.000\n')

    def test_save_boundary_ends_reflection_loop(self):
        for boundary in (b'save_other\n', b'save_\n', b'SAVE_other\n', b'SAVE_\n'):
            with self.subTest(boundary=boundary), self.assertRaises(ValueError):
                equivalent_phase_rows(
                    LEFT + boundary + b'1 2 3 4.0 180.000\n',
                    LEFT + boundary + b'1 2 3 4.0 -180.000\n')

    def test_stop_boundary_ends_reflection_loop(self):
        for boundary in (b'stop_\n', b'STOP_\n'):
            with self.subTest(boundary=boundary), self.assertRaises(ValueError):
                equivalent_phase_rows(
                    LEFT + boundary + b'1 2 3 4.0 180.000\n',
                    LEFT + boundary + b'1 2 3 4.0 -180.000\n')

    def test_valid_unchanged_wrappers(self):
        loop = HEADER.split(b'\n', 1)[1]
        wrappers = (
            (b'data_test\n', b''),
            (b'data_test\n', b'stop_\n'),
            (b'data_test\nsave_reflections\n', b'save_\n'),
            (b'data_test\nsave_reflections\n', b'stop_\nsave_\n'),
            (b'data_test\n', b'data_next\n_other_value 1\n'),
        )
        for prefix, suffix in wrappers:
            with self.subTest(prefix=prefix, suffix=suffix):
                self.assertEqual(equivalent_phase_rows(
                    prefix + loop + b'0 0 -1 3.591 180.000\n' + suffix,
                    prefix + loop + b'0 0 -1 3.591 -180.000\n' + suffix), 1)

    def test_opposite_phases_rejected(self):
        with self.assertRaises(ValueError):
            equivalent_phase_rows(HEADER + b'0 0 -1 3.591 90.000\n',
                                  HEADER + b'0 0 -1 3.591 -90.000\n')

    def test_reordered_reflections_rejected(self):
        with self.assertRaises(ValueError):
            equivalent_phase_rows(LEFT + b'1 2 3 4.0 0.000\n',
                                  HEADER + b'1 2 3 4.0 0.000\n0 0 -1 3.591 -180.000\n')


if __name__ == '__main__':
    unittest.main()
