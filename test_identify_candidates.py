"""
Tests for identifying candidates
"""

import os
import unittest
import numpy as np

from zero_points import zero_points
from identify_candidates_v2 import calculate_snr, write_region_file,\
    write_txt_file, write_output, update_candidate_red_list, remove_bad_values,\
    Inputs, MagCutSelection, ClassicSNR
class TestSNR(unittest.TestCase):
    """Testing that SNR function works."""
    def test_calculate_snr(self):
        """Calculating the snr"""
        snr_val = calculate_snr(0.1)
        self.assertEqual(snr_val, (2.5/np.log(10))/0.1)

class TestRegionFile(unittest.TestCase):
    """Testing that the output region file is being created."""
    ra = np.arange(10)
    dec = np.arange(10)

    def test_region_exists(self):
        """Testing that the region file actually exists and isn't empty."""
        file_name = 'reg_test.reg'
        write_region_file(self.ra, self.dec, file_name)
        val = os.system(f'ls {file_name}')
        self.assertEqual(val, 0)
        data = np.loadtxt(file_name, dtype=str)
        self.assertNotEqual(len(data), 0)
        os.system(f'rm {file_name}')

class TestTextfile(unittest.TestCase):
    """Testing that txt output is working."""
    ra = np.arange(10)
    dec = np.arange(10)

    def test_txt(self):
        """Testing file is created and is not empty."""
        file_name = 'txt_test.txt'
        write_txt_file(self.ra, self.dec, file_name)
        val = os.system(f'ls {file_name}')
        self.assertEqual(val, 0)
        data = np.loadtxt(file_name, dtype=str)
        self.assertNotEqual(len(data), 0)
        os.system(f'rm {file_name}')

class TestWriteOutput(unittest.TestCase):
    """Testing that the outputs are being written."""
    ra = np.arange(10)
    dec = np.arange(10)
    def test_exists(self):
        """Testing that the files exist."""
        suffix = 'output_test'
        write_output(self.ra, self.dec, suffix)
        val_reg = os.system(f'ls {suffix}.reg')
        val_txt = os.system(f'ls {suffix}.txt')
        self.assertEqual(val_reg, 0)
        self.assertEqual(val_txt, 0)
        os.system(f'rm {suffix}.reg')
        os.system(f'rm {suffix}.txt')

class TestUpdateRedList(unittest.TestCase):
    """Testing that the updating of the red list works"""
    ra = np.arange(10)
    dec = np.arange(10)
    file_name = 'test_red.txt'
    os.system(f'rm {file_name}')
    def test_file_exists(self):
        """Testing that the file is actually made."""
        update_candidate_red_list(self.ra, self.dec, self.file_name)
        val = os.system(f'ls {self.file_name}')
        self.assertEqual(val, 0)
        os.system(f'rm {self.file_name}')

    def test_append(self):
        """Testing that the append is not skipping lines"""
        update_candidate_red_list(self.ra, self.dec, self.file_name)
        update_candidate_red_list(self.ra, self.dec, self.file_name)
        # Running the command again should duplicate the results
        # Therefore data being 10 + 10 = 20 rows.
        data = np.loadtxt(self.file_name)
        self.assertEqual(len(data), 20)
        os.system(f'rm {self.file_name}')

class TestRemovingRedValues(unittest.TestCase):
    """Testing that removing from the redlist is working."""
    all_ra = np.arange(10)
    all_dec = np.arange(10)
    bad_ra = np.arange(3)
    bad_dec = np.arange(3)
    red_list_name = 'test_remove_red.txt'
    update_candidate_red_list(bad_ra, bad_dec,red_list_name)

    def test_removing_bad(self):
        """Testing that only the good remain."""
        good_ra, good_dec = remove_bad_values(self.all_ra, self.all_dec, self.red_list_name)
        self.assertEqual(list(good_ra), list(np.arange(3,10)))
        self.assertEqual(list(good_dec), list(np.arange(3, 10)))
        os.system(f'rm {self.red_list_name}')

class TestMagCutSelection(unittest.TestCase):
    """Testing that the mag selection criteria is functioning correctly."""
    inputs = Inputs(
        'test_red.txt', 'test_mag', 'test_catalogs/n964.cat', 'test_catalogs/n964_135.cat',
        'test_catalogs/i.cat', 'test_catalogs/z.cat', zero_points, images = ('1','1','1'),
        aperture_radii=1.
    )

    selection = MagCutSelection(inputs, 24.2, 24, 25.8, 25.6)

    def test_n964_data(self):
        """Testing tha reading in the n964 data is working."""
        self.assertAlmostEqual(self.selection.n964_data[0][0], 23.20313014)
        self.assertAlmostEqual(self.selection.n964_data[0][3], 29.00313014)
        self.assertAlmostEqual(self.selection.n964_data[1][0], 23.07263309)
        self.assertAlmostEqual(self.selection.n964_data[1][3], 28.77263309)

    def test_i_data(self):
        """Testing that the i data is converting correctly."""
        self.assertAlmostEqual(self.selection.i_data[0][0], 129.93983036)
        self.assertAlmostEqual(self.selection.i_data[0][3], 20.93983036)

    def test_z_data(self):
        """Testing that the z data is converting correctly."""
        self.assertAlmostEqual(self.selection.z_data[0][0], 25.56102665)
        self.assertAlmostEqual(self.selection.z_data[0][1],129.56102665)
        self.assertAlmostEqual(self.selection.z_data[0][3],-68.43897335)

    def test_n964_selection(self):
        """Testing the n964 selection."""
        good_n964_vals = self.selection.select_n964()
        self.assertEqual(list(good_n964_vals), [0, 1, 2, 10])

    def test_i_selection(self):
        """Testing that the i selection."""
        good_i_vals = self.selection.select_i()
        self.assertEqual(list(good_i_vals), [0, 1, 2, 9])

    def test_z_selection(self):
        """Testing z selection."""
        good_z_vals = self.selection.select_z()
        self.assertEqual(list(good_z_vals), [0, 1, 2, 5])

    def test_apply_selection(self):
        """Testing that the selection is being applied correctly."""
        reduced = self.selection.apply_selection_criteria()
        self.assertEqual(list(reduced), [0, 1, 2])

class TestClassicSNR(unittest.TestCase):
    """Testing that the selection criteria for snr is working correctly."""
    inputs = Inputs(
        'test_red.txt', 'test_mag', 'test_catalogs/n964.cat', 'test_catalogs/n964_135.cat',
        'test_catalogs/i.cat', 'test_catalogs/z.cat', zero_points, images = ('1','1','1'),
        aperture_radii=1.
    )
    selection = ClassicSNR(inputs, 5, 5, 3, 3)

    def test_n964_selection(self):
        """Testing that the selection criteria is working."""
        self.assertEqual(list(self.selection.select_n964()),[0, 1, 2])

    def test_i_selection(self):
        """Testing i selection working correctly"""
        self.assertEqual(list(self.selection.select_i()), [0, 1, 2, 9])

    def test_z_selection(self):
        """Testing that the z selection is working correctly."""
        self.assertEqual(list(self.selection.select_z()), [0, 1, 2, 5])

    def test_apply_selection(self):
        """Testing all combinations."""
        reduced = self.selection.apply_selection_criteria()
        self.assertEqual(list(reduced), [0, 1, 2])

if __name__ == '__main__':
    unittest.main()
