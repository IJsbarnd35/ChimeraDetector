#!/usr/bin/env python3
"""
unit_tests.py: created for the testing of the tool.
Author: IJbrand Pool
Version: 2.1.0
Date: 09-06-2023
"""

import subprocess
import unittest
import pandas as pd
import os
from unittest.mock import patch, call

from Bio.Seq import Seq

from Mapping_Method import create_df, find_chimeras
from Remove_Chimeras import remove_from_file
from Self_Aligned import check_chimera_flat, file_reader, self_aligned
from Create_Artificial_Reads import make_chimera, add_errors, make_normal


class TestCreateDF(unittest.TestCase):
    def setUp(self):
        # Create a sample overlaps.paf file
        paf_data = [
            "read1\t100\t0\t100\t+\tread1\t200\t0\t200",
            "read2\t200\t0\t200\t+\tread2\t300\t0\t300",
            "read3\t150\t0\t150\t+\tread3\t250\t0\t250",
        ]
        with open("overlaps.paf", "w") as f:
            f.write("\n".join(paf_data))

    def tearDown(self):
        # Clean up the sample overlaps.paf file
        os.remove("overlaps.paf")

    @patch('subprocess.run')
    def test_create_df(self, mock_subprocess_run):
        # Set up the mock subprocess.run to return a successful completion
        mock_subprocess_run.return_value.returncode = 0

        try:
            # Call the function to create the DataFrame
            df = create_df("overlaps.paf")

            # Check if the DataFrame is not empty
            self.assertFalse(df.empty)

            # Check if the DataFrame has the correct columns
            expected_columns = ['query_id', 'query_length', 'query_start', 'query_end', 'strand', 'target_id',
                                'target_length', 'target_start', 'target_end']
            self.assertEqual(list(df.columns), expected_columns)

            # Check if the DataFrame has the correct number of rows
            self.assertEqual(len(df), 3)

            # Check if the subprocess.run was called with the correct arguments
            mock_subprocess_run.assert_called_once_with(
                "minimap2 -x ava-ont reads.fa reads.fa > overlaps.paf",
                shell=True,
                check=True
            )

        except FileNotFoundError:
            self.fail("The test failed because the overlaps.paf file was not found.")

        except subprocess.CalledProcessError:
            self.fail("The test failed because an error occurred while executing the 'minimap2' command.")


class TestFindChimeras(unittest.TestCase):
    def test_find_chimeras(self):
        # Create a sample DataFrame for testing
        data = {
            'query_id': ['read=1', 'read=2', 'read=3', 'read=4', 'read=5', 'read=6'],
            'query_length': [1000, 2000, 3000, 4000, 5000, 6000],
            'query_start': [10, 10, 10, 10, 10, 10],
            'query_end': [990, 1990, 2990, 3990, 4990, 5990],
            'strand': ['-', '-', '-', '-', '-', '-'],
            'target_id': ['read=1', 'read=2', 'read=3', 'read=4', 'read=5', 'read=6'],
            'target_length': [1000, 2000, 3000, 4000, 5000, 6000],
            'target_start': [10, 10, 10, 10, 10, 10],
            'target_end': [990, 1990, 2990, 3990, 4990, 5990],
        }
        df = pd.DataFrame(data)

        # Call the function to find chimeric reads
        chimeric_reads = find_chimeras(df)

        # Check the expected output
        expected_chimeras = ['read=1', 'read=2', 'read=3', 'read=4', 'read=5', 'read=6']
        self.assertEqual(expected_chimeras, chimeric_reads)


class TestRemoveFromFile(unittest.TestCase):
    def setUp(self):
        self.test_file = "test_file.txt"
        with open(self.test_file, "w") as f:
            f.write("Line 1\n")
            f.write("Line 2\n")
            f.write("Chimera Line 1\n")
            f.write("Chimera Line 2\n")
            f.write("Line 5\n")

    def tearDown(self):
        if os.path.exists(self.test_file):
            os.remove(self.test_file)

    def test_remove_chimeras(self):
        chimeras = ["Chimera"]
        remove_from_file(chimeras, self.test_file)

        with open(self.test_file, "r") as f:
            lines = f.readlines()

        expected_lines = ["Line 1\n", "Line 2\n"]
        self.assertEqual(expected_lines, lines)

    def test_file_not_found(self):
        chimeras = ["Chimera"]
        nonexistent_file = "nonexistent_file.txt"
        with self.assertRaises(FileNotFoundError):
            remove_from_file(chimeras, nonexistent_file)


class MakeChimeraTestCase(unittest.TestCase):
    def test_make_chimera(self):
        flength = 1000
        mlength = 150
        rlength = 1200

        result = make_chimera(flength, mlength, rlength)

        # Check the length of the resulting sequence
        expected_length = flength + mlength + rlength
        self.assertEqual(len(result), expected_length)

        # Check if the middle sequence is correctly positioned
        expected_middle = result[mlength]
        self.assertTrue(expected_middle in result)

        # Check if the reverse complement sequence is correctly positioned
        expected_reverse = Seq(result[rlength]).reverse_complement()
        self.assertTrue(expected_reverse in result)


class AddErrorsTestCase(unittest.TestCase):

    def test_add_errors(self):
        sequence = Seq("ATCGATCGATCG")
        result = add_errors(sequence)

        # Check if the length of the resulting sequence is the same as the original
        self.assertEqual(len(result), len(sequence))

        # Check if the resulting sequence contains at least one substituted base
        self.assertTrue(result != sequence)


class MakeNormalTestCase(unittest.TestCase):

    def test_make_normal(self):
        length = 10
        result = make_normal(length)

        # Check if the generated sequence has the correct length
        self.assertEqual(len(result), length)

        # Check if the generated sequence contains only valid bases
        bases = set(result)
        expected_bases = {"A", "C", "T", "G"}
        self.assertTrue(bases.issubset(expected_bases))


class CheckChimeraFlatTestCase(unittest.TestCase):
    def test_check_chimera_flat(self):
        # Test case 1
        read = "ATCGTACGCTAGCTAGCTAGCTAGCTACGATCGATCGATCGATCGATCGATCGATCGATCG"
        interval = 10

        result = check_chimera_flat(read, interval)

        # Assert that the result is within the valid range of alignment scores
        self.assertGreaterEqual(result, 0)
        self.assertLessEqual(result, 1.25)

        # Test case 2
        read = "CTCGCTCGCTCGCTCGCTCGCTCGCTCGCTCGCTCGCTCGCTCGCTCGCTCGCTCGCTCGCT"
        interval = 5

        result = check_chimera_flat(read, interval)

        # Assert that the result is within the valid range of alignment scores
        self.assertGreaterEqual(result, 0)
        self.assertLessEqual(result, 1.25)

class FileReaderTestCase(unittest.TestCase):
    def test_file_reader(self):
        # Define a sample input file content
        input_file_content = """
        read=1 AAAATCGTAGCTGCTAGC
        read=2 CGTACGATCGTACGAT
        read=3 
        read=4 ATCGATCGATCGATCGA
        """

        # Create a temporary input file
        with open("temp_input_file.txt", "w") as file:
            file.write(input_file_content)

        # Call the file_reader function to read the temporary input file
        reads, sequences = file_reader("temp_input_file.txt")

        # Define the expected output
        expected_reads = ["read=1", "read=2", "read=4"]
        expected_sequences = ["AAAATCGTAGCTGCTAGC", "CGTACGATCGTACGAT", "ATCGATCGATCGATCGA"]

        # Assert the actual output matches the expected output
        self.assertEqual(expected_reads, reads)
        self.assertEqual(expected_sequences, sequences)

        # Remove the temporary input file
        os.remove("temp_input_file.txt")


class SelfAlignedTestCase(unittest.TestCase):
    @patch('Self_Aligned.file_reader')
    @patch('Self_Aligned.check_chimera_flat')
    def test_self_aligned(self, mock_check_chimera_flat, mock_file_reader):
        # Mock the return values for file_reader
        mock_file_reader.return_value = (['read=1', 'read=2', 'read=3'],
                                         ['AGCTATCGATCGATCG', 'ACGTACGTATCGATCGATCG', 'ATCGATCGATCG'])

        # Mock the return values for check_chimera_flat
        mock_check_chimera_flat.side_effect = [0.5, 0.7, 0.4]

        # Expected chimeras based on the mock values
        expected_chimeras = ['read=2']

        # Call the function with test input
        chimeras = self_aligned('test_file.txt', interval=75, minlength=4)

        # Assert the expected chimeras
        self.assertEqual(expected_chimeras, chimeras)

        # Assert that file_reader was called with the correct arguments
        mock_file_reader.assert_called_once_with('test_file.txt')

        # Assert that check_chimera_flat was called with the correct arguments for each sequence
        expected_calls = [
            call('AGCTATCGATCGATCG', 75),
            call('ACGTACGTATCGATCGATCG', 75),
            call('ATCGATCGATCG', 75)
        ]
        mock_check_chimera_flat.assert_has_calls(expected_calls)

if __name__ == "__main__":
    unittest.main()
