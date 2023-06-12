#!/usr/bin/env python3

"""
Mapping_Method.py script to find chimeras with help from a PAF file.
Required: Minimap PAF file.
Author: IJsbrand Pool
Version: 1.2.1
Date: 09-06-2023
"""

# Imports
import pandas as pd
import subprocess
import os


# Call minimap
# Only when used seperately
def create_df(fastq):
    """
    Create a DataFrame from a FASTQ file using minimap2 overlaps.

    Args:
        fastq (str): The path to the FASTQ file.

    Returns:
        pandas.DataFrame: The DataFrame containing the important columns.

    Raises:
        subprocess.CalledProcessError: If an error occurs while executing the 'minimap2' command.
        FileNotFoundError: If the file path given does not exist.
    """
    if os.path.exists(fastq):
        try:
            # Find the minimap2 executable
            minimap_executable = subprocess.run('find / -type f -name "minimap2" -executable -print -quit 2>/dev/null',
                                                shell=True, capture_output=True, text=True)
            minimap_location = minimap_executable.stdout.strip()

            # Execute minimap2 command to generate overlaps.paf file
            minimap_command = f"{minimap_location} -x ava-ont {fastq} {fastq} > overlaps.paf"
            subprocess.run(minimap_command, shell=True, check=True)

            # Define column names for the DataFrame
            columns = ['query_id', 'query_length', 'query_start', 'query_end', 'strand', 'target_id',
                       'target_length', 'target_start', 'target_end']

            # Read the CSV file and create the DataFrame
            df = pd.read_csv("overlaps.paf", sep='\t', usecols=range(9), header=None, names=columns)
            return df
        except subprocess.CalledProcessError:
            raise subprocess.CalledProcessError("An error occurred while executing the 'minimap2' command.")
    else:
        print("The file path given does not exist. Please check if you entered the correct path.")


# Read paf file
def find_chimeras(df):
    """
    Find chimeric reads in the given DataFrame.

    Args:
        df (pd.DataFrame): DataFrame containing read information.

    Returns:
        list: List of chimeric read names.

    Raises:
        KeyError: If the required columns are missing in the DataFrame.

    """
    try:
        # Get unique read names
        unique_reads = df["query_id"].unique()

        # Initialize variables
        chimeras_count = 0
        chimeras_list = []

        for read in unique_reads:
            # Filter matching reads with the same query and target name
            matching_reads = df[(df['query_id'] == read) & (df['target_id'] == read)]

            if len(matching_reads["query_length"].unique()) >= 2:
                # Skip and remove the read if it has multiple sequence lengths, probably a minimap error
                df = df[df.query_name != read]
                df = df[df.target_name != read]
                unique_reads = unique_reads[unique_reads != read]
                continue

            if len(matching_reads) >= 1:
                # Ignore reads shorter than the threshold
                if matching_reads['query_length'].unique() >= 2000:
                    strand_values = matching_reads['strand'].unique()

                    if strand_values == "-":
                        # Check coverage for forward strand
                        for index, match_row in matching_reads.iterrows():
                            coverage_ratio = (match_row[3] - match_row[2]) / match_row[1]
                            if coverage_ratio >= 0.5:
                                chimeras_count += 1
                                chimeras_list.append(read)
                    elif len(strand_values) == 2:
                        # Check if the reverse complement match has full coverage
                        for index, match_row in matching_reads.iterrows():
                            if (match_row[3] - match_row[2]) / match_row[1] >= 0.5:
                                chimeras_count += 1
                                chimeras_list.append(read)

        print("Total chimeras found:", chimeras_count)
        print("Chimeric read list:", chimeras_list)
        return chimeras_list

    except KeyError as e:
        print("Error: The provided DataFrame is missing the required columns:", e)
        return 0


def mapping_method(fastq):
    """ Creates the DataFrame and finds the chimeras """
    df = create_df(fastq)
    chimeras = find_chimeras(df)
    return chimeras
