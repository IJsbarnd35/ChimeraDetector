#!/usr/bin/env python3
"""
Self_Aligned.py: Created for chimera detection.
Author: IJsbrand Pool
Version: 1.3.3
Date: 09-06-2023
"""

from Bio.Align import PairwiseAligner
import re
from Bio.Seq import Seq


def check_chimera_flat(read, interval):
    """
    Check for chimeric regions in a DNA read using pairwise alignment.

    This function divides the DNA read into intervals and performs pairwise alignment between
    subsequences to detect potential chimeric regions. It calculates alignment scores, identifies
    positions, and calculates scores for forward and reversed alignments. It also performs validation
    by checking the proximity of aligned positions.

    Args:
        read (str): The DNA read sequence.
        interval (int): The interval length for dividing the read.

    Returns:
        float: The maximum alignment score for forward, reversed, or validated alignments.

    """
    interval = interval

    old_cut = 0
    new_cut = interval
    mode = "forward"

    matches_forward = []
    matches_reversed = []

    positions_forward = []
    positions_reversed = []

    # Make the aligner
    aligner = PairwiseAligner()
    aligner.mode = 'local'

    # Assign costs
    aligner.mismatch_score = -0.5
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.5

    read = str(read)

    if len(read) < 500:
        x = int(len(read) * 0.5)
    else:
        x = int(len(read) * 0.25)
    fnscore = 0
    rnscore = 0
    count = 0
    while new_cut <= x:
        if mode == "forward":
            if fnscore > 0.3 and count >= 10:
                mode = "reverse"
                continue
            slice = str(Seq(read[old_cut:new_cut]).reverse_complement())
            temp_read = read[new_cut:]
        elif mode == "reverse":
            if rnscore > 0.3 and count >= 10:
                old_cut = new_cut
                new_cut += interval
                mode = "forward"
                continue
            if old_cut == 0:
                slice = str(Seq(read[-new_cut:]).reverse_complement())
                temp_read = read[:-new_cut]
            else:
                slice = str(Seq(read[-new_cut:-old_cut]).reverse_complement())
                temp_read = read[:-new_cut]
        # Align
        alignment = next(aligner.align(temp_read, slice))
        identity = alignment.score / len(slice)

        bpos = alignment.aligned[0][0][0]
        epos = alignment.aligned[0][-1][-1]

        if mode == "forward":
            bpos += new_cut
            epos += new_cut

        if mode == "forward":
            positions_forward.append([epos, bpos])
            if identity < 0.75:
                fnscore = fnscore + (0.75 - identity)
            elif fnscore > 0 and identity > 0.8:
                fnscore = fnscore - (identity - 0.8)
            matches_forward.append(identity)
            mode = "reverse"
        elif mode == "reverse":
            positions_reversed.append([bpos, epos])
            if identity < 0.75:
                rnscore = rnscore + (0.75 - identity)
            elif rnscore > 0 and identity > 0.8:
                rnscore = rnscore - (identity - 0.8)
            matches_reversed.append(identity)
            old_cut = new_cut
            new_cut += interval
            mode = "forward"
            count += 1

    forward_score = sum(matches_forward)/len(matches_forward)
    reversed_score = sum(matches_reversed)/len(matches_reversed)

    forward_valscore = 0
    reversed_valscore = 0
    if forward_score >= 0.6:
        #validate
        for match in range(len(matches_forward)-1):
            if abs(positions_forward[match][1] - positions_forward[match+1][0]) < 5:
                forward_valscore += 1
        forward_valscore = forward_valscore/(len(matches_forward)-1)
        forward_score = forward_score + forward_valscore/4
    if reversed_score >= 0.6:
        # validate
        for match in range(len(matches_reversed) - 1):
            if abs(positions_reversed[match][1] - positions_reversed[match + 1][0]) < 5:
                reversed_valscore += 1
        reversed_valscore = reversed_valscore/(len(matches_reversed)-1)
        reversed_score = reversed_score + reversed_valscore/4
    return max(forward_score, reversed_score, 0)


def file_reader(input_file):
    """
    Read an input file and extract reads and sequences.

    This function reads the input file line by line and extracts the read numbers and sequences.
    It uses regular expressions to isolate the sequences and read numbers from each line.
    The extracted reads and sequences are returned as lists.

    Args:
        input_file (str): The path to the input file.

    Returns:
        tuple: A tuple containing two lists - reads and sequences.

    """
    reads = []
    sequences = []

    # Iterate through each line in the input file
    for line in open(input_file):

        # Isolates the sequences using regular expression pattern matching
        sequence = re.compile("^[A-Z]{5,}")

        # Isolates the read numbers using regular expression pattern matching
        read = re.compile("read=\d*")

        # Finds all the matches for reads and sequences in the line
        match = read.findall(line.strip())
        seqmatch = sequence.findall(line.strip())

        # Remove empty matches and append the non-empty reads and sequences
        if match:
            reads.append(''.join(match))
        elif seqmatch:
            sequences.append(''.join(seqmatch))


        # Remove empty elements from the lists if present
        while match and seqmatch:
            reads.append(''.join(match))
            sequences.append(''.join(seqmatch))

    return reads, sequences


def self_aligned(file, interval=75, minlength=2000):
    """
    Identify self-aligned chimeric reads from an input file.

    This function reads an input file and identifies self-aligned chimeric reads by applying the following steps:
    1. Reads and sequences are extracted using the 'file_reader' function.
    2. For each sequence, if the length is greater than the specified 'minlength':
        - The 'check_chimera_flat' function is called to check the chimera score based on the given 'interval'.
        - If the chimera score is greater than 0.6, the read ID is appended to the 'chimeras' list.
    3. The 'chimeras' list containing the read IDs of self-aligned chimeric reads is returned.

    Args:
        file (str): The path to the input file.
        interval (int, optional): The interval used for checking chimera score. Defaults to 75.
        minlength (int, optional): The minimum length required for a sequence to be considered. Defaults to 2000.

    Returns:
        list: A list containing the read IDs of self-aligned chimeric reads.

    """
    reads, sequences = file_reader(file)

    chimeras = []
    count = 0
    for seq in sequences:
        if len(seq) > minlength:
            if check_chimera_flat(seq, interval) > 0.75:
                chimeras.append(reads[count])
        count += 1
    print(chimeras)
    return chimeras
