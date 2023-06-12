#!/usr/bin/env python3
"""
Create_Artificial_Reads.py: Created for artificial read creation.
Author: IJsbrand Pool
Version: 1.3.1
Date: 09-06-2023
"""

import random
from Bio.Seq import Seq


def make_chimera(flength, mlength, rlength):
    """
    Generate a chimera DNA sequence.

    Args:
        flength (int): Length of the forward sequence.
        mlength (int): Length of the middle sequence.
        rlength (int): Length of the reverse sequence.

    Returns:
        str: The generated chimera DNA sequence.

    Raises:
        None
    """

    # Available DNA bases
    bases = ["A", "C", "T", "G"]

    # Generate the middle sequence
    middle = ''.join([random.choice(bases) for _ in range(mlength)])

    if rlength > flength:
        # Generate the full forward sequence and extract the required length
        fullforward = ''.join([random.choice(bases) for _ in range(rlength)])
        forward = fullforward[-flength:]
        reverse_complement = Seq(fullforward).reverse_complement()
    else:
        # Generate the forward sequence and extract the required length
        forward = ''.join([random.choice(bases) for _ in range(flength)])
        reverse_complement = Seq(forward[-rlength:]).reverse_complement()

    # Concatenate the sequences to create the chimera DNA sequence
    read = Seq(forward + middle + reverse_complement)

    # Return the chimera DNA sequence with errors added
    return add_errors(read)


def add_errors(read):
    """
    Add errors to a DNA sequence by randomly substituting bases.

    This function takes a DNA sequence and introduces errors by randomly substituting
    bases at a rate of 10% of the sequence length. The substituted bases are chosen randomly
    from the set of DNA bases: A, C, T, G.

    Args:
        read (Bio.Seq.Seq): The DNA sequence to add errors to.

    Returns:
        Bio.Seq.Seq: The DNA sequence with errors.

    """
    # Calculate the number of bases to substitute (10% of sequence length)
    amount = int(len(read) * 0.1)
    read = list(read)
    positions = random.sample(range(0, len(read)), amount)

    for pos in positions:
        read[pos] = random.choice(["A", "C", "T", "G"])

    return Seq(''.join(read))


def make_normal(length):
    """
    Generate a random DNA sequence of the specified length.

    Args:
        length (int): The length of the DNA sequence to generate.

    Returns:
        str: The randomly generated DNA sequence.

    """
    bases = ["A", "C", "T", "G"]
    return ''.join(([random.choice(bases) for _ in range(length)]))
