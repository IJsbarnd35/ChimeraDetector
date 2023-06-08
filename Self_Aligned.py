import random
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import re
import time


def make_chimera(flength, mlength, rlength):
    bases = ["A", "C", "T", "G"]
    middle = ''.join([random.choice(bases) for x in range(mlength)])

    if rlength > flength:
        fullforward = ''.join([random.choice(bases) for x in range(rlength)])
        forward = fullforward[-flength:]
        reverse_complement = Seq(fullforward).reverse_complement()
    else:
        forward = ''.join([random.choice(bases) for x in range(flength)])
        reverse_complement = Seq(forward[-rlength:]).reverse_complement()

    read = Seq(forward + middle + reverse_complement)
    return add_errors(read)


def add_errors(read):
    amount = int(len(read) * 0.1)
    read = list(read)
    positions = random.sample(range(0, len(read)), amount)
    for pos in positions:
        read[pos] = random.choice(["A", "C", "T", "G"])
    return Seq(''.join(read))


def make_normal(length):
    bases = ["A", "C", "T", "G"]
    return ''.join(([random.choice(bases) for x in range(length)]))


def check_chimera_flat(read, interval):
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

    print(matches_forward)
    print(matches_reversed)

    forward_score = sum(matches_forward)/len(matches_forward)
    reversed_score = sum(matches_reversed)/len(matches_reversed)
    print("forward :", forward_score)
    print("reversed :", reversed_score)

    print(positions_forward)
    print(positions_reversed)

    if forward_score >= 0.75:
        #validate
        valscore = 0
        for match in range(len(matches_forward)-1):
            if abs(positions_forward[match][1] - positions_forward[match+1][0]) < 5:
                valscore += 1
        valscore = valscore/(len(matches_forward)-1)
    elif reversed_score >= 0.75:
        # validate
        valscore = 0
        for match in range(len(matches_reversed) - 1):
            if abs(positions_reversed[match][1] - positions_reversed[match + 1][0]) < 5:
                valscore += 1
        valscore = valscore/(len(matches_reversed)-1)

    return max(forward_score, reversed_score, 0)


def file_reader(input_file):
    reads = []
    sequences = []
    """ Reads through the Fastq file and extracts the sequences and read numbers. """
    for line in open(input_file):
        # Isolates the sequences.
        sequence = re.compile("^[A-Z]{5,}")
        # Isolates the read numbers.
        read = re.compile("read=\d*")
        # Finds all the matches.
        match = (read.findall(line.strip()))
        seqmatch = (sequence.findall(line.strip()))
        # Removes empty matches
        if match != [] or seqmatch != []:
            reads.append(''.join(match))
            sequences.append(''.join(seqmatch))
        while '' in reads and '' in sequences:
            reads.remove("")
            sequences.remove("")
    return reads, sequences


def self_aligned(file, interval=75, minlength=2000):
    reads, sequences = file_reader(file)

    chimeras = []
    count = 0
    for seq in sequences:
        if len(seq) > minlength:
            if check_chimera_flat(seq, interval) > 0.6:
                chimeras.append(reads[count])
        count += 1

    return chimeras


#if __name__ == '__main__':
#    doc = """
#    Usage:
#      Self_Aligned.py <file>
#      comtest.py (-h | --help)
#    """
#    arguments = docopt(__doc__)
#    check_file(arguments["<file>"])