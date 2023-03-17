import random
from biopython.Seq import Seq
from biopython.Align import PairwiseAligner
import re


bases = ["A", "C", "T", "G"]


#def translate(rc):
#    translation_dict = {"A":"T", "T":"A", "C":"G", "G":"C"}
#    for base in range(len(rc)):
#        rc[base] = translation_dict[rc[base]]
#    rc = ''.join(rc)
#    rc = rc[::-1]
#    return rc


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
    amount = int(len(read) * (random.randint(1, 10) / 100))
    read = list(read)
    positions = random.sample(range(0, len(read)), amount)
    for pos in positions:
        # Kan dezelfde zijn, volgens mij niet erg, voegt meer random toe??
        read[pos] = random.choice(["A", "C", "T", "G"])
    return Seq(''.join(read))


def make_normal(length):
    bases = ["A", "C", "T", "G"]
    return ''.join(([random.choice(bases) for x in range(length)]))


def check_chimera_old(read):
    # Take the last min% of bases
    minpercent = 0.90
    last = read[int(len(read)*minpercent):]
    last = last.reverse_complement()
    print(last)

    read = str(read)
    last = str(last)

    aligner = PairwiseAligner()
    aligner.mode = 'local'
    # Assign costs
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.5
    print(aligner.score(read, last)) # [:int(len(read)*0.9)]

    # Align
    alignment = next(aligner.align(read, last)) # [:int(len(read)*0.9)]
    print(alignment)
    identity = alignment.score / len(last)
    print(identity)

    # Check if they align enough
    while identity > 0.9:
        # Take the next x% or n amount of bases
        minpercent = minpercent - 0.01
        biggerlast = read[int(len(read)*minpercent):]
        biggerlast = str(Seq(biggerlast).reverse_complement())
        # Next?
        alignment = next(aligner.align(read, biggerlast))
        identity = alignment.score / len(biggerlast)
        print(alignment.score)
        print(alignment)
        print(identity)


def check_chimera(read):
    # Figure out how many and how long the parts should be
    parts = 4
    old_cut = 0.00
    new_cut = 0.05

    scores_forward = []
    scores_reversed = []

    for part in range(parts):
        forward_slice = str(Seq(read[int(len(read) * old_cut):int(len(read) * new_cut)]).reverse_complement())
        if old_cut == 0.00:
            reverse_slice = str(Seq(read[-int(len(read) * new_cut):]).reverse_complement())
        else:
            reverse_slice = str(Seq(read[-int(len(read) * new_cut):-int(len(read) * old_cut)]).reverse_complement())
        #print("forward: " + forward_slice)
        #print("reverse: " + reverse_slice)

        read = str(read)

        # Make the aligner
        aligner = PairwiseAligner()
        aligner.mode = 'local'
        # Assign costs
        aligner.mismatch_score = -0.5
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = -0.5
        #print(aligner.score(read, forward_slice))
        #print(aligner.score(read, reverse_slice))

        # Align
        alignment = next(aligner.align(read, forward_slice))
        identity = alignment.score / len(forward_slice)
        scores_forward.append(identity)
        #print(alignment)
        #print(identity)

        alignment = next(aligner.align(read, reverse_slice))
        identity = alignment.score / len(reverse_slice)
        scores_reversed.append(identity)
        #print(alignment)
        #print(identity)

        old_cut = new_cut
        new_cut += 0.05

    print(scores_forward)
    print(scores_reversed)


read = make_chimera(500, 300, 200)
#read = make_chimera(100, 150, 250)
#read = make_chimera(250, 0, 250)
#read = make_normal(500)
print(read)
#check_chimera(read)

reads = []
sequences = []

def file_reader(input_file):
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
    # Orders the matches with each other.
    for index in range(len(reads)):
        myread = reads[index]
        # Splits sequence in lines of 70 nucleotides.
        myseq = (re.sub("(.{70})", "\\1\n", sequences[index], 0, re.DOTALL))

        # Prints the sequence (REMOVE LATER)
        #print(f"{myread} \n {myseq}".replace(' ', ''))

file_reader("D:/School/stage2/data/b1_1.fq")

count = 0
for seq in sequences:
    if(len(seq)) < 500:
        count += 1
print(count)
#for seq in sequences:
#    print("Sequence number: " + str(count))
#    print(len(seq))
#    check_chimera(seq)
#    count += 1
#    print('--------------------------------')




