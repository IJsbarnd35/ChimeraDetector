import random
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from Bio.Align.AlignInfo import SummaryInfo
from docopt import docopt
import re
import time


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
    #amount = int(len(read) * (random.randint(1, 10) / 100))
    amount = int(len(read) * 0.1)
    read = list(read)
    positions = random.sample(range(0, len(read)), amount)
    for pos in positions:
        # Kan dezelfde zijn, volgens mij niet erg, voegt meer random toe??
        error = random.choice(["A", "C", "T", "G", "", "insertion"])
        if error == "insertion":
            read[pos] = random.choice(["A", "C", "T", "G"]) + random.choice(["A", "C", "T", "G"])
        else:
            read[pos] = error
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


def check_chimera_old2(read):
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
        print(alignment)
        #print(identity)

        alignment = next(aligner.align(read, reverse_slice))
        identity = alignment.score / len(reverse_slice)
        scores_reversed.append(identity)
        print(alignment)
        bpos = alignment.aligned[0][0][0]
        epos = alignment.aligned[0][0][-1]
        print(bpos, epos)
        #print(identity)

        old_cut = new_cut
        new_cut += 0.05

    print(scores_forward)
    print(scores_reversed)

    if all(x > 0.75 for x in scores_reversed) or all(x > 0.75 for x in scores_forward):
        return read

def check_chimera_5p(read):
    # Figure out how many and how long the parts should be
    old_cut = 0.00
    new_cut = 0.05
    mode = "forward"

    matches_forward = []
    matches_reversed = []

    # Make the aligner
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    # Assign costs
    aligner.mismatch_score = -0.5
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5

    read = str(read)

    #for part in range(parts):
    while new_cut <= 0.20:
        if mode == "forward":
            slice = str(Seq(read[int(len(read) * old_cut):int(len(read) * new_cut)]).reverse_complement())
        elif mode == "reverse":
            if old_cut == 0.00:
                slice = str(Seq(read[-int(len(read) * new_cut):]).reverse_complement())
            else:
                slice = str(Seq(read[-int(len(read) * new_cut):-int(len(read) * old_cut)]).reverse_complement())

        # Align
        alignment = next(aligner.align(read, slice))
        identity = alignment.score / len(slice)
        #print(alignment)

        bpos = alignment.aligned[0][0][0]
        epos = alignment.aligned[0][-1][-1]
        #print(bpos, epos)
        # print(identity)

        if mode == "forward":
            matches_forward.append(identity)
            #matches_forward.append([identity, bpos, epos])
            mode = "reverse"
        elif mode == "reverse":
            matches_reversed.append(identity)
            #matches_reversed.append([identity, bpos, epos])
            old_cut = new_cut
            new_cut += 0.05
            mode = "forward"

    #print(matches_forward)
    #print(matches_reversed)
    print(sum(matches_forward)/len(matches_forward))
    print(sum(matches_reversed)/len(matches_reversed))


def check_chimera_1p(read):
    # Figure out how many and how long the parts should be
    old_cut = 0.00
    new_cut = 0.01
    mode = "forward"

    matches_forward = []
    matches_reversed = []

    # Make the aligner
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    # Assign costs
    aligner.mismatch_score = -0.5
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5

    read = str(read)

    #for part in range(parts):
    while new_cut <= 0.20:
        if mode == "forward":
            slice = str(Seq(read[int(len(read) * old_cut):int(len(read) * new_cut)]).reverse_complement())
        elif mode == "reverse":
            if old_cut == 0.00:
                slice = str(Seq(read[-int(len(read) * new_cut):]).reverse_complement())
            else:
                slice = str(Seq(read[-int(len(read) * new_cut):-int(len(read) * old_cut)]).reverse_complement())

        # Align
        alignment = next(aligner.align(read, slice))
        identity = alignment.score / len(slice)
        #print(alignment)

        bpos = alignment.aligned[0][0][0]
        epos = alignment.aligned[0][-1][-1]
        #print(bpos, epos)
        # print(identity)

        if mode == "forward":
            matches_forward.append(identity)
            #matches_forward.append([identity, bpos, epos])
            mode = "reverse"
        elif mode == "reverse":
            matches_reversed.append(identity)
            #matches_reversed.append([identity, bpos, epos])
            old_cut = new_cut
            new_cut += 0.01
            mode = "forward"

    #print(matches_forward)
    #print(matches_reversed)
    print(sum(matches_forward)/len(matches_forward))
    print(sum(matches_reversed)/len(matches_reversed))


def check_chimera_flat(read):
    # Figure out how many and how long the parts should be
    if len(read) < 500:
        return 0

    if len(read) < 2000:
        interval = 20
    elif len(read) < 5000:
        interval = 100
    else:
        interval = 200
    interval = 75

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

    #for part in range(parts):
    if len(read) < 500:
        x = int(len(read) * 0.5)
    else:
        x = int(len(read) * 0.25)
    bp = len(read)-x
    fnscore = 0
    rnscore = 0
    count = 0
    while new_cut <= x:
        if mode == "forward":
            if fnscore > 0.3 and count >= 5:
                mode = "reverse"
                continue
            slice = str(Seq(read[old_cut:new_cut]).reverse_complement())
            temp_read = read[new_cut:]
        elif mode == "reverse":
            if rnscore > 0.3 and count >= 5:
                old_cut = new_cut
                new_cut += interval
                mode = "forward"
                continue
            if old_cut == 0:
                slice = str(Seq(read[-new_cut:]).reverse_complement())
                temp_read = read[:-new_cut]
                #slice = str(Seq(read[bp + old_cut:bp + new_cut]).reverse_complement())
            else:
                #slice = str(Seq(read[bp + old_cut:bp + new_cut]).reverse_complement())
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
            if identity < 0.75: # and count > 2:
                fnscore = fnscore + (0.75 - identity)
            elif fnscore > 0 and identity > 0.8:
                fnscore = fnscore - (identity - 0.8)
            matches_forward.append(identity)
            mode = "reverse"
        elif mode == "reverse":
            positions_reversed.append([bpos, epos])
            if identity < 0.75: # and count > 2:
                rnscore = rnscore + (0.75 - identity)
            elif rnscore > 0 and identity > 0.8:
                rnscore = rnscore - (identity - 0.8)
            matches_reversed.append(identity)
            old_cut = new_cut
            new_cut += interval
            mode = "forward"
            count += 1

    #print(matches_forward)
    #print(matches_reversed)
    forward_score = sum(matches_forward)/len(matches_forward)
    reversed_score = sum(matches_reversed)/len(matches_reversed)
    #print("forward :", forward_score)
    #print("reversed :", reversed_score)
#
    #print(positions_forward)
    #print(positions_reversed)

    total_val = 0
    print(forward_score, reversed_score)
    if forward_score >= 0.6:
        #validate
        valscore = 0
        for match in range(len(matches_forward)-1):
            if abs(positions_forward[match][1] - positions_forward[match+1][0]) < 5:
                valscore += 1
        valscore = valscore/(len(matches_forward)-1)
        total_val += valscore
        forward_score = forward_score + valscore/4
    if reversed_score >= 0.6:
        # validate
        valscore = 0
        for match in range(len(matches_reversed) - 1):
            if abs(positions_reversed[match][1] - positions_reversed[match + 1][0]) < 5:
                valscore += 1
        valscore = valscore/(len(matches_reversed)-1)
        total_val += valscore
        reversed_score = reversed_score + valscore/4

    print("valscore: ", total_val/2)
    chimeric_score = (forward_score + reversed_score) / 2
    print(chimeric_score)
    print("-------------------")
    return chimeric_score



def testbig():
    start = time.time()
    for x in range(10):
        read = make_chimera(3000, 2000, 1000)
        #read = make_chimera(500, 300, 200)
        #read = make_chimera(100, 75, 25)
        #read = make_normal(200)
        check_chimera_flat(read)
        print("---------------------------")
    end = time.time()

    print(end - start)
    print("==================================================")
    start = time.time()
    for x in range(10):
        read = make_chimera(3000, 2000, 1000)
        #read = make_chimera(500, 300, 200)
        #read = make_chimera(100, 75, 25)
        #read = make_normal(200)
        check_chimera_5p(read)
        print("---------------------------")
    end = time.time()

    print(end - start)
    print("==================================================")
    start = time.time()
    for x in range(10):
        read = make_chimera(3000, 2000, 1000)
        # read = make_chimera(500, 300, 200)
        #read = make_chimera(100, 75, 25)
        # read = make_normal(200)
        check_chimera_1p(read)
        print("---------------------------")
    end = time.time()

    print(end - start)


#read = make_chimera(4500, 1000, 4500)
#read = make_chimera(40, 0, 60)
#read = make_chimera(100, 150, 250)
#read = make_chimera(50, 0, 50)
#read = make_normal(6000)
#print(read)
#start = time.time()
##testbig()
#check_chimera_flat(read)
#end = time.time()
#
#print(end - start)


def make_test_file():
    file = open("D:/School/stage2/data/testdata.fq", "w")
    for x in range(20):
        file.write(str(make_chimera(100, 0, 20)) + "\n")
    for x in range(10):
        file.write(str(make_chimera(60, 0, 60)) + "\n")
    for x in range(20):
        file.write(str(make_chimera(20, 0, 100)) + "\n")
    for x in range(50):
        file.write(str(make_normal(120)) + "\n")
#make_test_file()

def test_test_file():
    start = time.time()
    chimcount = 0
    for line in open("D:/School/stage2/data/testdata.fq", "r"):
        chimcount += check_chimera_flat(line)
    print(chimcount)
    end = time.time()
    print(end - start)
#test_test_file()


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
    # Orders the matches with each other.
    for index in range(len(reads)):
        myread = reads[index]
        # Splits sequence in lines of 70 nucleotides.
        myseq = (re.sub("(.{70})", "\\1\n", sequences[index], 0, re.DOTALL))
    return reads, sequences
#file_reader("D:/School/stage2/data/b1_1.fq")


#chimcount = 0
#seqcount = 0
#chimlist = []
def test_file(file):
    #reads, sequences = file_reader("D:/School/stage2/data/testdingfile.fq")
    #reads, sequences = file_reader("D:/School/stage2/data/tests/test1.fq")
    reads, sequences = file_reader("D:/School/stage2/data/testdingfile.fq")
    start = time.time()

    chimeras = []
    count = 0
    chimscore = 0
    normscore = 0
    for seq in sequences:
        score = check_chimera_flat(seq)
        if count % 2 == 0:
            chimscore += score
        else:
            normscore += score
        if score > 0.6:
            chimeras.append(reads[count])
        count += 1

    print(chimeras)
    print("chimeras: ", chimscore)
    print("normals:  ", normscore)
    end = time.time()
    print(end - start)
test_file("iets")


# with open("D:/School/stage2/data/testdingfile.fq", "r+") as f:
#    d = f.readlines()
#    f.seek(0)
#    count = 0
#    while count < len(d):
#        if any(chimera in d[count] for chimera in chimeras):
#            count += 4
#            continue
#        else:
#            f.write(d[count])
#        count += 1
#    f.truncate()



#for seq in sequences[:100]:
#    thing = check_chimera_flat(seq)
#    if thing == 1:
#        chimcount += 1
#        chimlist.append(seqcount)
###    if len(seq) > 500:
###        print("seqn = " + str(count))
###        print(len(seq))
###        check_chimera(seq)
###        print("-----------------------")
#    seqcount += 1
#


def check_chimera_by_readlength(length):
    reads, sequences = file_reader("D:/School/stage2/data/b1_1.fq")
    sding = ""
    coun = 0
    for seq in sequences:
        if len(seq) == length:
            print(reads[coun])
            print(check_chimera_flat(seq))
        coun += 1
#check_chimera_by_readlength(6736)




seqn1 = "AGTATACTTCGTTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTGGTGCTGAAGAAAGTTATCGGTATTTCTTTGTGATCAGCCGTAG"
seqn3 = "AAGTGTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCCAGCACCTAAGGTTAACACAAAGACACAAACAACTTTCTTCAGCACCTGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGAACATAGCTGTATGGCG"
#check_chimera(seqn1)


if __name__ == '__main__':
    doc = """
    Usage:
      comtest.py <file>
      comtest.py (-h | --help)
    Options:
      -h, --help
    """

    #arguments = docopt(doc)
    #file_reader(arguments["file"])
    #print(arguments)

def get_seq_by_readname(read):
    seqdict = {}
    for x in range(len(sequences)):
        seqdict[reads[x]] = sequences[x]
    print(len(seqdict))
    print(seqdict[read])

#get_seq_by_readname("read=96519")


def get_seq_by_len(leng):
    for x in sequences:
        if len(x) == leng:
            print(x)

#get_seq_by_len(2683)

seq = "ATGTACTTCGTTCAAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTCCATCTCCTCGTGGTATATCCTGTCAAATATACCTCACCATCAACATCAGTAATACTAAATCCAGTAGATTTAACCGTTCCACCACCTTCTAAAATGTGAAATTGGTTACCAAAACATAATTCATACTGAGCTCAGACTGATTTAATAGTGCTTTTGGGTCTACGGATGGTAACTCGTGTAATATTGATAATCGCAGGATCTGTTTCATCAACAATTCGTTGACTCTCAGAATACTTAAATCTACCCCAAATGCATTCAATTTAGCAGATTTGCTGTATTCTGTAAGAGAATTAACTACTTCTGCCTTAAGTTGGTCTTTATCGTCAATAATATGACATTATAGTAGACTGAAATGTCTAATTCTACATACAAAATCTTCAAATCTTCAATTCTTTTGGTTAATACCAGCAATTGCATATTGTTTTAAGTCATTAAGGATGTTCTGCTTGGTAAAGTCAGATAGAAATGTACCATTTTTCGGTTTATACTCAAAACTACCGTTCCAAACTCAGGAGGATCAATTCTTCACCTCCAACTACGAAACGGACTCCGTATTTGGATAAACAGACTGAATTGTTGCCTCATAATCCTAGGTGTAACTGCTCTGTTCTGTGCAGAGTAGATTCTAGAGTGCAAAATATCAGATAGAGTCAATTTGCTCAATATCAGACCCATTTCCTTGGTCTGAGTAGCAGTTACTATTACTACATTAGTAGATGAAAGAATTTCCTGCATCATCAACTATATATTTCAGAAAAAGTGAAGAATTTTCCTTCATTCCCTCTTTTCCGTCTGTTACTATATGAATATATCAATTTCATCGCCAATTTCTAGTTTTTTGCGCAAAAAGACCATCTCCAAAGAGCAATTCATAGGTTTCATTCATCGCTTCTTGCACCAAATAGACATTTGATACGAAGTAACATCAATAATGTTGTCAATTTTTTGAAAATTTAAGTCAGCAGAGGATCCAGACTTCTTAACAGTGACCCTAAAAGCGTATCTAAATCGACAAATGAGTTCTCAAGAATGAATCTTTGGTCAGAACTACCATTTACAGTGAATTTTTAGTTAAAAATTGTACCTTGATGAACAGTTACACCTGAAAATGTTGCTGTTCTAGGTGAATTATTACCTTATACCTCCAGAATCTAAAGGACTGGTGACAGATACATCTTCTGGGACAGAAAACTACTAACCTTGTTTTTCACGGCACCAACTAACGCTAATCCTTGCTTTAAAGTGACAGTAATACTGCTGAGTGAATTTATGTTGAAAATCAACGACTGCTTGAGCACATTTCTTGATCTAGGTACATATCCTATGTTTCTTGCAAGAGAAACTACATTCTCTCGTAGTGTTGCAGAATCAAAAAGGACTCGTGCTTAATATATTGGTATTAAATGCTACGCAATATATACAATTATGCAAAATATGTTAAAATCGACATATTTGACCCTTCAAAGTCAAAATCAAATCGGAAGTTAGCTCCTCAAATTAATCTTCTTTAGATTCATAATTTCATTGAAATTTAAATTAGTAAACTTAGTAAATGGCATTCTTATTACCTAGAAGGTTTCAACATAAAGTCAAATTCCTGACTTGGGAAGTCTTGACCTATAATATCGTATGTAATGAAGACTTCAAATGTATTATCGTTGGGTTTAGGTTCAACTCTTACTGCAGTATTATCAATACGACCTCAAACCCTTTAAAACATCAAAAATTTGCTGAGTTGTAAGCACTTGCAGTACCCATAATCAACAAAATCGAATAAAGATCTAGTCACATCTGTACCAGATTAGACGAGAATGGTCTTTCTCCACCTATAGTCTGAATCAAATTCCTCACAGCTCTCTTAATCGAATCCTCATTCTGCTTAAGGACTCGGCCACATCTCCTGAAACTGGGTGGCATAAAAAGAAGATCTATGTCCTTAAATGCCTTGGAATTCGATTCAATTTGGCCATACAGTGAAACCGGTATCGAGATTATTTATACCTTATGCACGGGATCATCATTGTGAATTACCCATATCCACTTTCCTGCTATCGGATCGTTTTCTGAAGAACCTCCCTTCTGGCATATCGAGCACCTTTTGACACCTATTAGCAGCACCTTCGTGATATCTTCCTCTATCTCTAATACAATTCATTGCAGTATTGTATATGTCAAACTCTGATGCCTCTGATTAATGACATCAAAGATTGCATCGTCTAATGCAGATAGGGCATAGTCTTCCCATTTGGTCTTATCGTCTATCAGGGGTATTACTTTATTACTACTCATTTTCCCTGACCTCTATATGGTTTCCGAGCAGCGTTTCGGGCGGTGCTAGAGAAAATTACAGGTTCTTTCCTTGTCCTTGTCGGAGTCTTTTTTTGAGGTAGATTCAACCCAACTGCCTCCGAGCAACGATTTTTTAATCTTAGCCCATTAGTCTCCTTCATAGTGTATGATATTGTAGAAGGATCTGGGGAACCACACTCATGGAATTCCTAGAGCATAATCCCCATTAGATCCAGAAATACAGTTTCACTTATATCGGTGTGTTTCTCAACTCTTGATGGGTGTGTAAACTGTCTTTGACATTAGATAATTCTCATCTTCTCGTGACCTACTCTAACACGAGGATCGCACCATATCTCGAATCCTGCTTCTATAGCATCGAGACAGAAACTCACATCCTCTCCACACATATCTTGAACATCACCTGATTCAAAGACTTGCATCTTAGAACCAGACCAAGGATACTTCATATCTGCATGTTCAAATACACGCTTCTTGATCATACCATCCAAATTATAACGAATAAAGTAAGAAGTAGGCTCTTACGCTTGGTCATTGTCTCACAATTTCGCTGATTCATGACTCCACCATTATTACGGAAATCATACCTGTCTAACCAATGAGCAACAGAAGTAGTACGACCATCTTCAGTCATATACCAACCTGCAGCAATATCCTTGTCATAAGGACAGTTGTAGAGACTTATCCAGTACCGTAACATATATTAGTTGATAGTCATAGTTTTAATTTACCATCCCAGAATTGCTTGATCACAGACCCCTTAATACATTTGCTCAAGACACCCACATCTAGCAAAGTTTTACCATACTAGGTAGTCTTGCGATATCTGAATACTTGCACCATGTTGTACCAAATCAAAACAGGTTGGCGAAACTCTTTAAGAACTGAAAGACATCCACATTAGGCATACAGAATACAATTGCCTTACCTTTTAACATCTCCCATGCTGCGTCATGGTCCCATTCAGGTTCTTGAACTTTTTCGGCGTTTTCGCCTTTACAGTAAATCAACTATGGCCATAATGATTAATGTACTTCCAGTTATTATAACAGTATTATGCGCTGAGTCAATCTTCTTCTTTTTCTTCGAGGAATACACCATCACTTTCAAGTGTCATAGTTATCTCACTACCCCATACACCAATCGAAGTCGTTATATACTCTCGGGCAGATCTAATACAAGTTGATCCTCTACAGGGTCGAACCTTGTAGTAACTTTTATATTATGGAAATTTTTTCACTCAAACGAAACCTGTGTGTCGTTTTTATATATTGGATATTTTTTCTAGAGATATAGCAAGGTCGATCTGGGTCGTTTATAGCTTACAAAGGTTCCTTCGATTTAAACCGCATCCCCAACACACCACGATATAACATAAGGGCGAACATACTGCCAATTAGTGTTAATTAGTGAGTGTAATATGCATGATGCTGTTACATAGTGCCTCTAATCTGTCTTGGCGTGTAGGTAACAATCAAGGTCATCATATTCATCATAATTGTTATTAGTAACTGTGTGCTCGCTATATCATCGAGCAGTGTAATCTTCCATACCATTATTGATAAGAATCTGACATACGATTGCAGTTCCTCGGTGTTATTAACCGTATTATAGCATGATGTTAGTCAAGGTGTCAATAATTGCCGCTTTGAGTTGTTTATAAGAACTGTGAGCATATTGTGACCTGCGTATTTGACATTTTGTTCACTCCATTATTGCGCTCGCTAAGCGAGTATACGGCATGGAGTGAACAAAATGTCAACATGCAGATCTAATATACTCACGATTCTTATAGAGTAACTCAAAAACAGTATTGACACCTGACTAACATCATGCTATAATACAGTTAATAACACCGAGGGGAATACAATCATATATCAGGTATCATCATGGTATGGAAGATTGGCACTGCTCGATGATACATAGTAACACAGTGCAATATGATATGATGACCTTGATTATTTTTATGATTAGAGGCACTATGCTAACAAGCATCATATATTACATGTTAACGCTGAGTGGCGTATGTCGCGCTATGTGAGTGATGGTTTAAATCGAAGGAATTCTTGTAAATATAGGCGACCCGGATCGACTGCTGCTATATCTCTGAAATAAAGGCGACACGCAGGCGTTTAGTGAAAATTTCCGTAATATAAGAAGTTACTATAAGGGTGACCTATGAGATCAACTGTATTAAGATTTGCGAAAATATGCAACGACCTTGATTGAGTAATTGGACCTTGAAAGTGATAGTTATTCCTCGAAGAAAAGAGAGATTGACTCCACATACATATAATACTATTATAATAACTGAAATACATTAATCATTATGGCACAAAATTTAAACGAGAGCATAGCGCTGAACCCAGATGGGGACTATGACACGAATATAGGAAGATGGGCAATTGCATTCTGCATACCTGTCATGGATGTTTTTCAGTTCTTAAGGATTTCATTACAATATGTTTGATTTGGTACAACATGGTGCGAAGTATTCGGATATCACCAAGACTCTCTAACATGTAAACTGCTAGATGTAAGTGTTTCTTTAGGTACAAATGCTAAGGAATCCCGAATCAAGTACCTGGGATAATGGGAGTGAAACTATGACTATCGGCTATGGATTGATAGTGATATAGTATACGAGTACTGATAATTCTTACAACTTGTTCTTATGGACAAGGATAATTGCTGCGGAACAGTATATGACGAGATGAAGTCGTACTACTCTGTTGCTCATTGGTAATTTCCATGACAACTGAGTAGTCATGAATCACGAAGACTGGTTTTAATGATGGCAGTAAAGCCATTTACTGTGGGCTGCAGGTTTCGGATGGTGCTTGTTGGAAGGTGTGAACATGCGGATATGGTAACAGTTCACTCTAAAGATGCAAGTCTTTAGTCAGGCGATATTCAAAGATATGTGTGGGGAGAGTTTGCTGTCTCGATGCTGTAAATTTCGAGATATGTGGGAATCCTAAATGTTAAATATCACGAGATGGCGATATCAATGTCAAGACAGTTTACACAGTCTATCAGGAAGCGTTGAAAACTTGCACCGATATAAGTATGAGCTATGATTTCTGGATCTGAGAGGATTATAACTCAGTTTATGAGTGTGAAGTTCCCCAGATCCTTCTGTAAATATCATGTGCTATGAAGGACTAATGAAATAAAAGTCGTTGCAGCCGAAAGGAGGCAGGTGAATCTACTTTTAAGACTCGTGCAGGACGGAAGAACTCAAGACTCTAGCTTTGCCGCCCGAACTTGTCCTCATAACCATATACAGAGGTCGGGAAAGGTAATTGAGTAATGAAATGTGCCCACGATGATGGACCAGATGGGAGACTATGCCTATCTAAATTGGTGATGCGATCTTTGATGTCATTAATTCAGGGCAATCAGTACAATACTATAATGAATTGTATTGAGACAGAGGAAGATATCACGAAGGCTGCTGCTGTAATGAGGTGGCTCCGTGCTGGAGCTTTTTAGGGGGGCAGTCCGTGCAGAGAGGAAATTGAACATGAGTGGGTATATAATGTGTCCTGCATAGAAGAAATTATCGTATCTCATTACCTAAGAGCCCCTTTCTAATAGGCAGATCGAATTCCAAGGCGACTAGCATAGATCTTTTTCTTTACGCGAATTTCAGAGATATAGGATCACTTAAGAATGAGTTCGGTAAAGCGGCAATTTGATTCGAACATAGGTGGAAAACCATTCTTTAATTGCAGATTTAGATCTTTTATTCATTTTATTGATATGTACTGCGAATGTCATAACTCGGCTGGAATTTTGATGTTTTAAAAGGGTTTGAAAATTAAGTATTAATACAGTAAAGTTGAACTACAAATGTGACTGTATTTGAGTCTTCATTACACATGCGGTAATATGGTAGAACTTCAAGTCGGAATTTGACTTTATATTTGAGTCTTCTAAACAAATGCCATTTACTAAGTTTACTAATTTAAAATTTCAATGAAATCTATTAAAAGATACCTGAAGCTGCTCTGATTTTCAAGATTTTGGCTTTTTGAAGGGTCAAGATATGTCGATTTTTAGTCGATATTTTGGCATATAGTTCGTATATTACAGCATTTAATACCAATATGAATTGCTAATGCTTCTCTCTGAAAGAATGTAGTTTCTCTTGCAGAAACATAGGATATGTGCCTAGATCAGAAAGTGTAACTGAAACTGTCATTGATTTTTCGAAAAATTCACTGGAAATAATAATACTGTCACTTTTAAAGCAGGTGTAATTGAATAACTATAAAAAATACAGGATTATATGTTTCTGTCCCGAAGATATATCTGTCATGATCCTTTAGATTCTGGAGACAAAATGATAATATTCGCACAGAACAACGACATTTCAGGTATGGCTATCTATCAAGATACACTTTGAACTAAAATTCACTGTAAATAGTTCTATTGAATTCGTCTTAATTTATCATTTAGATACAACCCCAAATTAAGTCTGGATCCTATAACAGACTAAATTTTCAAAATTGACAACATTATTGATATTGCCCCTATTGTATCAAATGTCTATTTGATTACGAAGGTAAGAATGAGAACCTATGAATTGCTCTTTTGGGAATGGTCTTTGGTAAGGCGATGAGTATTATATCTTACATAATGGCTGGACGGAAAGAGAAAGGAAGTAATTCACTCTAGATATAGTTGATGATGCAGGAGAACTATAACAAATCTACTAATGTAACAACGGCAACTGCTTCTCGAACTTATGGGAATAGTGCGGTGGTGAACGGAGATTGACTCTATTCGATATTTTGCCTAGAATGTACTGCGCAGAACAGAGCAGTACACAGGGTTGTGGGCAATGAATTCGGTCTGTTTTATCAAATGCGGAGTCCCGTCTTTCTGTAGTTGAAAATTGAAGAATTGGATCCTCTGAACTGGAACGACAGTTTGGTATAAAACCGAGAGAATGGTACACTGTGACTTTTACCAAGCGGGGTATCTTAAATGACTTAAGACAATAATGCAATTGCTATTAACAAAGAGTGAGATTTGAACAAGGAGTAGACACTTCAGTCTACTATAACGGCGAGTATTGACGACGAACTGGCTTAGACGAGTAGTTAATTCTCTTGCAAGAATACGGCAAGATCTAGAGTGAATGCATTTGAAGTAGACTAAGTATTCTGAGTCAACAGATATTGATGAAACAGATCCTGCGATTATAATGTACACGAGTTACAAATCCGTAGAACTAAAGCACTATTAGATGTTCTCTTCGGTATGAATTATGTTTTGGTAACCAATTTCACATTTTAGAAGATGAAGTGGAGATTTGGTGGAAATCTACTGGACATTACTGATATTTGATAATTGTGACTACGACAGATATACGAGAACTGGAGGTGCTGAAGAAGTTGTCGGTAACAACTGTGTTAACCTTAGCAATACGCTGAAAAAA"
#heck_chimera_flat(seq)