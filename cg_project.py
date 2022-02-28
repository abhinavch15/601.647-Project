# Libraries
import numpy as np
import itertools
import timeit
import sys
import sqlite3

# Connect to the data base
conn = sqlite3.connect(r'C:\Users\jai2m\Downloads\protein.db')

# Take file name arguments
input_file_name_dna = sys.argv[1]

# Load SQL data

start = timeit.default_timer()

cursor = conn.execute("SELECT ID, NAME, TYPE, SEQUENCE from ProteinDB")
proteins = {}
for row in cursor:
    proteins[row[1] + ': ' + row[0]] = row[3]

stop = timeit.default_timer()

print('Time taken to load SQL DB: ', stop - start) 

# Generates the 6 reading frames of a DNA sequence

def forward_seq(seq, init_pos=0):
    """Splitting the DNA sequence in groups of 3"""
    return [[seq[pos:pos+3]] for pos in range(init_pos,len(seq)-2,3)]

DNA_ReverseComplement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
def reverse_complement(seq):
    """
    Swapping adenine with thymine and guanine with cytosine.
    Reversing newly generated string
    """
    return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]


def gen_reading_frames(seq):
    """Generate the six reading frames of a DNA sequence"""
    """including reverse complement"""
    
    frames = []
    frames.append(forward_seq(seq, 0))
    frames.append(forward_seq(seq, 1))
    frames.append(forward_seq(seq, 2))
    frames.append(forward_seq(reverse_complement(seq), 0))
    frames.append(forward_seq(reverse_complement(seq), 1))
    frames.append(forward_seq(reverse_complement(seq), 2))
    return frames

# Translating DNA reading frames to proteins 

DNA_Nucleotides = ['A', 'C', 'G', 'T']
DNA_ReverseComplement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

def reverse_complement(seq):
    """
    Swapping adenine with thymine and guanine with cytosine.
    Reversing newly generated string
    """
    return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]

DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}

def translate_seq(seq, init_pos=0):
    """Translates a DNA sequence into an aminoacid sequence"""
    return [DNA_Codons[seq[pos:pos+3]] for pos in range(init_pos,len(seq)-2,3)]

def gen_translated_reading_frames(seq):
    """Generate the six reading frames of translated DNA sequence"""
    """including reverse complement"""
    
    frames = []
    frames.append(translate_seq(seq, 0))
    frames.append(translate_seq(seq, 1))
    frames.append(translate_seq(seq, 2))
    frames.append(translate_seq(reverse_complement(seq), 0))
    frames.append(translate_seq(reverse_complement(seq), 1))
    frames.append(translate_seq(reverse_complement(seq), 2))
    return frames

"""
Sample Test Case
"""
#dna = 'AGCCCTCCAGGACAGGCTGCATCAGAAGAGGCCATCAAGCAGATCACTGTCCTTCTGCCATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAGACGCAGCCCGCAGGCAGCCCCACACCCGCCGCCTCCTGCACCGAGAGAGATGGAATAAAGCCCTTGAACCAGC'
#protein = 'MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN'

with open(input_file_name_dna) as f:
    dna = f.read()

#with open(input_file_name_protein) as f:
#    protein = f.read()

translated = [''.join(translation) for translation in gen_translated_reading_frames(dna)]

# K-mer algorithm

def kmer_match(k, translated, protein):
    """
    k: Size of kmer
    translated: Translated DNA frames
    protein: Protein to match against
    """
    n = len(protein)
    output = []
    for j in range(len(translated)):
        s = translated[j]
        kmer = {}
        for i in range(len(s)-k):
            if s[i:i+k] not in kmer: kmer[s[i:i+k]] = []
            kmer[s[i:i+k]].append(i)
        count = 0
        match = []
        if protein[:k] not in kmer: continue
        for i in kmer[protein[:k]]:
            count += 1
            if i + n < len(s) and s[i:i+n] == protein: match.append(i)
        output.append([count, match, j+1])
    return output

# Smith Waterman algorithm

def matrix(a, b, match_score=3, gap_cost=2):
    H = np.zeros((len(a) + 1, len(b) + 1),)
    for i, j in itertools.product(range(1, H.shape[0]), range(1, H.shape[1])):
        match = H[i - 1, j - 1] + (match_score if a[i - 1] == b[j - 1] else - match_score)
        delete = H[i - 1, j] - gap_cost
        insert = H[i, j - 1] - gap_cost
        H[i, j] = max(match, delete, insert, 0)
    return H

def traceback(H, b, b_='', old_i=0):
    H_flip = np.flip(np.flip(H, 0), 1)
    i_, j_ = np.unravel_index(H_flip.argmax(), H_flip.shape)
    i, j = np.subtract(H.shape, (i_ + 1, j_ + 1))
    if H[i, j] == 0:
        return b_, j
    b_ = b[j - 1] + '-' + b_ if old_i - i > 1 else b[j - 1] + b_
    return traceback(H[0:i, 0:j], b, b_, i)

def smith_waterman(a, b, match_score=3, gap_cost=2):
    a, b = a.upper(), b.upper()
    H = matrix(a, b, match_score, gap_cost)
    b_, pos = traceback(H, b)
    return pos, pos + len(b_)

def alignment(a,b):
    start, end = smith_waterman(a, b)
    if b[start:end] == a:
        return start
    return None

# Boyer Moore

def boyerMoore(txt, pat):
    def badCharHeuristic(string, size):
        badChar = [-1]*256
        for i in range(size):
            badChar[ord(string[i])] = i;
        return badChar

    m = len(pat)
    n = len(txt)
    badChar = badCharHeuristic(pat, m)
 
    s = 0
    while(s <= n-m):
        j = m-1
        while j>=0 and pat[j] == txt[s+j]:
            j -= 1
        if j<0:
            return s
        else:
            s += max(1, j-badChar[ord(txt[s+j])])
          
    return None


# FM Index

class FM_index():
    def __init__(self, firstColumn, lastColumn, SA, RefTable):
        self.lastColumn = lastColumn  # store the last column of FM index
        self.firstColumn = firstColumn  # store the first column of FM index
        self.SA = SA
        # self.Tally = Tally
        self.RefTable = RefTable
        self.match = []

    def find_index(self, letter, indexNow):
        nearestCheckPoint = indexNow // 5
        nearestCheckPoint = nearestCheckPoint * 5
        location_in_list = self.RefTable[letter]
        checkpointValue = self.Tally[nearestCheckPoint][location_in_list]

        if (nearestCheckPoint > indexNow):
            distance = self.lastColumn[indexNow: nearestCheckPoint]
            occurance = distance.count(letter)
            return checkpointValue - occurance
        elif (nearestCheckPoint < indexNow):
            distance = self.lastColumn[nearestCheckPoint: indexNow]
            occurance = distance.count(letter)
            return checkpointValue + occurance
        else:
            return checkpointValue

    def recur_search(self, p, f, i, index):
        # print("P now is ",p)
        # print("f now is",f)
        # print("i now is",i)
        # print("Index is ",index)
        if len(p) == 0:
            if i in self.SA:
                self.match.append(self.SA[i])
            else:
                self.match.append(index)
            return
        L_iter = f[0]
        next_index = f[1]
        if L_iter == p[0]:
            start_point = self.firstColumn[L_iter][0]
            if i in self.SA:
                index = self.SA[i]
            elif index != None:
                index -= 1
            self.recur_search(p[1: len(p)], self.lastColumn[start_point + next_index], start_point + next_index, index)
        else:
            return

    def search(self, p):
        p = p.replace("_", "")
        p = p[::-1]
        P_length = len(p)
        firstLetter = p[0]
        p = p[1: len(p)]
        # print(p)
        if firstLetter in self.firstColumn:
            lowerbond = self.firstColumn[firstLetter][0]
            upperbond = self.firstColumn[firstLetter][1]
            for i in range(lowerbond, upperbond):
                # print("This is the i",i)
                self.recur_search(p, self.lastColumn[i], i, None)

def suffixArray(s):
    satups = sorted([(s[i:], i) for i in range(len(s))])
    # Extract and return just the offsets
    return map(lambda x: (x[0], x[1]), satups)

# Tally attemp failed. No longer under consideration.
def Tally_matrix(Tally, bw, refTable):
    checkPoint5 = {}
    arraysize = len(refTable)
    total_length = len(bw)
    for i in bw:
        if Tally[i][0] == False:
            Tally[i][0] = True;
            Tally[i].append(0)
        else:
            Tally[i].append(Tally[i][len(Tally[i]) - 1] + 1)
            for key in Tally:
                if key != i:
                    if Tally[key][0] == False:
                        Tally[key].append(0)
                    else:
                        Tally[key].append(Tally[key][len(Tally[key]) - 1])

        # print("A  is ",len(Tally['A']))
        iteration_length = total_length // 5
        for iter in range(0, iteration_length + 1):
            checkPoint5[iter * 5] = [0] * arraysize
            for key, value in Tally.items():
                if value[0] == True:
                    value.pop(0)
                else:
                    checkPoint5[iter * 5][refTable[key]] = value[iter * 5 - 1]

        return checkPoint5
            # Run all algorithms on given DNA and protein data



def buildFM(t):
    """ Given T, returns FM Index Object. """
    bw = []
    fr = {}
    refTable = {}
    SA = []
    new_SA = {}
    last_column = {}
    t = t.replace("_", "")
    t = t+'$'
    current = 0
    counter = 0
    for si in suffixArray(t):
        SA.append(si[1])
        if si[1] == 0:
            bw.append(("$", 0))

        else:
            l = t[si[1] - 1]
            if l in last_column:
                last_column[l] = last_column[l] + 1

            else:
                last_column[l] = 0
            # print("L is",l)
            bw.append((l, last_column[l]))
        letter = si[0][0]
        # print("letter is ", letter)
        # print(current)
        if letter in fr:
            fr[letter] = (fr[letter][0], fr[letter][1] + 1)
        else:
            if current == 0:
                fr[letter] = (current, current + 1)
            elif current == 1:
                fr[letter] = (current, current)
            else:
                fr[letter] = (current, current + 1)
            refTable[letter] = counter
            counter += 1
        # print("fr now is",fr)
        current += 1

    # print("SA before whipe ", SA)

    # Saving only the  even  number
    for idx, val in enumerate(SA):
        if val % 2 == 0:
            new_SA[idx] = val

    # print(len(bw))
    # Tally = Tally_matrix(Tally,bw,refTable)

    return FM_index(fr, bw, new_SA, refTable)  # return FM index object

# Build  FM index for all 6  possibilities
FM_List =  []
for j in range(len(translated)):
    FM_List.append(buildFM(translated[j]))

times = [0, 0, 0]
count = 0
for key, protein in proteins.items():
    count += 1
        
    start = timeit.default_timer()
    kmer_output = kmer_match(5, translated, protein)
    stop = timeit.default_timer()
    
    times[0] += stop - start
    
    if kmer_output:
        if kmer_output[0][1] == []: continue
    else: continue

    print(f"Aligning against protein {key}")
    
    print(f"K-mer output (hits, match index, frame number) {kmer_output}")
    print('Time taken: ', stop - start)
    
    start = timeit.default_timer()
    
    sw_output = []
    for j in range(len(translated)):
        result = alignment(protein, translated[j])
        if result: sw_output.append(result)
    
    stop = timeit.default_timer()
    
    times[1] += stop - start
    
    print(f"Smith Waterman output (match index) {sw_output}")
    print('Time taken: ', stop - start)
    
    start = timeit.default_timer()
    
    bm_output = []
    for j in range(len(translated)):
        result = boyerMoore(translated[j], protein)
        if result: bm_output.append(result)
    
    stop = timeit.default_timer()
    
    times[2] += stop - start
    
    print(f"Boyer Moore output (match index) {bm_output}")
    print('Time taken: ', stop - start)


    FM_output = []
    start = timeit.default_timer()
    for i in range(len(FM_List)):
        FM_List[i].search(protein)
        result = FM_List[i]
        if not result: FM_output.append(result)

    print(f"FM Index output (match index) {FM_output}")
    print('Time taken: ', stop - start)

print([time / count for time in times])