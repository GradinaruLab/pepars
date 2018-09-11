import numpy

# Function to convert nucleotide to amino acid sequence

# Function to get list of amino acid characteristics for a specific amino acid

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'#', 'TAG':'#', 'TGC':'C', 'TGT':'C', 'TGA':'#', 'TGG':'W'}

IUPAC = {
    'A':'A',
    'C':'C',
    'G':'G',
    'T':'T',
    'M':'AC',
    'R':'AG',
    'W':'AT',
    'S':'CG',
    'Y':'CT',
    'K':'GT',
    'V':'ACG',
    'H':'ACT',
    'D':'AGT',
    'B':'CGT',
    'N':'ACGT'
}

DNA_complements = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}

def get_nucleotides():
    return ['A','C','G','T']

def get_amino_acids():
    return ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
 
# a function to translate a dna sequence in a single frame
def translate_dna_single(dna):

    amino_acids = ''
    for start_index in range(0, len(dna) - 2, 3):
        amino_acids += gencode[dna[start_index:start_index+3]]

    return amino_acids


def translate_amino_acid_to_random_codons(amino_acid_sequence, template="NNN"):

    nucleotide_sequence = ""

    for amino_acid_to_generate in amino_acid_sequence:

        possible_codons = []

        for codon in gencode:

            if amino_acid_to_generate != gencode[codon]:
                continue

            is_valid_codon = True

            for nucleotide_index, nucleotide in enumerate(codon):

                if nucleotide not in IUPAC[template[nucleotide_index]]:
                    is_valid_codon = False
                    break

            if not is_valid_codon:
                continue

            possible_codons.append(codon)

        if len(possible_codons) == 0:
            raise ValueError("Template does not allow for this amino acid sequence!")

        codon = numpy.random.choice(possible_codons)

        nucleotide_sequence += codon

    return nucleotide_sequence


def get_complement(sequence):

    complement_sequence = ''

    for character in sequence:
        complement_sequence += DNA_complements[character]

    return complement_sequence

def get_reverse_complement(sequence):

    reverse_sequence = sequence[::-1]

    reverse_complement_sequence = get_complement(reverse_sequence)

    return reverse_complement_sequence

def get_sequence_distance(sequence_1, sequence_2):

    min_length = min(len(sequence_1), len(sequence_2))
    max_length = max(len(sequence_1), len(sequence_2))

    distance = 0

    for character_index in range(0, min_length):
        if sequence_1[character_index] != sequence_2[character_index]:
            distance += 1

    distance += max_length - min_length

    return distance
