import numpy

CODON_AA_MAP = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M", "ACA": "T", "ACC": "T",
    "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K", "AGC": "S", "AGT": "S",
    "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L", "CCA": "P", "CCC": "P",
    "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q", "CGA": "R", "CGC": "R",
    "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V", "GCA": "A", "GCC": "A",
    "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E", "GGA": "G", "GGC": "G",
    "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S", "TTC": "F", "TTT": "F",
    "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "#", "TAG": "#", "TGC": "C", "TGT": "C",
    "TGA": "#", "TGG": "W"}

IUPAC_GRAMMAR_MAP = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "N": "ACGT"
}

DNA_COMPLEMENT_MAP = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G"
}

NUCLEOTIDE_INDEX_MAP = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3
}

AMINO_ACID_INDEX_MAP = {
    "A": 0,
    "C": 1,
    "D": 2,
    "E": 3,
    "F": 4,
    "G": 5,
    "H": 6,
    "I": 7,
    "K": 8,
    "L": 9,
    "M": 10,
    "N": 11,
    "P": 12,
    "Q": 13,
    "R": 14,
    "S": 15,
    "T": 16,
    "V": 17,
    "W": 18,
    "Y": 19
}

NUCLEOTIDES = ["A", "C", "G", "T"]

AMINO_ACIDS = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y"
]


def get_nucleotides():
    return NUCLEOTIDES


def get_amino_acids():
    return AMINO_ACIDS


def translate_DNA_to_AA(DNA_sequence):

    amino_acids = ""

    for start_index in range(0, len(DNA_sequence) - 2, 3):
        amino_acids += CODON_AA_MAP[DNA_sequence[start_index:start_index + 3]]

    return amino_acids


def translate_AA_to_DNA(amino_acid_sequence, template="NNN"):

    nucleotide_sequence = ""

    for amino_acid_to_generate in amino_acid_sequence:

        possible_codons = []

        for codon in CODON_AA_MAP:

            if amino_acid_to_generate != CODON_AA_MAP[codon]:
                continue

            is_valid_codon = True

            for nucleotide_index, nucleotide in enumerate(codon):

                if nucleotide not in \
                        IUPAC_GRAMMAR_MAP[template[nucleotide_index]]:
                    is_valid_codon = False
                    break

            if not is_valid_codon:
                continue

            possible_codons.append(codon)

        if len(possible_codons) == 0:
            raise ValueError(
                "Template does not allow for this amino acid sequence!")

        codon = numpy.random.choice(possible_codons)

        nucleotide_sequence += codon

    return nucleotide_sequence


def translate_complement(sequence):

    complement_sequence = ""

    for character in sequence:
        complement_sequence += DNA_COMPLEMENT_MAP[character]

    return complement_sequence


def translate_reverse_complement(sequence):

    reverse_sequence = sequence[::-1]

    reverse_complement_sequence = translate_complement(reverse_sequence)

    return reverse_complement_sequence


def is_template_match(template, sequence):

    if len(template) != len(sequence):
        return False

    for character_index, character in enumerate(sequence):
        if character not in IUPAC_GRAMMAR_MAP[template[character_index]]:
            return False

    return True
