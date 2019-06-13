import itertools

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

# Codon frequencies from http://hive.biochemistry.gwu.edu/review/codon
CODON_FREQUENCIES = {
    "mus_musculus": {
        "AAA": 21336462,
        "AAC": 17634886,
        "AAG": 30304539,
        "AAT": 14225820,
        "ACA": 14874852,
        "ACC": 16170180,
        "ACG": 5110332,
        "ACT": 12423691,
        "AGA": 11643187,
        "AGC": 18543196,
        "AGG": 11531432,
        "AGT": 12613419,
        "ATA": 6519548,
        "ATC": 18116609,
        "ATG": 19360204,
        "ATT": 13112138,
        "CAA": 11725522,
        "CAC": 13569509,
        "CAG": 32740441,
        "CAT": 10021378,
        "CCA": 16999974,
        "CCC": 16305633,
        "CCG": 5442893,
        "CCT": 17872896,
        "CGA": 6143838,
        "CGC": 7654491,
        "CGG": 9299948,
        "CGT": 4129594,
        "CTA": 7175663,
        "CTC": 16747656,
        "CTG": 33194957,
        "CTT": 11998369,
        "GAA": 27149113,
        "GAC": 23457466,
        "GAG": 36979161,
        "GAT": 19954082,
        "GCA": 14997741,
        "GCC": 22410169,
        "GCG": 5227733,
        "GCT": 18003486,
        "GGA": 14936244,
        "GGC": 17622820,
        "GGG": 13233138,
        "GGT": 9877142,
        "GTA": 6680252,
        "GTC": 12937644,
        "GTG": 23653667,
        "GTT": 9664997,
        "TAA": 346998,
        "TAC": 12845107,
        "TAG": 310192,
        "TAT": 9950306,
        "TCA": 11896395,
        "TCC": 16330278,
        "TCG": 3829599,
        "TCT": 15524659,
        "TGA": 663764,
        "TGC": 9741288,
        "TGG": 10179088,
        "TGT": 9509181,
        "TTA": 6532569,
        "TTC": 16734250,
        "TTG": 11852515,
        "TTT": 14238031
    }
}


def get_nucleotides():
    return NUCLEOTIDES


def get_amino_acids():
    return AMINO_ACIDS[:]


def get_degenerate_nucleotides():
    return set(IUPAC_GRAMMAR_MAP.keys()).difference(set(NUCLEOTIDES))


def translate_DNA_to_AA(DNA_sequence):

    amino_acids = ""

    for start_index in range(0, len(DNA_sequence) - 2, 3):
        amino_acids += CODON_AA_MAP[DNA_sequence[start_index:start_index + 3]]

    return amino_acids


def get_optimal_codon(
        amino_acid_sequence,
        template="NNN",
        species="mus_musculus"):

    nucleotide_sequence = ""

    for amino_acid in amino_acid_sequence:

        possible_codons = []

        for possible_codon in CODON_AA_MAP:

            if CODON_AA_MAP[possible_codon] != amino_acid:
                continue

            is_possible_codon = True

            for nucleotide_index, nucleotide in enumerate(possible_codon):
                if nucleotide not in \
                        IUPAC_GRAMMAR_MAP[template[nucleotide_index]]:
                    is_possible_codon = False
                    break

            if not is_possible_codon:
                continue

            possible_codons.append(possible_codon)

        optimal_codon = None
        optimal_codon_frequency = -numpy.infty

        for codon in possible_codons:
            if CODON_FREQUENCIES[species][codon] > optimal_codon_frequency:
                optimal_codon_frequency = CODON_FREQUENCIES[species][codon]
                optimal_codon = codon

        nucleotide_sequence += optimal_codon

    return nucleotide_sequence


def get_all_possible_nucleotide_seqs(amino_acid_sequence, template="NNN"):
    '''Returns all the possible nucleotide sequences of an input amino acid
    sequence

    Inputs:
        - amino_acid_sequence: The query amino acid sequence
        - template: Which template to follow to generate sequences
    Outputs:
        - all_pos_seqs: All the possible sequences for given amino_acid_sequence

    '''
    all_pos_codons = []

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

        all_pos_codons.append(possible_codons)

    all_pos_seqs = [''.join(i) for i in itertools.product(*all_pos_codons)]        
    return all_pos_seqs


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


def is_sequence_match(sequence_1, sequence_2, allowable_distance=0):

    num_mismatches = 0

    for character_index, character in enumerate(sequence_1):
        if character != sequence_2[character_index]:
            num_mismatches += 1
            if num_mismatches > allowable_distance:
                return False

    return True


def is_template_match(template, sequence, allowable_distance=0):

    if allowable_distance == 0:
        if len(template) != len(sequence):
            return False

        for character_index, character in enumerate(sequence):
            if character not in IUPAC_GRAMMAR_MAP[template[character_index]]:
                return False
    else:
        num_nucleotides_off = 0

        for character_index, character in enumerate(sequence):
            if character not in IUPAC_GRAMMAR_MAP[template[character_index]]:
                num_nucleotides_off += 1
                if num_nucleotides_off > allowable_distance:
                    return False

    return True


def get_variant_indices_from_template(template_sequence):

    template_indices = []

    degenerate_nucleotides = get_degenerate_nucleotides()

    for character_index, character in enumerate(template_sequence):
        if character in degenerate_nucleotides:
            template_indices.append(character_index)

    return template_indices


def get_template_distance(template, sequence):

    template_distance = 0

    for character_index, character in enumerate(sequence):
        if character not in IUPAC_GRAMMAR_MAP[template[character_index]]:
            template_distance += 1

    return template_distance
