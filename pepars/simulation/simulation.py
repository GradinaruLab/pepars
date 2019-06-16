import numpy

from ..utils import DNA


def generate_random_sequences(
        template_sequence="NNKNNKNNKNNKNNKNNKNNK",
        num_sequences=1000000,
        allow_stop_codon=False):

    sequences = set()

    while len(sequences) < num_sequences:

        nucleotide_sequence = ""

        for template_nucleotide in template_sequence:
            possible_nucleotides = list(DNA.IUPAC_GRAMMAR_MAP[template_nucleotide])

            nucleotide_sequence += numpy.random.choice(possible_nucleotides)

        amino_acid_sequence = DNA.translate_DNA_to_AA(nucleotide_sequence)

        if not allow_stop_codon and amino_acid_sequence.find("#") != -1:
            continue

        sequences.add(amino_acid_sequence)

    sequences = list(sequences)

    return sequences


def generate_simulated_dataset(
        template_sequence="NNKNNKNNKNNKNNKNNKNNK",
        num_sequences=1000000,
        starting_sequencing_depth=1000000,
        recovered_sequencing_depth=1000000,
        default_transduction_probability=0.01,
        default_production_rate=100,
        mean_num_viruses_starting=100,
        simulation_type="production"):

    sequences = []
    sequence_true_starting_counts = numpy.zeros((num_sequences, ))
    sequence_true_recovered_counts = numpy.zeros((num_sequences, ))

    default_transduction_x = -numpy.log((1-default_transduction_probability)/
                                        default_transduction_probability)

    sequence_index = 0

    while len(sequences) < num_sequences:

        if len(sequences) % 100000 == 0:
            print("%i sequence(s) generated so far" % len(sequences))

        nucleotide_sequence = ""

        for template_nucleotide in template_sequence:

            possible_nucleotides = list(DNA.IUPAC_GRAMMAR_MAP[template_nucleotide])

            nucleotide_sequence += numpy.random.choice(possible_nucleotides)

        amino_acid_sequence = DNA.translate_DNA_to_AA(nucleotide_sequence)

        if amino_acid_sequence.find("#") != -1:
            continue

        # Pick the original sequence count; that is, the number of virus
        # sequences that get made, by the Poisson. The scale doesn't matter here
        # since sampling will be relative, but perhaps the shape matters
        sequence_count = 1

        sequences.append(nucleotide_sequence)
        sequence_true_starting_counts[sequence_index] = sequence_count

        if simulation_type == "transduction":

            sequence_transduction_z_score = numpy.random.normal(
                loc=default_transduction_x)

            if amino_acid_sequence[0] == "T":
                sequence_transduction_z_score += 5

            sequence_transduction_probability = 1/(1+numpy.exp(-sequence_transduction_z_score))

            sequence_true_recovered_counts[sequence_index] = \
                numpy.random.binomial(
                    sequence_count,
                    p=sequence_transduction_probability
                )

        elif simulation_type == "production":

            mean_rate = default_production_rate

            if (amino_acid_sequence[0] == "T" or amino_acid_sequence[0] == "S")\
                    and (amino_acid_sequence[1] == "L" or amino_acid_sequence[1] == "I"):
                mean_rate *= 5
            elif (amino_acid_sequence[1] == "T" or amino_acid_sequence[1] == "S")\
                    and (amino_acid_sequence[2] == "L" or amino_acid_sequence[2] == "I"):
                mean_rate *= 5
            elif (amino_acid_sequence[2] == "T" or amino_acid_sequence[2] == "S")\
                    and (amino_acid_sequence[3] == "L" or amino_acid_sequence[3] == "I"):
                mean_rate *= 5

            if amino_acid_sequence[6] == "Q" and amino_acid_sequence[5] != "H" \
                    and amino_acid_sequence[5] != "S" and \
                    amino_acid_sequence[4] != "Y":
                mean_rate *= 5

            if amino_acid_sequence[2:6].find("P") != -1 and (
                    amino_acid_sequence[2:6].find("F") > amino_acid_sequence[2:6].find("P") or
                    amino_acid_sequence[2:6].find("G") > amino_acid_sequence[2:6].find("P") or
                    amino_acid_sequence[2:6].find("V") > amino_acid_sequence[2:6].find("P")):
                mean_rate *= 5

            # if amino_acid_sequence[1] == "T":
            #     mean_rate *= 2
            # if amino_acid_sequence[4] == "N":
            #     mean_rate *= 7

            sequence_production_rates = numpy.random.poisson(
                mean_rate,
                size=sequence_count)

            sequence_true_recovered_counts[sequence_index] = \
                sum(sequence_production_rates)

        sequence_index += 1

    sequence_starting_probabilities = \
        sequence_true_starting_counts / sequence_true_starting_counts.sum()

    sequence_recovered_probabilities = \
        sequence_true_recovered_counts / sequence_true_recovered_counts.sum()

    sequence_start_reads = numpy.random.choice(
        sequences,
        size=starting_sequencing_depth,
        p=sequence_starting_probabilities)

    sequence_recovered_reads = numpy.random.choice(
        sequences,
        size=recovered_sequencing_depth,
        p=sequence_recovered_probabilities)

    sequence_starting_counts = {}

    for sequence_start_read in sequence_start_reads:

        read_sequence = ""
        for nucleotide in sequence_start_read:
            if numpy.random.random() < 0.001:
                read_nucleotide = numpy.random.choice(DNA.get_nucleotides())
                read_sequence += read_nucleotide
            else:
                read_sequence += nucleotide

        amino_acid_sequence = DNA.translate_DNA_to_AA(read_sequence)

        if amino_acid_sequence.find("#") != -1:
            continue

        if amino_acid_sequence not in sequence_starting_counts:
            sequence_starting_counts[amino_acid_sequence] = 0
        sequence_starting_counts[amino_acid_sequence] += 1

    sequence_recovered_counts = {}

    for sequence_recovered_read in sequence_recovered_reads:

        read_sequence = ""
        for nucleotide in sequence_recovered_read:
            if numpy.random.random() < 0.001:
                read_nucleotide = numpy.random.choice(DNA.get_nucleotides())
                read_sequence += read_nucleotide
            else:
                read_sequence += nucleotide

        amino_acid_sequence = DNA.translate_DNA_to_AA(read_sequence)

        if amino_acid_sequence.find("#") != -1:
            continue

        if amino_acid_sequence not in sequence_recovered_counts:
            sequence_recovered_counts[amino_acid_sequence] = 0
        sequence_recovered_counts[amino_acid_sequence] += 1

    return sequence_starting_counts, sequence_recovered_counts
