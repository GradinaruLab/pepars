
class AminoAcid:

	
	def __init__(self,amino_acid):
		 
		if amino_acid == 'R':
			self.amino_acid_letter = 'R'
			self.amino_acid_name = 'Arginine'
			self.charge = 1.00
			self.hydropathy = 'hydrophilic'
			self.solubility = 71.8
			self.phosphorylation = 0.00
			self.average_flexibility_idx = 0.530
			self.hydrophobicity = 0.00
			self.ionic_bond = 1.00
			self.hydrogen_bond = 1.00


		if amino_acid== 'N':
			self.amino_acid_letter = 'N'
			self.amino_acid_name = 'Asparagine'
			self.charge = 0.00
			self.hydropathy = 'hydrophilic'
			self.solubility = 2.4
			self.phosphorylation = 0.00
			self.average_flexibility_idx = 0.46
			self.hydrophobicity = 0.448
			self.ionic_bond = 0.00
			self.hydrogen_bond = 1.00

		if amino_acid == 'D':
			self.amino_acid_letter = 'D'
			self.amino_acid_name = 'Aspartate'
			self.charge = -1.00
			self.hydropathy = 'hydrophilic'
			self.solubility = 0.42
			self.phosphorylation = 0.00
			self.average_flexibility_idx = 0.510
			self.hydrophobicity = 0.417
			self.ionic_bond = 1.00
			self.hydrogen_bond = 1.00

		if amino_acid == 'E':
			self.amino_acid_letter = 'E'
			self.amino_acid_name = 'Glutamate'
			self.charge = -1.00
			self.hydropathy = 'hydrophilic'
			self.solubility = 0.72
			self.phosphorylation = 'none'
			self.average_flexibility_idx = 0.500
			self.hydrophobicity = 0.458
			self.ionic_bond = 1.00
			self.hydrogen_bond = 1.00

		if amino_acid== 'Q':
			self.amino_acid_letter = 'Q'
			self.amino_acid_name = 'Glutamine'
			self.charge = 0.00
			self.hydropathy = 'hydrophilic'
			self.solubility = 2.6
			self.phosphorylation = 'none'
			self.average_flexibility_idx = 0.490
			self.hydrophobicity = 0.430
			self.ionic_bond = 0.00
			self.hydrogen_bond = 1.00

		if amino_acid== 'K':
			self.amino_acid_letter = 'K'
			self.amino_acid_name = 'Lysine'
			self.charge = 1.00
			self.hydropathy = 'hydrophilic'
			self.solubility = ''
			self.phosphorylation = 'none'
			self.average_flexibility_idx = 0.470
			self.hydrophobicity = 0.263
			self.ionic_bond = 1.00
			self.hydrogen_bond = 1.00

		if amino_acid== 'S':
			self.amino_acid_letter = 'S'
			self.amino_acid_name = 'Serine'
			self.charge = -1.00
			self.hydropathy = 'hydrophilic'
			self.solubility = 36.2
			self.phosphorylation = 1.00
			self.average_flexibility_idx = 0.510
			self.hydrophobicity = 0.601
			self.ionic_bond = 0.00
			self.hydrogen_bond = 1.00

		if amino_acid== 'T':
			self.amino_acid_letter = 'T'
			self.amino_acid_name = 'Threonine'
			self.charge = -1.00
			self.hydropathy = 'hydrophilic'
			self.solubility = 'free'
			self.phosphorylation = 1.00
			self.average_flexibility_idx = 0.440
			self.hydrophobicity = 0.634
			self.ionic_bond = 0.00
			self.hydrogen_bond = 1.00

		if amino_acid== 'C':
			self.amino_acid_letter = 'C'
			self.amino_acid_name = 'Cysteine'
			self.charge = -1.00
			self.hydropathy = 'moderate'
			self.solubility = 'free'
			self.phosphorylation = 'none'
			self.average_flexibility_idx = 0.350
			self.hydrophobicity = 0.721
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		if amino_acid== 'H':
			self.amino_acid_letter = 'H'
			self.amino_acid_name = 'Histidine'
			self.charge = 1.00
			self.hydropathy = 'moderate'
			self.solubility = 4.19
			self.phosphorylation = 'none'
			self.average_flexibility_idx = 0.320
			self.hydrophobicity = 0.548
			self.ionic_bond = 1.00
			self.hydrogen_bond = 1.00

		if amino_acid== 'M':
			self.amino_acid_letter = 'M'
			self.amino_acid_name = 'Methionine'
			self.charge = 0.00
			self.hydropathy = 'moderate'
			self.solubility = 5.14
			self.phosphorylation = 'none'
			self.average_flexibility_idx = 0.300
			self.hydrophobicity = 0.811
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00


		if amino_acid== 'A':
			self.amino_acid_letter = 'A'
			self.amino_acid_name = 'Alanine'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 15.8
			self.phosphorylation = 'none'
			self.average_flexibility_idx = 0.360
			self.hydrophobicity = 0.806
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		if amino_acid== 'V':
			self.amino_acid_letter = 'V'
			self.amino_acid_name = 'Valine'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 5.6
			self.phosphorylation = 'none'
			self.average_flexibility_idx = 0.390
			self.hydrophobicity = 0.923
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		if amino_acid== 'G':
			self.amino_acid_letter = 'G'
			self.amino_acid_name = 'Glycine'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 22.5
			self.phosphorylation = 'none'
			self.average_flexibility_idx = 0.540
			self.hydrophobicity = 0.770
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		if amino_acid== 'I':
			self.amino_acid_letter = 'I'
			self.amino_acid_name = 'Isoleucine'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 3.36
			self.phosphorylation = 'none'
			self.average_flexibility_idx = 0.460
			self.hydrophobicity = 1.00
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		if amino_acid== 'L':
			self.amino_acid_letter = 'L'
			self.amino_acid_name = 'Leucine'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 2.37
			self.phosphorylation = 'none'
			self.average_flexibility_idx = 0.370
			self.hydrophobicity = 0.918
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		if amino_acid== 'F':
			self.amino_acid_letter = 'F'
			self.amino_acid_name = 'Phenylalanine'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 2.7
			self.phosphorylation = 'none'
			self.average_flexibility_idx = 0.310
			self.hydrophobicity = 0.951
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		if amino_acid=='P':
			self.amino_acid_letter = 'P'
			self.amino_acid_name = 'Proline'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 1.54
			self.phosphorylation = 'none'
			self.average_flexibility_idx = 0.510
			self.hydrophobicity = 0.678
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		if amino_acid== 'W':
			self.amino_acid_letter = 'W'
			self.amino_acid_name = 'Tryptophan'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 1.00
			self.phosphorylation = 'none'
			self.average_flexibility_idx = 0.310
			self.hydrophobicity = 0.854
			self.ionic_bond = 0.00
			self.hydrogen_bond = 1.00

		if amino_acid== 'Y':
			self.amino_acid_letter = 'Y'
			self.amino_acid_name = 'Tyrosine'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 0.038
			self.phosphorylation = 1.00
			self.average_flexibility_idx = 0.420
			self.hydrophobicity = 0.714
			self.ionic_bond = 0.00
			self.hydrogen_bond = 1.00

	





		








