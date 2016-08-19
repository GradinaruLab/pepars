from Bio.SeqUtils.ProtParam import ProteinAnalysis
from analysis.Sequence_Library import Sequence_Library


def compare_to(val1, val2):
	
	if val1<val2:
		return -1
	elif val1>val2:
		return 1
	return 0

class AminoAcid:
	
	def __init__(self,amino_acid):


		self.is_valid = True
		self.molecular_weight = ProteinAnalysis(amino_acid).molecular_weight()
		self.hydrophobicity = ProteinAnalysis(amino_acid).gravy()

		if amino_acid == 'R':
			self.is_valid = True
			self.amino_acid_letter = 'R'
			self.amino_acid_name = 'Arginine'
			self.charge = 1.00
			self.hydropathy = 'hydrophilic'
			self.solubility = 71.8
			self.phosphorylation = 0.00
			self.average_flexibility_idx = 0.530
			self.ionic_bond = 1.00
			self.hydrogen_bond = 1.00


		elif amino_acid== 'N':
			self.is_valid = True
			self.amino_acid_letter = 'N'
			self.amino_acid_name = 'Asparagine'
			self.charge = 0.00
			self.hydropathy = 'hydrophilic'
			self.solubility = 2.4
			self.phosphorylation = 0.00
			self.average_flexibility_idx = 0.46
			self.ionic_bond = 0.00
			self.hydrogen_bond = 1.00

		elif amino_acid == 'D':
			self.is_valid = True
			self.amino_acid_letter = 'D'
			self.amino_acid_name = 'Aspartate'
			self.charge = -1.00
			self.hydropathy = 'hydrophilic'
			self.solubility = 0.42
			self.phosphorylation = 0.00
			self.average_flexibility_idx = 0.510
			self.ionic_bond = 1.00
			self.hydrogen_bond = 1.00

		elif amino_acid == 'E':
			self.is_valid = True
			self.amino_acid_letter = 'E'
			self.amino_acid_name = 'Glutamate'
			self.charge = -1.00
			self.hydropathy = 'hydrophilic'
			self.solubility = 0.72
			self.phosphorylation = 0.0
			self.average_flexibility_idx = 0.500
			self.ionic_bond = 1.00
			self.hydrogen_bond = 1.00

		elif amino_acid== 'Q':
			self.is_valid = True
			self.amino_acid_letter = 'Q'
			self.amino_acid_name = 'Glutamine'
			self.charge = 0.00
			self.hydropathy = 'hydrophilic'
			self.solubility = 2.6
			self.phosphorylation = 0.0
			self.average_flexibility_idx = 0.490
			self.ionic_bond = 0.00
			self.hydrogen_bond = 1.00

		elif amino_acid== 'K':
			self.is_valid = True
			self.amino_acid_letter = 'K'
			self.amino_acid_name = 'Lysine'
			self.charge = 1.00
			self.hydropathy = 'hydrophilic'
			self.solubility = 0.0
			self.phosphorylation = 0.0
			self.average_flexibility_idx = 0.470
			self.ionic_bond = 1.00
			self.hydrogen_bond = 1.00

		elif amino_acid== 'S':
			self.is_valid = True
			self.amino_acid_letter = 'S'
			self.amino_acid_name = 'Serine'
			self.charge = -1.00
			self.hydropathy = 'hydrophilic'
			self.solubility = 36.2
			self.phosphorylation = 1.00
			self.average_flexibility_idx = 0.510
			self.ionic_bond = 0.00
			self.hydrogen_bond = 1.00

		elif amino_acid== 'T':
			self.is_valid = True
			self.amino_acid_letter = 'T'
			self.amino_acid_name = 'Threonine'
			self.charge = -1.00
			self.hydropathy = 'hydrophilic'
			self.solubility = 0.0
			self.phosphorylation = 1.00
			self.average_flexibility_idx = 0.440
			self.ionic_bond = 0.00
			self.hydrogen_bond = 1.00

		elif amino_acid== 'C':
			self.is_valid = True
			self.amino_acid_letter = 'C'
			self.amino_acid_name = 'Cysteine'
			self.charge = -1.00
			self.hydropathy = 'moderate'
			self.solubility = 0.0
			self.phosphorylation = 0.0
			self.average_flexibility_idx = 0.350
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		elif amino_acid== 'H':
			self.is_valid = True
			self.amino_acid_letter = 'H'
			self.amino_acid_name = 'Histidine'
			self.charge = 1.00
			self.hydropathy = 'moderate'
			self.solubility = 4.19
			self.phosphorylation = 0.0
			self.average_flexibility_idx = 0.320
			self.ionic_bond = 1.00
			self.hydrogen_bond = 1.00

		elif amino_acid== 'M':
			self.is_valid = True
			self.amino_acid_letter = 'M'
			self.amino_acid_name = 'Methionine'
			self.charge = 0.00
			self.hydropathy = 'moderate'
			self.solubility = 5.14
			self.phosphorylation = 0.0
			self.average_flexibility_idx = 0.300
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00


		elif amino_acid== 'A':
			self.is_valid = True
			self.amino_acid_letter = 'A'
			self.amino_acid_name = 'Alanine'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 15.8
			self.phosphorylation = 0.0
			self.average_flexibility_idx = 0.360
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		elif amino_acid== 'V':
			self.is_valid = True
			self.amino_acid_letter = 'V'
			self.amino_acid_name = 'Valine'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 5.6
			self.phosphorylation = 0.0
			self.average_flexibility_idx = 0.390
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		elif amino_acid== 'G':
			self.is_valid = True
			self.amino_acid_letter = 'G'
			self.amino_acid_name = 'Glycine'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 22.5
			self.phosphorylation = 0.0
			self.average_flexibility_idx = 0.540
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		elif amino_acid== 'I':
			self.is_valid = True
			self.amino_acid_letter = 'I'
			self.amino_acid_name = 'Isoleucine'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 3.36
			self.phosphorylation = 0.0
			self.average_flexibility_idx = 0.460
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		elif amino_acid== 'L':
			self.is_valid = True
			self.amino_acid_letter = 'L'
			self.amino_acid_name = 'Leucine'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 2.37
			self.phosphorylation = 0.0
			self.average_flexibility_idx = 0.370
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		elif amino_acid== 'F':
			self.is_valid = True
			self.amino_acid_letter = 'F'
			self.amino_acid_name = 'Phenylalanine'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 2.7
			self.phosphorylation = 0.0
			self.average_flexibility_idx = 0.310
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		elif amino_acid=='P':
			self.is_valid = True
			self.amino_acid_letter = 'P'
			self.amino_acid_name = 'Proline'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 1.54
			self.phosphorylation = 0.0
			self.average_flexibility_idx = 0.510
			self.ionic_bond = 0.00
			self.hydrogen_bond = 0.00

		elif amino_acid== 'W':
			self.is_valid = True
			self.amino_acid_letter = 'W'
			self.amino_acid_name = 'Tryptophan'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 1.00
			self.phosphorylation = 0.0
			self.average_flexibility_idx = 0.310
			self.ionic_bond = 0.00
			self.hydrogen_bond = 1.00

		elif amino_acid== 'Y':
			self.is_valid = True
			self.amino_acid_letter = 'Y'
			self.amino_acid_name = 'Tyrosine'
			self.charge = -1.00
			self.hydropathy = 'hydrophobic'
			self.solubility = 0.038
			self.phosphorylation = 1.00
			self.average_flexibility_idx = 0.420
			self.ionic_bond = 0.00
			self.hydrogen_bond = 1.00

		else:
			self.is_valid = False
			print "Invalid Amino Acid "+amino_acid

	

	def compare_features(self, acid, feature):
		if feature == 'charge':
			return compare_to(self.charge, acid.charge)
		if feature == 'hydropathy':
			return compare_to(self.hydropathy, acid.hydropathy)
		if feature == 'solubility':
			return compare_to(self.solubility, acid.solubility)
		if feature == 'phosphorylation':
			return compare_to(self.phosphorylation, acid.phosphorylation)
		if feature == 'average_flexibility_idx':
			return compare_to(self.average_flexibility_idx, acid.average_flexibility_idx)
		if feature == 'hydrophobicity':
			return compare_to(self.hydrophobicity, acid.hydrophobicity)
		if feature == 'ionic_bond':
			return compare_to(self.ionic_bond, acid.ionic_bond)
		if feature == 'hydrogen_bond':
			return compare_to(self.hydrogen_bond, acid.hydrogen_bond)
		raise Exception('Feature {} not defined'.format(feature))


def sort_by_features(acids, features, rules=1):
	order=range(len(acids))
	for acid in acids:
		if not acid.is_valid:
			del acid
	if not type(rules).__name__=='list':
		rules= [rules for i in range(len(features))]
	for rule in rules:
		assert(rule in [-1, 1])
	assert(len(rules)==len(features))
	n = len(acids)
	for feature, rule in zip(features[::-1], rules[::-1]):
		#insertion sort for each feature
		for j in range(1, n):
			i = j-1
			key = acids[j]
			key2 = order[j]
			while i>=0 and rule*key.compare_features(acids[i],feature)>0:
				#print key.hydrophobicity,'is bigger than',acids[i].hydrophobicity
				acids[i+1] = acids[i]
				order[i+1] = order[i]
				#print 'Swapping',order[i],order[i+1]
				i -= 1
			acids[i+1] = key
			order[i+1] = key2
	return order