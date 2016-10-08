from Bio.SeqUtils.ProtParam import ProteinAnalysis
from analysis.Sequence_Library import Sequence_Library


def compare_to(val1, val2):
	
	if val1<val2:
		return -1
	elif val1>val2:
		return 1
	return 0

class AminoAcid:

	properties = ['charge','hydropathy','solubility','phosphorylation','average_flexibility_idx','ionic_bond','molecular_weight','hydrophobicity','typically_helix','typically_turn','typically_sheet']

	
	def __init__(self,amino_acid):


	self.properties = {}

		self.is_valid = True



		if amino_acid == 'R':
			self.is_valid = True
			self.amino_acid_letter = 'R'
			self.amino_acid_name = 'Arginine'
			self.properties['charge'] = 1.00
			self.properties['hydropathy'] = 1
			self.properties['solubility'] = 71.8
			self.properties['phosphorylation'] = 0.00
			self.properties['average_flexibility_idx'] = 0.530
			self.properties['ionic_bond'] = 1.00
			self.properties['hydrogen_bond'] = 1.00


		elif amino_acid== 'N':
			self.is_valid = True
			self.amino_acid_letter = 'N'
			self.amino_acid_name = 'Asparagine'
			self.properties['charge'] = 0.00
			self.properties['hydropathy'] = 1
			self.properties['solubility'] = 2.4
			self.properties['phosphorylation'] = 0.00
			self.properties['average_flexibility_idx'] = 0.46
			self.properties['ionic_bond'] = 0.00
			self.properties['hydrogen_bond'] = 1.00

		elif amino_acid == 'D':
			self.is_valid = True
			self.amino_acid_letter = 'D'
			self.amino_acid_name = 'Aspartate'
			self.properties['charge'] = -1.00
			self.properties['hydropathy'] = 1
			self.properties['solubility'] = 0.42
			self.properties['phosphorylation'] = 0.00
			self.properties['average_flexibility_idx'] = 0.510
			self.properties['ionic_bond'] = 1.00
			self.properties['hydrogen_bond'] = 1.00

		elif amino_acid == 'E':
			self.is_valid = True
			self.amino_acid_letter = 'E'
			self.amino_acid_name = 'Glutamate'
			self.properties['charge'] = -1.00
			self.properties['hydropathy'] = 1
			self.properties['solubility'] = 0.72
			self.properties['phosphorylation'] = 0.0
			self.properties['average_flexibility_idx'] = 0.500
			self.properties['ionic_bond'] = 1.00
			self.properties['hydrogen_bond'] = 1.00

		elif amino_acid== 'Q':
			self.is_valid = True
			self.amino_acid_letter = 'Q'
			self.amino_acid_name = 'Glutamine'
			self.properties['charge'] = 0.00
			self.properties['hydropathy'] = 1
			self.properties['solubility'] = 2.6
			self.properties['phosphorylation'] = 0.0
			self.properties['average_flexibility_idx'] = 0.490
			self.properties['ionic_bond'] = 0.00
			self.properties['hydrogen_bond'] = 1.00

		elif amino_acid== 'K':
			self.is_valid = True
			self.amino_acid_letter = 'K'
			self.amino_acid_name = 'Lysine'
			self.properties['charge'] = 1.00
			self.properties['hydropathy'] = 1
			self.properties['solubility'] = 0.0
			self.properties['phosphorylation'] = 0.0
			self.properties['average_flexibility_idx'] = 0.470
			self.properties['ionic_bond'] = 1.00
			self.properties['hydrogen_bond'] = 1.00

		elif amino_acid== 'S':
			self.is_valid = True
			self.amino_acid_letter = 'S'
			self.amino_acid_name = 'Serine'
			self.properties['charge'] = -1.00
			self.properties['hydropathy'] = 1
			self.properties['solubility'] = 36.2
			self.properties['phosphorylation'] = 1.00
			self.properties['average_flexibility_idx'] = 0.510
			self.properties['ionic_bond'] = 0.00
			self.properties['hydrogen_bond'] = 1.00

		elif amino_acid== 'T':
			self.is_valid = True
			self.amino_acid_letter = 'T'
			self.amino_acid_name = 'Threonine'
			self.properties['charge'] = -1.00
			self.properties['hydropathy'] = 1
			self.properties['solubility'] = 0.0
			self.properties['phosphorylation'] = 1.00
			self.properties['average_flexibility_idx'] = 0.440
			self.properties['ionic_bond'] = 0.00
			self.properties['hydrogen_bond'] = 1.00

		elif amino_acid== 'C':
			self.is_valid = True
			self.amino_acid_letter = 'C'
			self.amino_acid_name = 'Cysteine'
			self.properties['charge'] = -1.00
			self.properties['hydropathy'] = 0
			self.properties['solubility'] = 0.0
			self.properties['phosphorylation'] = 0.0
			self.properties['average_flexibility_idx'] = 0.350
			self.properties['ionic_bond'] = 0.00
			self.properties['hydrogen_bond'] = 0.00

		elif amino_acid== 'H':
			self.is_valid = True
			self.amino_acid_letter = 'H'
			self.amino_acid_name = 'Histidine'
			self.properties['charge'] = 1.00
			self.properties['hydropathy'] = 0
			self.properties['solubility'] = 4.19
			self.properties['phosphorylation'] = 0.0
			self.properties['average_flexibility_idx'] = 0.320
			self.properties['ionic_bond'] = 1.00
			self.properties['hydrogen_bond'] = 1.00

		elif amino_acid== 'M':
			self.is_valid = True
			self.amino_acid_letter = 'M'
			self.amino_acid_name = 'Methionine'
			self.properties['charge'] = 0.00
			self.properties['hydropathy'] = 0
			self.properties['solubility'] = 5.14
			self.properties['phosphorylation'] = 0.0
			self.properties['average_flexibility_idx'] = 0.300
			self.properties['ionic_bond'] = 0.00
			self.properties['hydrogen_bond'] = 0.00


		elif amino_acid== 'A':
			self.is_valid = True
			self.amino_acid_letter = 'A'
			self.amino_acid_name = 'Alanine'
			self.properties['charge'] = -1.00
			self.properties['hydropathy'] = -1
			self.properties['solubility'] = 15.8
			self.properties['phosphorylation'] = 0.0
			self.properties['average_flexibility_idx'] = 0.360
			self.properties['ionic_bond'] = 0.00
			self.properties['hydrogen_bond'] = 0.00

		elif amino_acid== 'V':
			self.is_valid = True
			self.amino_acid_letter = 'V'
			self.amino_acid_name = 'Valine'
			self.properties['charge'] = -1.00
			self.properties['hydropathy'] = -1
			self.properties['solubility'] = 5.6
			self.properties['phosphorylation'] = 0.0
			self.properties['average_flexibility_idx'] = 0.390
			self.properties['ionic_bond'] = 0.00
			self.properties['hydrogen_bond'] = 0.00

		elif amino_acid== 'G':
			self.is_valid = True
			self.amino_acid_letter = 'G'
			self.amino_acid_name = 'Glycine'
			self.properties['charge'] = -1.00
			self.properties['hydropathy'] = -1
			self.properties['solubility'] = 22.5
			self.properties['phosphorylation'] = 0.0
			self.properties['average_flexibility_idx'] = 0.540
			self.properties['ionic_bond'] = 0.00
			self.properties['hydrogen_bond'] = 0.00

		elif amino_acid== 'I':
			self.is_valid = True
			self.amino_acid_letter = 'I'
			self.amino_acid_name = 'Isoleucine'
			self.properties['charge'] = -1.00
			self.properties['hydropathy'] = -1
			self.properties['solubility'] = 3.36
			self.properties['phosphorylation'] = 0.0
			self.properties['average_flexibility_idx'] = 0.460
			self.properties['ionic_bond'] = 0.00
			self.properties['hydrogen_bond'] = 0.00

		elif amino_acid== 'L':
			self.is_valid = True
			self.amino_acid_letter = 'L'
			self.amino_acid_name = 'Leucine'
			self.properties['charge'] = -1.00
			self.properties['hydropathy'] = -1
			self.properties['solubility'] = 2.37
			self.properties['phosphorylation'] = 0.0
			self.properties['average_flexibility_idx'] = 0.370
			self.properties['ionic_bond'] = 0.00
			self.properties['hydrogen_bond'] = 0.00

		elif amino_acid== 'F':
			self.is_valid = True
			self.amino_acid_letter = 'F'
			self.amino_acid_name = 'Phenylalanine'
			self.properties['charge'] = -1.00
			self.properties['hydropathy'] = -1
			self.properties['solubility'] = 2.7
			self.properties['phosphorylation'] = 0.0
			self.properties['average_flexibility_idx'] = 0.310
			self.properties['ionic_bond'] = 0.00
			self.properties['hydrogen_bond'] = 0.00

		elif amino_acid=='P':
			self.is_valid = True
			self.amino_acid_letter = 'P'
			self.amino_acid_name = 'Proline'
			self.properties['charge'] = -1.00
			self.properties['hydropathy'] = -1
			self.properties['solubility'] = 1.54
			self.properties['phosphorylation'] = 0.0
			self.properties['average_flexibility_idx'] = 0.510
			self.properties['ionic_bond'] = 0.00
			self.properties['hydrogen_bond'] = 0.00

		elif amino_acid== 'W':
			self.is_valid = True
			self.amino_acid_letter = 'W'
			self.amino_acid_name = 'Tryptophan'
			self.properties['charge'] = -1.00
			self.properties['hydropathy'] = -1
			self.properties['solubility'] = 1.00
			self.properties['phosphorylation'] = 0.0
			self.properties['average_flexibility_idx'] = 0.310
			self.properties['ionic_bond'] = 0.00
			self.properties['hydrogen_bond'] = 1.00

		elif amino_acid== 'Y':
			self.is_valid = True
			self.amino_acid_letter = 'Y'
			self.amino_acid_name = 'Tyrosine'
			self.properties['charge'] = -1.00
			self.properties['hydropathy'] = -1
			self.properties['solubility'] = 0.038
			self.properties['phosphorylation'] = 1.00
			self.properties['average_flexibility_idx'] = 0.420
			self.properties['ionic_bond'] = 0.00
			self.properties['hydrogen_bond'] = 1.00

		else:
			self.is_valid = False
			print "Invalid Amino Acid "+amino_acid

		if (self.is_valid):
			self.properties['molecular_weight'] = ProteinAnalysis(str(amino_acid)).molecular_weight()
			self.properties['hydrophobicity'] = ProteinAnalysis(str(amino_acid)).gravy()
			secondary_struct = ProteinAnalysis(str(amino_acid)).secondary_structure_fraction()
			self.properties['typically_helix'] = secondary_struct[0]
			self.properties['typically_turn'] = secondary_struct[1]
			self.properties['typically_sheet'] = secondary_struct[2]


	

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