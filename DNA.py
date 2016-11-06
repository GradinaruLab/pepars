
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