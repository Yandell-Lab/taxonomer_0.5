import collections

#codon table using numbers
bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))
rev_table = dict(zip(amino_acids, codons))
aa_table = [23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    		23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    		23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    		23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    		23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
    		14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
    		23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
    		14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23]
rev_aa_table = [23]*128
for a in rev_table:
	rev_aa_table[aa_table[ord(a)]] = ord(a)
print rev_aa_table

codon_num_table = {} #{codon}->AA number
output_relationships = []
output_let_relationships = []
output_let_num = []
for cdn in codon_table:
	codon_num_table[cdn] = aa_table[ord(codon_table[cdn])]
	output_relationships.append( "{\"%s\",%d}"%(cdn, aa_table[ord(codon_table[cdn])]) )
	output_let_relationships.append( "{\"%s\",\'%s\'}"%(cdn, codon_table[cdn]) )
	output_let_num.append( "{%d,\'%s\'}"%(aa_table[ord(codon_table[cdn])], codon_table[cdn]) )

print ",".join(output_relationships)
print ",".join(output_let_relationships)
print ",".join(output_let_num)

nt_table = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]

#reverse complement lookup for nucleotides
rev_comp = ['N']*128
pairing = {}
complements = {'A':'T','T':'A','C':'G','G':'C','a':'t','t':'a','c':'g','g':'c'}
for c in complements:
	rev_comp[ord(c)] = complements[c]
	rev_comp[ord(c.lower())] = complements[c.lower()]

print rev_comp


print rev_comp[ord('C')]
print rev_comp[ord('a')]

#look at collapsed alphabet

#collapsed_alphabet = {'S':'S','G':'S','V':'V','M':'V','L':'V','I':'V','C':'V','F':'F','Y':'F','N':'N','E':'N','Z':'N','B':'N','D':'N','K':'N','Q':'N','R':'N','T':'T','*':'*','A':'A','W':'W','X':'X','P':'P','H':'H'}
collapsed_alphabet = {'A': 'A', 'C': 'C', 'B': 'B', 'E': 'E', 'D': 'B', 'G': 'G', 'F': 'F', 'I': 'I', 'H': 'H', 'K': 'K', 'M': 'M', 'L': 'I', 'N': 'N', 'Q': 'E', 'P': 'P', 'S': 'A', 'R': 'K', 'T': 'A', 'W': 'W', 'V': 'I', 'Y': 'Y', 'X': 'A', 'Z': 'E', '*':'*'} #kmeans 14

aa_lets = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','*']

collapsed_table = [23]*128
for i,a in enumerate(aa_lets):
	collapsed_table[ord(a)] = aa_table[ord(collapsed_alphabet[a])]
	collapsed_table[ord(a.lower())] = aa_table[ord(collapsed_alphabet[a].lower())]

print collections.Counter(aa_table)
print len(collections.Counter(aa_table))
print collections.Counter(collapsed_table)
print len(collections.Counter(collapsed_table))

print aa_table
print collapsed_table


output_relationships = []
for cdn in codon_table:
	output_relationships.append( "{\"%s\",%d}"%(cdn, collapsed_table[ord(codon_table[cdn])]) )

print ",".join(output_relationships)
