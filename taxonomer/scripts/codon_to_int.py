import collections

bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

def codon_to_int(cdn):
    ncdn = 0
    for l in cdn:
        ncdn = ncdn << 2
        if l == 'A':
            continue
        elif l == 'T':
            ncdn = ncdn | 3
            continue
        elif l == 'C':
            ncdn = ncdn | 1
            continue
        elif l == 'G':
            ncdn = ncdn | 2
            continue
        else:
            return -1
    return ncdn

#collapsed_alphabet = {'S':'S','G':'S','V':'V','M':'V','L':'V','I':'V','C':'V','F':'F','Y':'F','N':'N','E':'N','Z':'N','B':'N','D':'N','K':'N','Q':'N','R':'N','T':'T','*':'*','A':'A','W':'W','X':'X','P':'P','H':'H'} #original collapse from mark's perl script
collapsed_alphabet = {'A': 'A', 'C': 'C', 'E': 'E', 'D': 'D', 'G': 'G', 'F': 'F', 'I': 'I', 'H': 'H', 'K': 'K', 'M': 'M', 'L': 'I', 'N': 'N', 'Q': 'E', 'P': 'P', 'S': 'S', 'R': 'K', 'T': 'S', 'W': 'W', 'V': 'I', 'Y': 'F', '*':'*'} #kmeans 14

codon_vals = []
aa_map = {}
codon_mapped_vals = [0]*64
seen_aa = {} #{AA}->(codon num)
for cdn,aa in codon_table.items():
    cdn_val = codon_to_int(cdn)
    codon_vals.append(cdn_val)
    if aa in seen_aa:
        codon_mapped_vals[cdn_val] = seen_aa[aa]
    else:
        aa_map[aa] = cdn
        codon_mapped_vals[cdn_val] = cdn_val
        seen_aa[aa] = cdn_val

collapsed_codon_mapped_vals = [0]*64
for cdn,aa in codon_table.items():
    cdn_val = codon_to_int(cdn)
    collapsed_codon_mapped_vals[cdn_val] = seen_aa[collapsed_alphabet[aa]]

print aa_map
print "aa_map"
print codon_vals
print "codon_vals"
print codon_mapped_vals
print "codon_mapped_vals"
print collapsed_codon_mapped_vals
print "collapsed_codon_mapped_vals"
print seen_aa
print "seen_aa"
print collections.Counter(codon_mapped_vals)
print ""
print collections.Counter(collapsed_codon_mapped_vals)
print ""
