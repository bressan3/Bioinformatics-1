""" Codons Dictionary """

code = {}

code['UUU'] = 'F'
code['UUC'] = 'F'
code['UUA'] = 'L'
code['UUG'] = 'L'

code['UCU'] = 'S'
code['UCC'] = 'S'
code['UCA'] = 'S'
code['UCG'] = 'S'

code['UAU'] = 'Y'
code['UAC'] = 'Y'
code['UAA'] = '*'
code['UAG'] = '*'

code['UGU'] = 'C'
code['UGC'] = 'C'
code['UGA'] = '*'
code['UGG'] = 'W'

code['CUU'] = 'L'
code['CUC'] = 'L'
code['CUA'] = 'L'
code['CUG'] = 'L'

code['CCU'] = 'P'
code['CCC'] = 'P'
code['CCA'] = 'P'
code['CCG'] = 'P'

code['CAU'] = 'H'
code['CAC'] = 'H'
code['CAA'] = 'Q'
code['CAG'] = 'Q'

code['CGU'] = 'R'
code['CGC'] = 'R'
code['CGA'] = 'R'
code['CGG'] = 'R'

code['AUU'] = 'I'
code['AUC'] = 'I'
code['AUA'] = 'I'
code['AUG'] = 'M'

code['ACU'] = 'T'
code['ACC'] = 'T'
code['ACA'] = 'T'
code['ACG'] = 'T'

code['AAU'] = 'N'
code['AAC'] = 'N'
code['AAA'] = 'K'
code['AAG'] = 'K'

code['AGU'] = 'S'
code['AGC'] = 'S'
code['AGA'] = 'R'
code['AGG'] = 'R'

code['GUU'] = 'V'
code['GUC'] = 'V'
code['GUA'] = 'V'
code['GUG'] = 'V'

code['GCU'] = 'A'
code['GCC'] = 'A'
code['GCA'] = 'A'
code['GCG'] = 'A'

code['GAU'] = 'D'
code['GAC'] = 'D'
code['GAA'] = 'E'
code['GAG'] = 'E'

code['GGU'] = 'G'
code['GGC'] = 'G'
code['GGA'] = 'G'
code['GGG'] = 'G'







hydrophobicityScale={}
#Using the POPC bilayer interface hydrophobicity scale found on http://blanco.biomol.uci.edu/hydrophobicity_scales.html  
hydrophobicityScale['W'] = -1.85
hydrophobicityScale['F'] = -1.13
hydrophobicityScale['Y'] = -0.94
hydrophobicityScale['L'] = -0.56
hydrophobicityScale['I'] = -0.31
hydrophobicityScale['C'] = -0.24
hydrophobicityScale['M'] = -0.23
hydrophobicityScale['G'] = -0.01
hydrophobicityScale['V'] = 0.07
hydrophobicityScale['S'] = 0.13
hydrophobicityScale['T'] = 0.14
hydrophobicityScale['A'] = 0.17
    #H = 0.17
hydrophobicityScale['N'] = 0.42
hydrophobicityScale['P'] = 0.45
hydrophobicityScale['Q'] = 0.58
hydrophobicityScale['H'] = 0.96
hydrophobicityScale['R'] = 0.81
hydrophobicityScale['K'] = 0.99
hydrophobicityScale['D'] = 1.23
hydrophobicityScale['E'] = 2.02
