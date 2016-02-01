"""Library containing the needed functions for assignemnt 1 for the Bioinformatics class.
Please read: project_1_transmembrane_regions_2016.pdf for more info.
"""
from codes import *

"""function to reads to contents of the file and output a list containing the nucleic acid sequences """
def readInput(filename):
    file = [x.strip() for x in open(filename).readlines()]
    return file
		

def translate(dna):
    """Function that takes a DNA string converts it to RNA and then translates and returns its aminoacid sequence
    Args:
        dna (string): DNA sequence
    Return:
        string: Aminoacid sequence held by the DNA string given as a parameter to the function
    """
    rna = dna.replace('T', 'U')
    startIndex = dna.find('AUG') + 1
    aminoAcidsSeq = ""
    for i in range(startIndex, len(rna), 3):
        aminoAcidsSeq += code[rna[i: i+3]]
        if aminoAcidsSeq[len(aminoAcidsSeq) - 1] == '*':
            aminoAcidsSeq = aminoAcidsSeq[:-1]
            break
    return aminoAcidsSeq



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
    
    
def mostHydrophobicRegion(aaSeq,winSize):
    
    global hydrophobicityScale    
    index=0
    aaSeqHydReg=0
    mostHydReg=0
    while index<len(aaSeq):
        if (len(aaSeq)-index)!=winSize-1:
            aaSeqwin=aaSeq[index:index+winSize]
        else:
            break
        count=0
        aaSeqReg=0
        while count<winSize:
            aaSeqReg=aaSeqReg+hydrophobicityScale[aaSeqwin[count]]
            count=count+1 
        if aaSeqHydReg<aaSeqReg:
            aaSeqHydReg=aaSeqReg 
            mostHydReg=index
        index=index+1
    print(mostHydReg)
    
mostHydrophobicRegion('KDFEEE',3)
