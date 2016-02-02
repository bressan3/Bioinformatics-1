"""Library containing the needed functions for assignemnt 1 for the Bioinformatics class.
Please read: project_1_transmembrane_regions_2016.pdf for more info.
"""
from codes import *
import HelperFunctions


def readInput(filename):
    """Function that reads a file's contents outputs a list containing a nucleic acid sequence
    Args:
        filename (string): File path
    Return:
        string: Sequence read from the file
    """
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


def mostHydrophobicRegion(aaSeq, winSize):
    """Function that takes a amino acid sequence and a window size and returns the start position of
    a window which is the most hydrophobic region of the given sequence
    Args:
        aaSeq (string): Amino Acid sequence
        winSize (int): Window size
    Returns:
        int: Start location of a window of the given size (Most hidrophobic region of that size)
    """
    global hydrophobicityScale
    index = 0
    aaSeqHydReg = 0
    mostHydReg = 0
    while index < len(aaSeq):
        if (len(aaSeq) - index) != winSize - 1:
            aaSeqwin = aaSeq[index:index + winSize]
        else:
            break
        count = 0
        aaSeqReg = 0
        while count < winSize:
            aaSeqReg = aaSeqReg+hydrophobicityScale[aaSeqwin[count]]
            count = count + 1
        if aaSeqHydReg < aaSeqReg:
            aaSeqHydReg = aaSeqReg
            mostHydReg = index
        index = index + 1
    return {'mostHydReg': mostHydReg, 'score': aaSeqHydReg}


def bestWindowSize(listOfAASeqs, lowWinSize, highWinSize):
    """Takes in a list of amino acid sequences and two numbers representing low and high values for
     possible window sizes. Returns the window size which has the highest total hydrophibcity across
      all of the amino acid sequences
    Args:
        listOfAASeqs (string): List of amino acid sequences
        lowWinSize (int): low value for possible window size
        highWinSize (int): low value for possible window size
    Returns:
        int: window size
    """
    index = 0
    wholesum = [0 for x in range(highWinSize - lowWinSize + 1)]
    while index <= (highWinSize-lowWinSize):
        tempsum = 0
        for i in listOfAASeqs:
            hydrophobicityInfo = mostHydrophobicRegion(i, lowWinSize + index)
            tempsum = tempsum + hydrophobicityInfo['score']
        wholesum[index] = tempsum
        index = index + 1
    bestWinSize = lowWinSize + wholesum.index(max(wholesum))
    print(bestWinSize)
    return bestWinSize


def gatherContributions(aminoSeqList):
    countsList = HelperFunctions.gatherCounts(aminoSeqList)
        
    finalDict = [0 for x in range(len(countsList))]
    index =0
    for oneCountList in countsList:
      #  print("onecount at a time : \n",oneCountList)
        
        probabilityPerCount = HelperFunctions.calcProbs(oneCountList)  

        probabilityList = [value for key, value in probabilityPerCount.items()]
        print("Probability list : \n", probabilityList)

        entropy = HelperFunctions.entropy(probabilityList)
        print("Entropy = ", entropy)
        

        info = HelperFunctions.information(entropy, 4)
        print("Information = ",info)

        for key,value in probabilityPerCount.items():
            if (probabilityPerCount[key] == value): 
                probabilityPerCount[key] = value * info
        
        finalDict[index] = probabilityPerCount
        index = index+1
        
    return finalDict
