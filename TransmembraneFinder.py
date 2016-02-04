"""Library containing the needed functions for assignemnt 1 for the Bioinformatics class.
Please read: project_1_transmembrane_regions_2016.pdf for more info.
"""
from codes import *
import HelperFunctions
import matplotlib as plt


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
        aminoAcidsSeq += code[rna[i: i + 3]]
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
    mostHydReg = ''
    aaSeqwin=[]
    while index < len(aaSeq):
        if (len(aaSeq) - index) != winSize - 1:
            aaSeqwin = aaSeq[index:index + winSize]
        else:
            break
        count = 0
        aaSeqReg = 0
        while count < winSize:
            aaSeqReg = aaSeqReg + hydrophobicityScale[aaSeqwin[count]]
            count = count + 1
        if aaSeqHydReg < aaSeqReg:
            aaSeqHydReg = aaSeqReg
            mostHydReg = aaSeqwin
        index = index + 1

    return {'mostHydReg': mostHydReg, 'score': aaSeqHydReg}


'''
def mostHydrophobicRegion1(aaSeq, winSize):
    global HydrophobicityScale
    aaSeqHydReg = 0
    mostHydRegIndex = 0
    
    for index in range(0,(len(aaSeq) - winSize)):
        
        aminoSegment = aaSeq[index:index + winSize]
        aaSegmentHReg = 0
        for i in range(0, len(aminoSegment)):
            aaSegmentHReg = aaSegmentHReg+hydrophobicityScale[aminoSegment[i]]
        if aaSeqHydReg < aaSegmentHReg:
            aaSeqHydReg = aaSegmentHReg
            mostHydRegIndex = index
            
    return {'mostHydReg': mostHydRegIndex, 'score': aaSeqHydReg}
'''


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
    while index <= (highWinSize - lowWinSize):
        tempsum = 0
        for aaSeq in listOfAASeqs:
            hydrophobicityInfo = mostHydrophobicRegion(aaSeq, lowWinSize + index)
            tempsum = tempsum + hydrophobicityInfo['score']
        wholesum[index] = tempsum
        index = index + 1
    bestWinSize = lowWinSize + wholesum.index(max(wholesum))
    return bestWinSize


def gatherContributions(aminoSeqList):
    """ Takes a list of amino acids and returns a list of dictionaries where every entry in the dictionary
     is the base and the contribution of that base to the information for that position
    Args:
        aminoSeqList (string): List of amino acids
    Returns:
        Dictionary: List of Dictionaries
    """
    countsList = HelperFunctions.gatherCounts(aminoSeqList)

    listOfDicts = [0 for x in range(len(countsList))]
    index = 0
    for oneCountList in countsList:
      #  print("onecount at a time : \n",oneCountList)

        probabilityPerCount = HelperFunctions.calcProbs(oneCountList)

        probabilityList = [value for key, value in probabilityPerCount.items()]
        print("Probability list : \n", probabilityList)

        entropy = HelperFunctions.entropy(probabilityList)
        print("Entropy = ", entropy)

        info = HelperFunctions.information(entropy, 4)
        print("Information = ", info)

        for key, value in probabilityPerCount.items():
            if (probabilityPerCount[key] == value):
                probabilityPerCount[key] = value * info

        listOfDicts[index] = probabilityPerCount
        index = index + 1

    return listOfDicts


def findHydrophobicRegions(listOfDicts, aaSeq):
    """ Takes in a list of dictionaries containing the contributions of
     the base to the information (output from the gatherContributions() function) and an amino acid
     sequence. Returns a list of strings where every string is a region of the amino
     acid that had a high match to the model of contributions
    Args:
        listOfDicts (Dictionary): List of dictionaries
        aaSeq (string): Amino acid sequence
    Returns:
        String: List of regions
    """
    windowSize = len(listOfDicts)
    contributionValueList = [0 for x in range(0, len(aaSeq) - windowSize)]
    for index in range(0, len(aaSeq) - windowSize):
        partialAASeq = aaSeq[index:index + windowSize]
        partialAASeqSum = 0
        for i in range(0, len(partialAASeq)):
            contributionValue = listOfDicts[0].get(partialAASeq[i])
            partialAASeqSum = partialAASeqSum + contributionValue
        # keeping track of aminoSeq used and their information value
        contributionValueList[index] = {'index': index, 'sum': partialAASeqSum}
        index = index + 1
    
    print("\n\n\n\n Contribution value list :\n", contributionValueList)

    hydrophobicRegionIndex = []

    i = 0
    print("\n\n\n\n\n")
    for x in range(0, len(contributionValueList)):
        print("checking contribution", contributionValueList[x])
        if contributionValueList[x]['sum'] > 2:
            hydrophobicRegionIndex[i] = contributionValueList[x]['index']
            i = i + 1

    return hydrophobicRegionIndex

    '''
    needs a checking point. if greater than that add amino seq to reutrning list
    '''
#    return contributionValueList


def constructGraph(listOfDicts, aaSeq):
    """ The function should create a graph where the X axis is the position in the
    amino acid sequence and the Y axis is how well the model matches at that location, graphically
    representing the same information as on findHydrophobicRegions().
    Args:
        listOfDicts (Dictionary): List of dictionaries containing the contributions of
        the base to the information (output from the gatherContributions() function)
        aaSeq (string): An amino acid sequence
    Returns:
        None: No Return
    """
    hydrophobicRegions = findHydrophobicRegions(listOfDicts, aaSeq)
    freqCounts = gatherCounts(hydrophobicRegions)
    probs = []
    for i in range(i, len(freqCounts)):
        probs[i] = calcProbs(freqCounts[i])
    informations = []
    for i in range(0, len(probs)):
        informations[i] = information(entropy(probs[i]), len(freqCounts[0]))
    plt.xlabel('Position in the amino acid sequence')
    plt.ylabel('Matches')
    plt.plot(informations, range(0, len(hydrophobicRegions) - 1))
    plt.show()
