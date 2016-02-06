"""Library containing the needed functions for assignemnt 1 for the Bioinformatics class.
Please read: project_1_transmembrane_regions_2016.pdf for more info.
Authors:
    Suman, Lucas, Stephane, Magadi
"""
from codes import *
import HelperFunctions
import matplotlib.pyplot as plt


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
        probabilityPerCount = HelperFunctions.calcProbs(oneCountList)
        probabilityList = [value for key, value in probabilityPerCount.items()]
        
        entropy = HelperFunctions.entropy(probabilityList)
        info = HelperFunctions.information(entropy, 4)
      
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
    
    cutOffValue = -0.25
    
    windowSize = len(listOfDicts)
    contributionValueList = [0 for x in range(0, len(aaSeq) - windowSize)]
    for index in range(0, len(aaSeq) - windowSize):
        #taking segment of amino acid
        partialAASeq = aaSeq[index:index + windowSize]
        partialAASeqSum = 0
        
        #calculating contribution value of the segment
        for i in range(0, len(partialAASeq)):
            contributionValue = listOfDicts[i].get(partialAASeq[i])
            partialAASeqSum = partialAASeqSum + contributionValue
        
        #storing start index of segment and its total contribution
        contributionValueList[index] = {'index': index, 'sum': partialAASeqSum}
        #print("Index = ",index, "value = ", partialAASeqSum)
    
    #holds index of segments which have higher contribution value than cut off value
    hydrophobicRegionIndex = []
    

    for x in range(0, len(contributionValueList)):
        
        if contributionValueList[x]['sum'] > cutOffValue:
            hydrophobicRegionIndex.append(contributionValueList[x]['index'])

    if (hydrophobicRegionIndex == []):
        print("NO hydrophobic region")
        return "NO hydrophobic Region"

    #contains start and end index of hydrophobic part alternatingly. Example: for 3 region [start, end, start, end, start, end]
    refinedHRIndex = [] 
    
    refinedHRIndex.append(hydrophobicRegionIndex[0])
    
    for index in range(1, len(hydrophobicRegionIndex)):
        #if start index of following sequence lies within windowSize range then combine the sequence with previous one
        #if start index of following sequence is farther than windowSize range then that following index is start index of new hydrophobic sequence
        if ((hydrophobicRegionIndex[index -1] + windowSize) > hydrophobicRegionIndex[index]):
            continue
        else:
            refinedHRIndex.append(hydrophobicRegionIndex[index -1] + windowSize)
            refinedHRIndex.append(hydrophobicRegionIndex[index])
    
    #holds list of hydrophobic segment region 
    hydropgobicRegionList = []
    index = 1
    while(index < len(refinedHRIndex)):
        x = refinedHRIndex[index - 1]   #start index
        y = refinedHRIndex[index]       #end index
        hydropgobicRegionList.append(aaSeq[x : y])
        index = index +2
        
    #if there is extra start region then gather all amino acid from that region to length of window size    
    if (len(refinedHRIndex) % 2 != 0):
        lastx = refinedHRIndex[-1]
        lasty = lastx + windowSize
        hydropgobicRegionList.append(aaSeq[lastx : lasty] )        
    
    return hydropgobicRegionList
    


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
    
    
    cutOffValue = -0.25
    windowSize = len(listOfDicts)
    
    contributionValueList = [0 for x in range(0, len(aaSeq) - windowSize)]

    for index in range(0, len(aaSeq) - windowSize):
        #taking segment of amino acid
        partialAASeq = aaSeq[index:index + windowSize]
        partialAASeqSum = 0
        
        #calculating contribution value of the segment
        for i in range(0, len(partialAASeq)):
            contributionValue = listOfDicts[i].get(partialAASeq[i])
            partialAASeqSum = partialAASeqSum + contributionValue
        
        #storing start index of segment and its total contribution
        contributionValueList[index] = {'index': index, 'sum': partialAASeqSum}
        
    contributionInformation = [0 for x in range(0, len(contributionValueList))]
    
    
    for x in range(0, len(contributionValueList)):
        # Neglecting contribution of segment whose contribution value is less than the cutoff value
        # Otherwise subtracting it from negative of cutoff value to make it easier to plot along positive y axis
        if contributionValueList[x]['sum'] < cutOffValue : 
            contributionInformation[x] = 0
        else:
            contributionInformation[x] = -cutOffValue - contributionValueList[x]['sum'] 
        
       
    plt.xlabel('Position in the amino acid sequence')
    plt.ylabel('Contribution')
    xaxis = [x for x in range(0, len(aaSeq) - windowSize)]
    yaxis = contributionInformation

    
    plt.plot( xaxis, yaxis)
    plt.show()
