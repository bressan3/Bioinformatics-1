"""Library containing the required helper functions for assignemnt 1 for the Bioinformatics class.
Please read: project_1_transmembrane_regions_2016.pdf for more info.
"""
import math


def gatherCounts(seqList):
    """Accepts a list of sequences and returns a list of the frequence
        counts of each amino acids at each position
    Args:
        seqList (char): List of aminoacid sequences
    Returns:
        Int: List of the frequency counts of each aminoacids at each position
    """
    returnList = [0 for x in seqList[0]]
    # print("Length of return list = ",len(returnList))
    for index in range(0, len(seqList[0])):
        aminoSeq = {'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0, 'E': 0, 'Q': 0, 'G': 0, 'H': 0,
                    'I': 0, 'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0}

        for i in range(0, len(seqList)):
            aminoSeq[seqList[i][index]] = aminoSeq[seqList[i][index]] + 1

            returnList[index] = aminoSeq
        # print(returnList[index])

    return returnList

'''
def gatherCountsDNATest(seqList):
    """Accepts a list of sequences and returns a list of the frequence
        counts of each amino acids at each position
    Args:
        seqList (char): List of aminoacid sequences
    Returns:
        Int: List of the frequency counts of each aminoacids at each position
    """
    returnList = [0 for x in seqList[0]]
    print("Length of return list = ", len(returnList))
    for index in range(0, len(seqList[0])):
        a = 0
        c = 0
        g = 0
        t = 0
        for i in range(0, len(seqList)):
            if (seqList[i][index] == 'A'):
                a = a + 1
            if (seqList[i][index] == 'C'):
                c = c + 1
            if (seqList[i][index] == 'G'):
                g = g + 1
            if (seqList[i][index] == 'T'):
                t = t + 1
        print("A = ", a, " C = ", c, " G = ", g, "T = ", t)
        returnList[index] = [a, c, g, t]
    return returnList
'''


def entropy(probList):
    """Calculate and returns the entropy based on a given list of probabilities
    Args:
        probList (List): List of probabilities
    Returns:
        float: Return a number representing the entropy
    """
    return -(sum(p * math.log(p, 2) for p in probList if p > 0))


def information(entropy, numBases):
    """ Accepts two numbers representing the entropy
        and the number of bases and returns the information
    Args:
        entropy (float): Entropy's value
        numBases (int): Number of Bases
    Returns:
        float: Returns the information
    """
    return (math.log(numBases, 2) - entropy)


def calcProbs1(freqCounts):
    """Takes a list of frequence counts of each base and returns a list
        containing the frequence for each one of them
    Args:
        freqCounts (int): A list of frequence counts
    Returns:
        Float: List of probabilities
    """
    return [freqCounts[x] / len(freqCounts) for x in range(0, len(freqCounts))]


def calcProbs(freqCounts):
    """Takes a list of frequence counts of each base and returns a list
        containing the frequence for each one of them
    Args:
        freqCounts (int): A list of frequence counts
    Returns:
        Float: List of probabilities
    """
    total = 0
    for key, value in freqCounts.items():
        total = total + value

    for key, value in freqCounts.items():
        freqCounts[key] = value / total

    return freqCounts
