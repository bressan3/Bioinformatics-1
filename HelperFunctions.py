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


def entropy(probList):
    """Calculate and returns the entropy based on a given list of probabilities
    Args:
        probList (List): List of probabilities
    Returns:
        float: Return a number representing the entropy
    """
    # entropy = 0
    # for p in probList:
    # if p > 0:
    # entropy += p * math.log(p, 2)
    # return - entropy
    return -[sum(p * math.log(p, 2) for p in probList if p > 0)]


def information(entropy, numBases):
    """ Accepts two numbers representing the entropy
        and the number of bases and returns the information
    Args:
        entropy (float): Entropy's value
        numBases (int): Number of Bases
    Returns:
        float: Returns the information
    """
    return []


def calcProbs(freqCounts):
    """Takes a list of frequence counts of each base and returns a list
        containing the frequence for each one of them
    Args:
        freqCounts (int): A list of frequence counts
    Returns:
        Float: List of probabilities
    """
    return [freqCounts[x]/len(freqCounts) for x in range(0, 4)]
