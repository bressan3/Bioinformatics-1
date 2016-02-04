"""Library containing the needed functions for assignemnt 1 for the Bioinformatics class.
Please read: project_1_transmembrane_regions_2016.pdf for more info.
"""

import TransmembraneFinder
import HelperFunctions


def main():
    
    nucleicAcidSeqList = TransmembraneFinder.readInput("Sequences/Group2TMseqs.txt")    
    aminoAcidSeqList = [0 for x in range(0, len(nucleicAcidSeqList))]
    
    for index in range(0,len(nucleicAcidSeqList)):
        amino = TransmembraneFinder.translate(nucleicAcidSeqList[index])
        aminoAcidSeqList[index] = amino
                
                
    bestWinSize = TransmembraneFinder.bestWindowSize(aminoAcidSeqList, 6, 10)
    
    print("Best window size ---------------------------------------------------",bestWinSize)
    
    aaSeqSegmentList = ["" for i in aminoAcidSeqList]
    
    for aminoIndex in range(0, len(aminoAcidSeqList)):
        print("Amino seq ======================================================", aminoAcidSeqList[aminoIndex])
        hydrophobic = TransmembraneFinder.mostHydrophobicRegion(aminoAcidSeqList[aminoIndex], bestWinSize)        
        hydrophobicSeq = hydrophobic['mostHydReg']
        aaSeqSegmentList[aminoIndex] = hydrophobicSeq    
    
    print("Hydrophobic segment of best window size :\n",aaSeqSegmentList)
    
    #dictionaryList = TransmembraneFinder.gatherContributions(aminoAcidSeqList)
    

    #test = 'ATGTTTGGGCCCCCAGCCGAGCAATTGCCTCCCGTATTGCTGGTCGAAATGGGATTCGTGAGTTCCATGCCGGTTCTGGGGCTCGTTATACTCATTTTCACCCAAATCCCCGTAACGGATCTCATGTTTGGTATGGGTCGCCACGTTGATCTGTTTTGA'
    #newSeq = TransmembraneFinder.findHydrophobicRegion(dictionaryList, test)
    #print(newSeq)
    
    

def test():
    #count = HelperFunctions.gatherCounts(['FLS','LSS','LSY','LSF'])
    
    #count = HelperFunctions.gatherCountsDNATest(['ACG','CGG','CGT','CGA'])
    #print(count)
    
    
    #entropy = HelperFunctions.entropy([0.33,0.14,0.41,0.12])
    #print(entropy)
    
    #info = HelperFunctions.information(entropy,4)
    #print(info)
    
    #contribution = TransmembraneFinder.gatherContributions(['ACGGTCGT', 'CCGTGCGC', 'ACGGACGT', 'GGCCGCAA'])
    
    #contribution = TransmembraneFinder.gatherContributions(['FLS','LSS','LSY','LSF'])
    
    #print("Contribution = ",contribution)
    print("nothing")

main()



