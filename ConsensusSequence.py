"""Library containing the needed functions for assignemnt 1 for the Bioinformatics class.
Please read: project_1_transmembrane_regions_2016.pdf for more info.
"""

import TransmembraneFinder



def main():
    
    nucleicAcidSeqList = TransmembraneFinder.readInput("Sequences/Group2TMseqs.txt")    
    aminoAcidSeqList = [0 for x in range(0, len(nucleicAcidSeqList))]
    
    for index in range(0,len(nucleicAcidSeqList)):
        amino = TransmembraneFinder.translate(nucleicAcidSeqList[index])
        aminoAcidSeqList[index] = amino
                
                
    bestWinSize = TransmembraneFinder.bestWindowSize(aminoAcidSeqList, 6, 10)
    
    #print("Best window size ---------------------------------------------------",bestWinSize)
    
    aaSeqSegmentList = ["" for i in aminoAcidSeqList]
    
    for aminoIndex in range(0, len(aminoAcidSeqList)):
        #print("Amino seq ======================================================", aminoAcidSeqList[aminoIndex])
    
        hydrophobic = TransmembraneFinder.mostHydrophobicRegion(aminoAcidSeqList[aminoIndex], bestWinSize)        
        #print("Hydrophobic region is ==========================================", hydrophobic)
        
        hydrophobicSeq = hydrophobic['mostHydReg']
        aaSeqSegmentList[aminoIndex] = hydrophobicSeq    
    
    #print("Hydrophobic segment of best window size :\n",aaSeqSegmentList)
    
    dictionaryList = TransmembraneFinder.gatherContributions(aaSeqSegmentList)
    #print("\n\n\n\n List of Dictionary : \n",len(dictionaryList))
    
    

    test = 'ATGGCCGCCCCCTGCGAAGGCCTGCCCCCTCCACCACGGCCTGTCCGGCAGGATGCCCTCATCCGCCACAGTGCCCTTCCTGTCGCGGAGTCGAGCGTGAACCTCGAAAAAGTCGTCTTTGCCTTCTTCACGATCCTGGCCTGTACCCTGAACTTCGGGTTCTTTCTGGGCGAGATCGACCGTGCCGACTTCCACCACCCGGCCGAGCTGTTCATCGCCGTGGTCATCAACCTGATCACGCTGATCATCAAGTTCGGCGACCGTACCCAGATGGGCGCCACGCACCTGGCCACCAGCCTGGTGGCGACGCTGCAGCTGCTCTTTGCCTCGCTGGTCTGGATGTGGGTCGAACAGTTCAACAACACCCCGCTGGACGGCCACACGGTCAGCATCATCGTGTCGCTGTCGGGCGGTGCGCTGCTGGCCAACCTGGTCTCGGTCATCCTGCTGATCGGCGAGACGCTGCGCCAGACGCGCTGA'
    print("test amino Seq length = ",len(test))
    aatest = TransmembraneFinder.translate(test)
    print("aatest = ", aatest);
    newSeq = TransmembraneFinder.findHydrophobicRegions(dictionaryList, aatest)
    print("new hydro phobic sequence is ", newSeq)
    
    

main()



