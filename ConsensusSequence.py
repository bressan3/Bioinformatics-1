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
                
    #calculating best window size           
    bestWinSize = TransmembraneFinder.bestWindowSize(aminoAcidSeqList, 6, 10)
    
    aaSeqSegmentList = ["" for i in aminoAcidSeqList]
    
    for aminoIndex in range(0, len(aminoAcidSeqList)):
        #finding most hydrophobic region
        hydrophobic = TransmembraneFinder.mostHydrophobicRegion(aminoAcidSeqList[aminoIndex], bestWinSize)               
        hydrophobicSeq = hydrophobic['mostHydReg']
        aaSeqSegmentList[aminoIndex] = hydrophobicSeq    
    
    #creating a dictionary list
    dictionaryList = TransmembraneFinder.gatherContributions(aaSeqSegmentList)
    
    #loading sequences from FinalTestSeqs.txt
    testSeqList = TransmembraneFinder.readInput("Sequences/FinalTestSeqs.txt")    
       
    for index in range(0,len(testSeqList)):
        amino = TransmembraneFinder.translate(testSeqList[index])
        print("Test data ",index," Amino acid :\n", amino)
        newSeq = TransmembraneFinder.findHydrophobicRegions(dictionaryList, amino)
        print("Hydro phobic sequence for test data",index, "is :", newSeq)
        
        #Constructing graph
        TransmembraneFinder.constructGraph(dictionaryList, amino)
    
main()



