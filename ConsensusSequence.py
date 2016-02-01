"""Library containing the needed functions for assignemnt 1 for the Bioinformatics class.
Please read: project_1_transmembrane_regions_2016.pdf for more info.
"""

import TransmembraneFinder
import HelperFunctions


def main():
    nucleicAcidSeqList = TransmembraneFinder.readInput("Sequences/Group2TMseqs.txt")    
    #print(nucleicAcidSeqList)
 
    for x in nucleicAcidSeqList:
        amino = TransmembraneFinder.translate(x)
        print(amino)       
        print("Amino seq Len =",len(amino))
        print("DNA seq len =",len(x))



    
    # mostHydrophobicRegion('MNQSTKRKHVLIPLPPVSTKRHD',8)
    # bestWindowSize(['MNQSTKRKHVLIPLPPVSTKRHD','MNQSTKRKHVLIPLPPVSTKRHD','MNQSTKRKHVLIPLPPVSTKRHD'],6,10)

count = HelperFunctions.gatherCounts(['FLS','LSS','LSY','LSF'])

#count = HelperFunctions.gatherCountsDNATest(['ACG','CGG','CGT','CGA'])
print(count)


prob = HelperFunctions.calcProbs(count[0])
print(prob)

entropy = HelperFunctions.entropy([0.33,0.14,0.41,0.12])
print(entropy)

info = HelperFunctions.information(entropy,4)
print(info)
#contribution = TransmembraneFinder.gatherContributions(['ACGGTCGT', 'CCGTGCGC', 'ACGGACGT', 'GGCCGCAA'])




#main()