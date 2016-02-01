"""Library containing the needed functions for assignemnt 1 for the Bioinformatics class.
Please read: project_1_transmembrane_regions_2016.pdf for more info.
"""

import TransmembraneFinder
import HelperFunctions


def main():
     nucleicAcidSeqList = TransmembraneFinder.readInput("Group2TMseqs.txt")    
    # print(nucleicAcidSeqList)
     print(TransmembraneFinder.translate('ATGTTGATTTAATAGTTT'))
    
  #  mostHydrophobicRegion('MNQSTKRKHVLIPLPPVSTKRHD',8)
  #  bestWindowSize(['MNQSTKRKHVLIPLPPVSTKRHD','MNQSTKRKHVLIPLPPVSTKRHD','MNQSTKRKHVLIPLPPVSTKRHD'],6,10)


main()


