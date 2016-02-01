# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 09:54:51 2016

@author: Suman
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