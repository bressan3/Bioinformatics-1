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

#count = HelperFunctions.gatherCounts(['FLS','LSS','LSY','LSF'])

#count = HelperFunctions.gatherCountsDNATest(['ACG','CGG','CGT','CGA'])
#print(count)
'''
list1 = {'x':1,'y':3,'z':0,'a':0}
prob = HelperFunctions.calcProbs(list1)
print(prob)

li = [value for key, value in prob.items()]
print(li) 

en = HelperFunctions.entropy(li)
print(en)

inf = HelperFunctions.information(en, 4)
print(inf)

#final = [value*inf for key,value in list1.items()]
#print(final)

for key,value in list1.items():
    if (list1[key] == value): 
        list1[key] = value * inf
        

print(list1)
'''
#entropy = HelperFunctions.entropy([0.33,0.14,0.41,0.12])
#print(entropy)

#info = HelperFunctions.information(entropy,4)
#print(info)

#contribution = TransmembraneFinder.gatherContributions(['ACGGTCGT', 'CCGTGCGC', 'ACGGACGT', 'GGCCGCAA'])

#contribution = TransmembraneFinder.gatherContributions(['FLS','LSS','LSY','LSF'])

#print("Contribution = ",contribution)


#main()


dict1 = TransmembraneFinder.gatherContributions(['FLS','LSS','LSY','LSF'])
print(dict1)
