# Bioinformatics-1
First Bioinformatics Assignment (Jan/2016)

Given a set of known transmembrane genes on the sequences folder
this program converts these DNA sequences into protein sequences 
and find the most likely candidate portion of each protein sequence 
(that is the same width across all sequences) that could be the transmembrane region of the protein. 

It builds a probabilistic consensus sequence for this region based upon your set of sequences and 
use this model to predict the transmembrane regions in an unknown protein that may have multiple transmembrane regions.

# Building it locally

 - Clone the repo: ``` $git clone https://github.com/bressan3/Bioinformatics-1.git ```.
 - Make sure you have [Matplotlib](http://matplotlib.org) graphing library installed.
 - Run it using your python3 interpreter. Ex.: ``` $python3 ConsensusSequence.py ```
