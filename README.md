# 601.647-Project
Final project for 601.647 - Computational Genomics

Goal
The goal of our project is to test out the best alignment method we can apply on aligning known Protein Sequence With DNA Reads. Our project used several methods for both global and local alignment including K-mer, Boyer-Wheeler, Smith-Waterman, and FM index.

________________
Prerequisites
Before running our code, please make sure to download the SQLite file(around 5GB) for known protein bases.

https://livejohnshopkins-my.sharepoint.com/:u:/g/personal/ssun49_jh_edu/EWxc3ZQDW4ZBohelBqqLhhAB1SdzWza4Q-AINFcl_-JvFQ?e=aVC2f3

After downloading the SQLite file, please set the path in our python file on Line 9
 
________________
Run The File
Use command  python3 cg_project.py dna.txt  to run the file. Dna.txt should contain only the DNA sequence and nothing else.

________________
Test cases

ACTTACACTATACGCATGTAAAATTAGAAATAAAGTACGATTTTTAGACTGATGTATAATATAATT
ATGTAGATGTGATGAGTTTCTTTTATATGCTTCACCTGTCGGATCGGTCTTGACATAGCGTAG
TAATACATAAGCAAATAATAATTAACTACTGTGTATTTGTTATAACATCAGACAGTTTAAGTTGG
GATAATAGGAGCCACAATATACAATTTATCACATTTTAACTTAATTGACATATTACCATAAATGAC
TAACCATATCACTAGTTTTTAGATAACCTGATATAGTGATTATGAAAGGTTTATGGAATATAGCT
ATGACTTACTTAACTACATATGAAAGAAAAAACTTTTGTGTATTTATATGTTCACTCGTCTATTA
CTCATGCTTGAAGATTATATAAGTTGTGAGATGCAGGAAAAGTTCTTAACTTTCTCATAGGAC
GTTAACTTATTCTTTAATAGAGCATTTTATTCGAGCATGACAATAAGTACGCTGGCAACGATCC
GCTGAGCCAGTTCTTAATTATATTAATTTTTATTCTTATTAAAGTTTAGAGTTAATGATCAGGATT
ATTGCTTTGCAATAAATTTCTTATTCACAGTCATTCATATTGAGTTACTCGATAGATTATCAACT
TGATTCAGTCTGTTAGGCCATAATTACGTAAGTTAGGATTCTGAACTGCGTTGTATAATCGAAT
CTAATTCGACTTTTATAACTGCAAACTCTAATTTATTTAGATAACATGATTAGTTGAAGTTATATG
GGATATTTACCGTAAAATCCTCTTCGGGTGTCCTTCCTTTATTTGATGATAAACAGCCGCTAC
CATTATTAATTAATATAAGAAATAGTGATATTATCATAAATTTAGCACATTATTTTTATAGATATAGA
ATCATTTAATTACGTGCCGAAGTCTTATGACAAAATTGATAGACAATAATTCAAGTAATACTAAA
AGTCCATAGCACATACATTCTAATGTGATATACGTATAATTTAACCA

Protein: GNDPLSQ


AGCATTAAAAATCTAAACCTTCAGAAAGTGAAGATCCCGAGTATAGACCTTTATCTGCGGTTC
AAGTTAAGCATAAGGCTGCATGCTATCTTGTTACACCTACACTGCTCGAAGTAAATATAGGAA
GCGTGCGACCTGGCTCCAGGTGTTCCGCATCGTCACGTATTCGTTAACTGTTAATTGGTGA
CACATAAATAATATTATAGTCTCTCAAATTCAGCTCAGTTATCTTGAGCGTTATGTCTTAAATGG
CGTAGAACAGCATTGACTGTTTGACACTAACTGGTGTTCGGTTCGGTAATGGAGAATCTATG
CGGCAATGTCATTAATACATTTGAAACACGCCATATCGATACTGAACAAATCAATGCAAACTTC
CATGTTAGAATAAGGATAAACATACAAGTCGATAGAAAATGGGTAGGGGCTTTTAATTCATCC
AACACTCTACGACTTCTTCAAGAGCTAGTAGAGCACCCTGTAGTTGGAAAAGAAATTATTTC
GTAAGGTAAGCTCATACCGTCTTTCTTGCGGAAGACTTAACACGATAAGAAGTTGGAATAGTT
TCAAACGATAGTTATTAATCCTAATAACGGAACGTTATCTGAAGAATAAGTCTGACGGAGTGTA
ACTCGATATGGGCGGCAAATGGAGCAAAAGCCAGTTACTCACTATTCGAACTGGGCGAAAG
ATCCCAGCGCTTATGCACTTAATCCCGAGACTTGACCCGATATATGAGCTTAGACTAAAGCG
GGACTGTTGACGTTTGGAGTTGAAAAAATCTATTATACTAATCGGCTTCAACGTGCTCTACAG
CAGGCACTTGACGAGGGGCCCACACCGAGGAAGTAAACTGTTATACGTTGGGGATAGTGGT
AACTAATTAAGATGCTTGCCACAATAATAGTATCAAACCCGTATAAAGAGAACATCCACACTTT
AGTGAATTGAAGCGCGGCATCAGAATTTCCTTTTGGATACCTGATACAAAGCCCATCGTGAT
TTTTAGATTTCGTATATTTAC

Protein: MGGKWSKS
