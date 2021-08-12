Part 1 – Quality Score Distribution per-nucleotide
1)  
```
   LABEL    |   Filename
--------------------------  
   Read 1   |   1294_S1_L008_R1_001.fastq.gz
   Read 2   |   1294_S1_L008_R4_001.fastq.gz
   Index 1  |   1294_S1_L008_R2_001.fastq.gz
   Index 2  |   1294_S1_L008_R3_001.fastq.gz
```

2) 
```
i) Look in respository
ii) For reads, because we are aligning a read to a reference genome, we don't need a strict cut off as high as 99.9-%, but a relatively good quality. For indexes, we want a individualized quality score instead of an average quality scores because a bad quality score for one base would mean the whole index is bad. Index reprsents what the sample read came from. So if there was one base in the index that had a bad quality score that could possibly be masked by the rest of quality scores. Thus the minimum cut off score for index quality is -log(0.001) = 30

iii) _ index have undetermined (N) base calls

Index 1
Command:
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | grep "N" | wc -l
3976613

Index 2
zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | grep "N" | wc -l
3328051
```

Part 2 – Develop an algorithm to de-multiplex the samples
```
1) Fix index swapping algorithmically. Index swapping happens in RNA-seq lab, my algorithm will be used to see if an index was properly ligated and see where the reads came from.
2) 2 Fastq files per bin, one contain r1 the other containing r2. A total of 24x2(2(swap/bad) x 2(r1/r2)) = 52 fastq output files
3) zcat <file> | sed -n '2~4p' to create unit files
4) Pseudocode located in responsitory in the form of a pdf doc
```

5) Read Below
```
A) Descriptions and Doc Strings tells readers what the codes does through out the script
B) Function headers
B) Test Examples

Hamming Function
1) Takes in two sequences and returns the number of incorrect matches.
2) Test. Input seq = AAA, seq2 = CCC

Reverse Complement
1) Takes a sequence and returns the reverse complement
2) Input seq = AAA, seq2 = TTT
```






