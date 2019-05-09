# stLFR reads simulator

## Introduction

A easy-to-use and easy-to-analysis simulator for stLFR reads.

## User's Guide

### Installation

```
git clone https://github.com/BGI-Qingdao/stLFR_reads_sim.git YOUR-INSTALL-DIR
cd YOUR-INSTALL-DIR/main && make && cd -
cd YOUR-INSTALL-DIR/bin && stLFRSim --help && cd -
```

### Usage

#### Distribution file 

Distribution file contain the expected-discrete-distribution information.
It is a three columns format and each column refers to :

* areaBegin
* areaLength
* areaWeight

The random engine will first choose an areaBegin by areaWeight, then it random choose a final result from areaBegin to areaBegin+areaLength -1 .


Here is an example of distribution file

```
300 50  50
350 100 30
450 200 20
```
And the corresponding distribution looks like :

```
      _
     |*|
     |*|   _
     |*|  |*|    _
     |*|  |*|   |*|
     |*|  |*|   |*|
 ------------------->
  300~~349~~~449~~~649
```

#### Baisc usage :

```
./bin/stLFRSim args 
                --ref                           [required]      [ string arg ]  reference fasta file
                --o_prefix                      [required]      [ string arg ]  output file prefix . print into o_prefix.1.fa.fq && o_prefix.2.fa.fq
                --lr_length_distribution        [required]      [ string arg ]  distribution file of long read length
                --pe_num_distribution           [required]      [ string arg ]  distribution file of number of read-pair in 1 long read
                --if_lenth_distribution         [required]      [ string arg ]  distribution file of insert fragment length
                --readpair_num                  [required]      [ long arg ]    total number of final generated read-pairs
                --mutation_rate                 [optional]      [ float arg ]   mutation rate    [ default= 0.005 ]
                --insert_percent                [optional]      [ float arg ]   insert percent   [ default= 0.005 ]
                --delete_percent                [optional]      [ float arg ]   delete percent   [ default= 0.005 ]
                --substitute_percent            [optional]      [ float arg ]   substitute percent   [ default= 0.99 ]
                --max_slr_cov                   [optional]      [ float arg ]   max single long read cov      [ default= 0.5 ]
                --read_len                      [optional]      [ int arg ]     read length      [ default= 100 ]
```

- running example

```
stLFRSim_Main --ref chr19.fa --o_prefix testsim --lr_length_distribution lr_length_dis.txt --pe_num_distribution pe_num_dis.txt --if_lenth_distribution pe_length_dis.txt --readpair_num 15000000 
```

#### Output format

- format details

To simplify the explanation , I need to define some symbol first :

*Ri* : a pair of "read1 and read2" , index by i .

*IFi* : the insert fragment that generate *Ri* .

*LFi* : the long sequence framgment that generate *IFi* .

*RSi* : the reference sequence that generate *IFi* .


The first line of the fastq format contain below information :


Column 1 : @read_name#barcode_name/[1 or 2 ] .

Column 2 : barcode_num .

Column 3 : the sequence name of *RSi* .

Column 4 : the sequence length of *RSi* .

Column 5 : the start postion of *LFi* in *RSi* , index start from 0 .

Column 6 : the length of *LFi* in *RSi* 

Column 7 : the start postion of *IFi* in *LFi* , index start from 0 .

Column 8 : the length of *IFi* in *LFi*  .

Column 9 : the CIGAR string about how read mutation .


- read1 example 

```
@stlfrsim_1#barcode_1/1 1       chr19   59128983        55629801        48983   11651   592     8=1X39=1X51=
AAAAAAAATAAATTAGCTGGGTGTGTTGGCACGGGCCTGTAATCCCAGATACTTGGGAGGCTGAGGCAGGAGAATCGCTTGAACCAGGGAGTGGGAGCTT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@stlfrsim_2#barcode_1/1 1       chr19   59128983        55629801        48983   27896   385     100=
GGCGGAGGGAGGGAAGGGTGGTCTTGGAGGTTGGGGCCCGAGGATATCGGGGGTCCCCCCGGGCCCCCGACATCGGTCTCGGGAAGCGAAGCAGCCGCGG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

```

- read2 example

```
@stlfrsim_1#barcode_1/2 1       chr19   59128983        55629801        48983   11651   592     100=
CTATGCTTTGAATGTATGTATCCCCTCAAAATTTACATGTTGAAACCTGATCATGAATGCAATAGTATTCCACGGGACTTTAGGAGGTGACTAGGTCATG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@stlfrsim_2#barcode_1/2 1       chr19   59128983        55629801        48983   27896   385     100=
TCTCCAGAGGAGGCTGCGGAGGAGGAGGAGGAAGGTGAGGTCCCGGACTCCGCAGGTCTGGAGCTGGGGGGTGGGGGGGCGGGGACGCTGGGCCCGGGAG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@stlfrsim_3#barcode_1/2 1       chr19   59128983        55629801        48983   35530   525     82=1X9=1X7=
CACATGTTGGTAAATTAGCTCAGGCACTGGCCAGGGAATTGTGATTTGCATGTAGCTGGACCAGGTTATGCCAGTGGTTTTGCGAGGTGAGGTTGGAGCA

```

