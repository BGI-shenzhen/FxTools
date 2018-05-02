# FxTools
<b>FxTools: a comprehensive toolkit for FASTA and FASTQ file manipulation</b>
 </br></br> FxTools is a full-featured toolkit for comprehensive analysis of both FASTA and FASTQ file, covering users' needs from sequence modification to data analysis. This tool consists of three parts: **Fatools** , **Fqtools** and **Formtools**, categorized by the types of files to deal with. In Formtools, we provide simple conversion and file manipulation for other NGS data files including BAM,SAM and SOAP.
FxTools is implemented in C/C++ language, available for Linux and Mac OS X operating system. 
###  1) Download and Install
------------
  <b> Download : git clone https://github.com/BGI-shenzhen/FxTools.git </b>
  </br> </br> Pre-installations of 4 libraries are required before installing FxTools
  </br> 1 htslib: [samtools-1.6/htslib-1.6](https://sourceforge.net/projects/samtools/files/samtools)
  </br> 2 boost : [boost](http://www.boost.org/) with [g++](https://gcc.gnu.org/) > 4.9 is recommended
  </br> 3 zlib  : [zlib](https://zlib.net/) > 1.2.8 is recommended
  </br> 4 ncurses: [ncurses](https://www.gnu.org/software/ncurses/)  >5.7 is recommended

 For <b>linux /Unix </b>  statics 
- you can use the statically compiled programs directly 
     
 <pre>
           git clone https://github.com/BGI-shenzhen/FxTools.git
           cd FxTools-XXX;  
	   chmod 775 bin/FxTools_Linux ; 
	   ./bin/FxTools_Linux

</pre>

- To compile FxTools, do [./configure] first and than [make]
- Final software can be found in the direcoty [bin/FxTools]

 For <b>linux /Unix </b>
<pre>
        git clone https://github.com/BGI-shenzhen/FxTools.git
        cd FxTools-XXX;                             
        ./configure
        make ;                    
        ../bin/FxTools                             
</pre>


### 2) Features 
------------

### Parameter description</b>
```php

Program: FxTools
Version: 0.16   hewm2008@gmail.com/xuxiaomin@bgi.com     2018-4-15

        Usage:

                Fatools        Tools For Fasta
                Fqtools        Tools For Fastq
                Formtools      Tools For Form convert

                Help           Show help in detail

```

#### Fatools

|Module |    Function   |       Description                                                |
|:-----:|:--------------|:-----------------------------------------------------------------|
|Summary|               |                                                                  |
|       |stat           |statistics of FASTA                                               |
|       |dict           |generate a header file for FASTA                                  |
| Split |               |                                                                  |
|       |split          |split FASTA. Default is splitting by se-quence ID                 |
|       |rand           |randomly sample FASTA by proportion                               |
| Search|               |                                                                  |
|       |findN          |find the regions of N in FASTA file                               |
|       |findSubSeq     |find the region containing the subse-quences                      |
|       |grep           |search for the target subsequence                                 |
|       |extractP       |extract sequences with specific ID                                |
|       |extractN       |extract sequences by specified order range                        |
|       |getCdsPep      |find CDS & peptide sequences (GFF re-quired)                      |
|       |sort           |sort the FASTA by sequence ID or length                           |
|Modify |                                                                                  |
|       |filter         |remove the sequences either too short or with too many missing (N)|
|       |reform         |edit the FASTA (reverse, comple-ment,etc.)                        |
|       |mergaSca       |reform current FASTA into new scaffolds                           |
|       |BaseModify     |modify a single base in FASTA                                     |
|       |ChangPosi      |change FASTA back to scaffolds                                    |
|       |Fa2Fq          |covert FASTA to FASTQ                                             |

#### Fqtools

|Module |    Function   |       Description                                                |
|:-----:|:--------------|:-----------------------------------------------------------------|
|Summary|               |                                                                  |
|       |valid          |check fq and change to valid Fq                                   |
|       |stat           |statistics of FASTQ                                               |
|       |fqcheck        |base and quality distribution                                     |
| Split |               |                                                                  |
|       |splitpool      |split pooling FASTQ to samples for RAD(GBS)                       |
|       |splitFq        |split FASTQ by sepcifying number of sequences in each file        |
|       |cut            |extract subsequence in FASTQ                                      |
|       |rand           |randomly sample FASTQ by proportion                               |
| Modify|               |                                                                  |
|       |filter         |filter FASTQ to clean dataset                                     |
|       |rmAdapter      |remove adapter of FASTQ                                           |
|       |reform         |edit the FASTQ file (reverse/complement)                          |
|       |Mul2Sin        |covert muti-lines FQ to single-line Fq                            |
|       |bubble         |filter regions with large number of N                             |
|       |changQ         |update the quality of FASTQ                                       |
|       |Fq2Fa          |covert FASTQ to FASTA                                             |

#### Formtools

|      Function                                                 |                      Description                  |
|:--------------------------------------------------------------|:--------------------------------------------------|
|                   CDS2Pep                                     |      convert CDS to Pep format                    |
|                   Soap2fq                                     |      convert SOAP to FASTQ format                 |
|                   Bam2Fq                                      |      convert BAM to FASTQ format                  |
|                   Soap2Sam                                    |      convert SOAP to Sam format                   |
|                   Xam2SOAP                                    |      convert BAM/SAM to SOAP format               |
|                   SF                                          | finding intersections or differences of two files |

### 3) Examples
------------
* 1) sort fa files
```
   # sort by seq length
   ./bin/FxTools  Fatools sort    -i  ref.fa      -s   length   -r    > ref.sort.fa
   #  sort by seq ID & gzip out 
   ./bin/FxTools  Fatools sort    -i  ref.fa      -s  name  -o  ref.sort.fa.gz
```

* 2) split the fa files
```
# split by one seq on file 
	./FxTools  Fatools split   -i   in.fa.gz    -o outDir/ -g
# split to fixed Number of subflie 
	./FxTools  Fatools split   -i   in.fa       -o outDir/  -f 12 
# split to fixed seq Number in one sub-file 
        ./FxTools  Fatools split   -i   in.fa       -o outDir/  -s 12 
```

* 3) Calculate qulity of fq files
```
# Give pdf of fastq  Base Q Distribute and stat result
	 ./FxTools  Fqtools   fqcheck  -i A_1.fq.gz    A_2.fq.gz  -o out1Prefix  out2Prefix 
# SE fqstq also can be supported 
	 ./FxTools  Fqtools   fqcheck  -i A.fq.gz   -o outPrefix 
```

* 4) change qulity of fq 
```
#  fstaq Q chang : by ASCII33-->ASCII64[+31] with ResetID & MaxQ:h
	./FxTools  Fqtools   changQ   -i in.fq.gz   -o out.fq  -s 4 
```

see more other Usage in the Documentation

### 4) Format
------------
Format Introduction
* [FASTA format](https://en.wikipedia.org/wiki/FASTA_format)
* [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format)
* [sam/bam format](https://samtools.github.io/hts-specs/SAMv1.pdf)
* [soap format](http://soap.genomics.org.cn/soapaligner.html)   
     [chinese soap introduction](http://blog.sina.com.cn/s/blog_70b2b6020101b609.html)


### 5）Discussing
------------
- [:email:](https://github.com/BGI-shenzhen/FxTools) hewm2008@gmail.com / hewm2008@qq.com / xuxiaomin@bgi.com
- join the<b><i> QQ Group : 125293663</b></i>


######################swimming in the sky and flying in the sea ########################### ##

