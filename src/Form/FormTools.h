#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <zlib.h>
#include "../ALL/kseq.h"
#include "../ALL/comm.h"

#include "Fq2Fa.h"
#include "FileSF.h"
#include "Fa2Fq.h"
#include "Soap2Fq.h"
#include "Bam2fq.h"
#include "bam2soap.h"
#include "soap2bam.h"
#include "FileMerge.h"

using namespace std;
int Xam2fq_main(int argc, char **argv) ;
int Form_Fq2Fa_main(int argc, char *argv[]) ;
int FA_CDS2Pep_main(int argc, char *argv[]) ;
int Form_Fa2Fq_main(int argc, char *argv[]) ;
int Soap2fq_main(int argc, char **argv) ;
int Soap2Xam_main(int argc, char *argv[]) ;
int Xam2Soap_main(int argc, char **argv);
int File_SF_main(int argc,char *argv[]);
int Other_Merge_main(int argc,char *argv[]);


static int  FormTools_usage ()
{
	cerr<<""
		"\n"
		"\tForm Tools Usage:\n\n"
		"\t\tCDS2Pep      CDS seq --> Pep Format\n"
		"\t\tBam2fq       Alignment bam/sam --> fq Format\n"
		"\t\tSoap2fq      Alignment Soap --> fq Format\n"
		"\t\tSoap2Bam     Alignment Soap -->  Bam/Sam\n"
		"\t\tBam2Soap     Alignment Bam/Sam --> Soap Format\n"
		"\t\tFq2Fa        Fastq --> Fasta Format\n"
		"\t\tFa2Fq        Fasta --> Fastq Format\n"
		"\t\tSF           same and Diff(commond &shuffle) two file\n"
		"\t\tMerge        Meger muti-soft file to a big sort File\n"

		"\n"
		"\t\tHelp         Show this help\n"
		"\n";
	return 1;
}


int Form_Tools_main(int argc, char *argv[])
{
	if (argc < 2) { return FormTools_usage(); }
	else if (strcmp(argv[1], "CDS2Pep") == 0) { return  FA_CDS2Pep_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Soap2fq") == 0) { return  Soap2fq_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Soap2Bam") == 0) { return  Soap2Xam_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "bam2fq") == 0) { return  Xam2fq_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Fq2Fa") == 0) { return  Form_Fq2Fa_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Fa2Fq") == 0) { return  Form_Fa2Fq_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Bam2fq") == 0) { return  Xam2fq_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "SF") == 0) { return  File_SF_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Merge") == 0) { return  Other_Merge_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Xam2Soap") == 0) { return Xam2Soap_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Bam2Soap") == 0) { return Xam2Soap_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Help") == 0  || strcmp(argv[1], "help") == 0)  {  FormTools_usage(); }
	else
	{
		cerr<<"Form Tools [main] unrecognized command "<<argv[1]<<endl;
		return 1;
	}
	return 0;	
}


