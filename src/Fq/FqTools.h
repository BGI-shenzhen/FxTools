#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <zlib.h>
#include "../ALL/kseq.h"
#include "../ALL/comm.h"

#include "FQ_Filter.h"
#include "FQ_Stat.h"
#include "FQ_ChangQ.h"
#include "FQ_RAD_GBS_SplitV2.h"
#include "FQ_RmAdapter.h"
#include "FQ_CutIndex.h"
#include "FQ_FilterV2.h"
#include "FQ_bubble.h"
#include "FQ_check.h"
#include "FQ_Split.h"
#include "FQ_Reform.h"
#include "FQ_Rand.h"
#include "FQ_RmDup.h"
#include "FQ_Valid.h"
#include "FQ_Mul2Sin.h"

//KSEQ_INIT(gzFile, gzread)


using namespace std;

int FQ_Valid_main( int argc,char *argv[] );
int FQ_RmDup_main(int argc, char **argv);
int FQ_Rand_main(int argc, char **argv);
int FQ_Split_main(int argc, char **argv);
int FQ_Stat_main( int argc,char *argv[] ) ;
int FQ_Filter_main(int argc, char **argv) ;
bool FQ_Pooling_main(int argc, char *argv[]) ;
bool FQ_Index_main(int argc, char *argv[]) ;
bool FQ_Se_main(int argc, char *argv[]) ;
int FQ_IndexCut_main(int argc, char **argv) ;
int FQ_FilterNotrim_main(int argc, char **argv);
int FQ_bubble_main(int argc, char **argv);
//int FQ_Split_Pooing_main(int argc, char *argv[]);
int FQ_ChangQ_main(int argc, char **argv) ;
int FQ_Check_main(int argc, char **argv) ;
int FQ_Reform_main(int argc, char *argv[]) ;
int FQ_RmAdpter_main(int argc, char **argv) ;
int FQ_Mul2Sin_main( int argc,char *argv[] ) ;

static int  FQ_usage ()
{
	cerr<<""
		"\n"
		"\tFqTools Usage:\n\n"
		"\t-- Summary Module:\n"
		"\t\tvalid          check fq and change to valid Fq\n"
		"\t\tstat           quick stat fastq's info\n"
		"\t\tfqcheck        fqcheck Base Q Distribute & adapter list\n"
		"\t-- Split Module:\n"
		"\t\tsplitpool      split pooling Fq to sample for RAD (GBS)\n"
		"\t\tsplitFq        split Fq to small File with max Read Num\n"
		"\t\tcut            cut the Read Length in the Fq\n"
		"\t\trand           rand out proportion of the Read in Fq\n"
//		"\t\tfilterV1       filter fastq for clean datas with trim\n"
		"\t-- Modify Module:\n"
		"\t\tfilter         filter fastq for clean datas select trim\n"
		"\t\trmdup          remove duplicatesthe by Block\n"
		"\t\tmul2sin        convert muti-lines FQ to single-line Fq\n"
		"\t\trmAdapter      remove adapter read of Fq\n"
		"\t\treform         reform/modify the Fq sequence\n"
		"\t\tbubble         filter the N bubble site Read\n"
		"\t\tchangeQ        change Fq seq Quality (+/- 31) and ID\n"
//		"\t\tFq2Fa          Fastq --> Fasta Format\n"
//		"\n"
//		"\t\tHelp           show more details for help\n"
		"\n";
	return 1;
}


int FQ_Tools_main(int argc, char *argv[])
{
	if (argc < 2) { return FQ_usage(); }
	else if (strcmp(argv[1], "stat") == 0) { return FQ_Stat_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "valid") == 0) { return FQ_Valid_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "fqcheck") == 0) { return FQ_Check_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "rmdup") == 0) { return FQ_RmDup_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Mul2Sin") == 0) { return FQ_Mul2Sin_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "mul2sin") == 0) { return FQ_Mul2Sin_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "filterV1") == 0) { return FQ_Filter_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "filterTrim") == 0) { return FQ_Filter_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "filterNoTrim") == 0) { return FQ_FilterNotrim_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "filter") == 0) { return FQ_FilterNotrim_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "rand") == 0) { return FQ_Rand_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "rmAdapter") == 0) { return FQ_RmAdpter_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "reform") == 0) { return FQ_Reform_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "cut") == 0) { return FQ_IndexCut_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "bubble") == 0) { return FQ_bubble_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "changQ") == 0) { return FQ_ChangQ_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "changeQ") == 0) { return FQ_ChangQ_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "splitpool") == 0) { return FQ_Split_Pooing_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "splitFq") == 0) { return FQ_Split_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Fq2Fa") == 0) { return    Form_Fq2Fa_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Help") == 0  || strcmp(argv[1], "help") == 0)  {  FQ_usage(); }
	else
	{
		cerr<<"FqTools [main] unrecognized command "<<argv[1]<<endl;
		return 1;
	}
	return 0;	
}


