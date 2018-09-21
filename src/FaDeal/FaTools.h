#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <zlib.h>
#include "../ALL/kseq.h"
#include "../ALL/comm.h"

#include "RefMatch.h"
#include "BaseModify.h"
#include "RefGetCDS.h"
#include "CDS2Pep.h"
#include "RefFindN.h"
#include "RefSubSeq.h"
#include "RefSample.h"
#include "RefSort.h"
#include "RefReform.h"
#include "RefSplit.h"
#include "RefHeader.h"
#include "ChangPosi.h"
#include "ReffindSubSeq.h"
#include "RefRegenerate.h"
#include "RefStat.h"
#include "RefRand.h"
#include "RefFilter.h"

//KSEQ_INIT(gzFile, gzread)

using namespace std;
int FA_Rand_main(int argc, char *argv[]);
int FA_BaseModify_main(int argc, char *argv[]) ;
int FA_FSubSeq_main(int argc, char *argv[]);
int FA_macth_main(int argc, char *argv[]) ;
int FA_NRegion_main(int argc, char *argv[]) ;
int FA_Order_main(int argc, char *argv[]) ;
int FA_stat_main(int argc, char *argv[]);
int FA_SubSeq_main(int argc, char *argv[]) ;
int FA_Filter_main(int argc, char *argv[]) ;
int FA_Regenerate_main(int argc, char *argv[]);
int FA_Reform_main(int argc, char *argv[]) ;
int FA_GetCDS_main(int argc, char *argv[]) ;
int FA_CDS2Pep_main(int argc, char *argv[]) ;
int FA_Sort_main(int argc, char *argv[]) ;
int FA_Split_main(int argc, char *argv[]) ;
int FA_SamHeader_main(int argc, char *argv[]) ;
int FA_ChangPosi_main(int argc, char *argv[]) ;
//int FA_Cut_main(int argc, char *argv[]) ;

static int  FA_usage ()
{
	cerr<<""
		"\n"
		"\tFaTools Usage:\n\n"
		"\t-- summary module:\n"
		"\t\tstat        quick stat fasta's Length GC N50 etc info\n"
		"\t\tdict        quick give out Ref.dict,the Samheader of fa\n"
		"\t-- split module:\n"
//		"\t\tsplit       split InFa, each Seq one File\n"
//		"\t\tcutF        cut fasta to fixed Num of subFile in total\n"
		"\t\tsplit       split InFa to subFile by fixed SeqNum or FileNum\n"
		"\t\trand        rand out proportion of the seq in Fa\n"
		"\t-- search module:\n"
		"\t\tfindN       quick find N region in fasta\n"
		"\t\tlocate      quick find fasta's one SubSeq region\n"
		//"\t\tfindSubSeq  quick find fasta's one SubSeq region\n"
		"\t\tgrep        search for the subsequence by target region\n"
		"\t\textractP    get Seq by specified ID or pattern\n"
		"\t\textractN    get Seq by specified order Num range\n"
		"\t\tgetCdsPep   find CDS & Pep sequences base on the GFF file\n"
		"\t\tsort        rank the seq By SeqID or Length\n"
		"\t-- modify module:\n"
		"\t\tfilter      filter the short & too many 'N' Seq\n"
		"\t\treform      reform/modify the sequence\n"
		"\t\tBaseModify  modify the muti single-base in seq\n"
		"\t\tJoinSca     joining scaffolds into pseudo chromosomes\n"
		"\t\tchangePosi  Change Position back to Scaffolds(regenerate)\n"
//		"\t\tFa2Fq       Fasta --> Fastq Format\n"
//		"\n"        
//		"\t\tHelp        show more details for help\n"
		"\n";
	return 1;
}


int Fa_Tools_main(int argc, char *argv[])
{
	if (argc < 2) { return FA_usage(); }
	else if (strcmp(argv[1], "findN") == 0) { return FA_NRegion_main(argc-1, argv+1) ; }
	//else if (strcmp(argv[1], "split") == 0) { return FA_Split_main(argc-1, argv+1) ;}
	else if (strcmp(argv[1], "stat") == 0) { return FA_stat_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "split") == 0) { return FA_Split_main(argc-1, argv+1) ;}
	else if (strcmp(argv[1], "findSubSeq") == 0) { return FA_FSubSeq_main(argc-1, argv+1) ;}
	else if (strcmp(argv[1], "locate") == 0) { return FA_FSubSeq_main(argc-1, argv+1) ;}
//	else if (strcmp(argv[1], "cutS") == 0) { return FA_CutSeq_main(argc-1, argv+1) ;}
	else if (strcmp(argv[1], "extractP") == 0) { return FA_macth_main(argc-1, argv+1) ;}
	else if (strcmp(argv[1], "extractN") == 0) { return FA_Order_main(argc-1, argv+1) ;}
	else if (strcmp(argv[1], "grep") == 0) { return FA_SubSeq_main(argc-1, argv+1); }
	else if (strcmp(argv[1], "rand") == 0) { return FA_Rand_main(argc-1, argv+1); }
	else if (strcmp(argv[1], "filter") == 0) { return FA_Filter_main(argc-1, argv+1);}
	else if (strcmp(argv[1], "reform") == 0) { return FA_Reform_main(argc-1, argv+1);}
	else if (strcmp(argv[1], "regenerate") == 0) { return FA_Regenerate_main(argc-1,argv+1) ;}
	else if (strcmp(argv[1], "JoinSca") == 0) { return FA_Regenerate_main(argc-1,argv+1) ;}
	else if (strcmp(argv[1], "getCdsPep") == 0) { return FA_GetCDS_main(argc-1, argv+1) ;}
	else if (strcmp(argv[1], "CDS2Pep") == 0) { return FA_CDS2Pep_main(argc-1, argv+1); }
	else if (strcmp(argv[1], "sort") == 0) { return  FA_Sort_main(argc-1, argv+1);}
	else if (strcmp(argv[1], "dict") == 0) { return  FA_SamHeader_main(argc-1, argv+1);}
	else if (strcmp(argv[1], "BaseModify") == 0) { return FA_BaseModify_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "ChangPosi") == 0) { return FA_ChangPosi_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "changePosi") == 0) { return FA_ChangPosi_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Fa2Fq") == 0) { return    Form_Fa2Fq_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Help") == 0  || strcmp(argv[1], "help") == 0) { FA_usage(); }
	else
	{
		cerr<<"FaTools [main] unrecognized command "<<argv[1]<<endl;
		return 1;
	}
	return 0;	
}


