#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "../ALL/comm.h"
#include <htslib/sam.h>
#include <cstdlib>
#include <gzstream.h>
using namespace std;
typedef long long llong ;



int  print_Ausage_Xam2fq()
{
	cout <<""
		"\n"
		"\tUsage: Bam2fq -i <in.bam> -o <out.fq>\n"
		"\n"
		"\t\t-i     <str>   InPut sam/bam File\n"
		"\t\t-o     <str>   OutPut fq File\n"
		"\n"
		"\t\t-u             only output UnMapp read [NA]\n"
		"\n"
		"\t\t-h             show this help\n" 
		"\n";
	//"\t\t-ShiftQ    <int>   Final phred quality [+31/0/-31]\n" 
	return 1;
}

int parse_Acmd_Xam2fq(int argc, char **argv, In3str1v * para_Xam2fq)
{
	if (argc <=2 ) {print_Ausage_Xam2fq();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","") ;

		if (flag  == "InPut"  || flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Xam2fq->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut"   || flag  == "o" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0; }
			i++;
			para_Xam2fq->InStr2=argv[i];
		}
		else if (flag  ==  "UnMap"  || flag  == "u")
		{
			para_Xam2fq->TF=false;
		}
		///*////
		else if (flag  == "help")
		{
			print_Ausage_Xam2fq();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para_Xam2fq->InStr2).empty() ||   (para_Xam2fq->InStr1).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;    
	}
	para_Xam2fq->InStr2=add_Asuffix(para_Xam2fq->InStr2);
	return 1 ;
}



//programme entry
///////// swimming in the sky and flying in the sea ////////////
//int main(int argc, char **argv)
int Xam2fq_main(int argc, char **argv)
{
	In3str1v * para_Xam2fq = new In3str1v;
	if (parse_Acmd_Xam2fq(argc, argv , para_Xam2fq )==0)
	{
		delete para_Xam2fq ; 
		return 1;
	}

	bam_hdr_t *header;
	bam1_t *aln = bam_init1();
	samFile *in = sam_open((para_Xam2fq->InStr1).c_str(), "r");
	header = sam_hdr_read(in);


	ogzstream OUT(para_Xam2fq->InStr2.c_str()) ;


	int Samp[256]={0};
	Samp['A']='T'; Samp['C']='G';
	Samp['T']='A'; Samp['G']='C';
	Samp['N']='N';

	int  Count=0;
	uint8_t Base[16] = {0,65,67,0,71,0,0,0,84,0,0,0,0,0,0,78};
	if ( para_Xam2fq->TF)
	{
		while  ( (sam_read1(in, header, aln) >= 0)    )
		{
			string readID=bam_get_qname(aln);

			if  ( (((aln)->core.flag&64) != 0))
			{
				readID= "@"+readID + "/1" ;
			}
			else if ((((aln)->core.flag&128) != 0) )
			{
				readID= "@"+readID + "/2" ;
			}
			else
			{
				readID= "@"+readID ;
			}


			uint8_t  * seqQ=bam_get_qual(aln);
			uint8_t   *seq=bam_get_seq(aln);
			string BaseSeq="";
			string BaseSeqQ="";
			for(int i=0; i < (aln->core).l_qseq ; i++)
			{
				BaseSeqQ=BaseSeqQ+(char)((seqQ[i])+31);
				BaseSeq=BaseSeq+(char)(Base[bam_seqi(seq, i)]);
			}
	



			if  ((((aln)->core.flag&16) != 0))
			{
				reverse(BaseSeq.begin(), BaseSeq.end());
				reverse(BaseSeqQ.begin(), BaseSeqQ.end());
				for (int i=0 ; i<(aln->core).l_qseq ; i++)
				{
					BaseSeq[i]=Samp[BaseSeq[i]];
				}
			}

			OUT<<readID<<"\n"<<BaseSeq<<"\n+\n"<<BaseSeqQ<<"\n";


		}




	}



	else
	{



		while  ( (sam_read1(in, header, aln) >= 0)    )
		{
			if ((aln->core).tid >= 0)
			{
				continue ;
			}

			string readID=bam_get_qname(aln);

			if  ( (((aln)->core.flag&64) != 0))
			{
				readID= "@"+readID + "/1" ;
			}
			else if ((((aln)->core.flag&128) != 0) )
			{
				readID= "@"+readID + "/2" ;
			}
			else
			{
				readID= "@"+readID ;
			}



			uint8_t  * seqQ=bam_get_qual(aln);
			uint8_t   *seq=bam_get_seq(aln);
			string BaseSeq="";
			string BaseSeqQ="";
			for(int i=0; i < (aln->core).l_qseq ; i++)
			{
				BaseSeqQ=BaseSeqQ+Int2Str(seqQ[i]);
				BaseSeq=BaseSeq+Int2Str(Base[bam_seqi(seq, i)]);
			}




			if  ((((aln)->core.flag&16) != 0))
			{
				reverse(BaseSeq.begin(), BaseSeq.end());
				reverse(BaseSeqQ.begin(), BaseSeqQ.end());
				for (int i=0 ; i<(aln->core).l_qseq ; i++)
				{
					BaseSeq[i]=Samp[BaseSeq[i]];
				}
			}

			OUT<<readID<<"\n"<<BaseSeq<<"\n+\n"<<BaseSeqQ<<"\n";


		}









	}
	sam_close(in);
	bam_destroy1(aln);
	bam_hdr_destroy(header);
	OUT.close();
	delete para_Xam2fq ;
	return 0;
}

///////// swimming in the sky and flying in the sea ////////////
