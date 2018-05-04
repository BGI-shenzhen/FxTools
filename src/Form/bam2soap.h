#ifndef xam2soap_H_
#define xam2soap_H_
#include <iostream>
#include <fstream>
#include <gzstream.h>
#include "../ALL/comm.h"
#include "../ALL/DataClass.h"
using namespace std;


int  print_Ausage_A14()
{
	cout <<""
		"\n"
		"\tUsage: bam2Soap  -i <in.bam>  -o <Out.soap>\n"
		"\n"
		"\t\t-i    <str>   Input Sam/Bam Format file\n"
		"\t\t-o    <str>   Output the SoapFormat file\n"
		"\n"
		"\t\t-Q    <int>   SeqQ shift Trans:(-31/+31/0)[0]\n"
		"\t\t              Sanger<->Illumina:! <-> @ \n"
		"\t\t-h            show this help\n" 
		"\n";
	return 1;
}


int parse_Acmd_A14(int argc, char **argv, In3str1v * para_A14)
{
	if (argc <=2 ) {print_Ausage_A14();return 0 ;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InFile" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_A14->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_A14->InStr3=argv[i];
		}
		else if (flag  ==  "QShift" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_A14->InInt=atoi(argv[i]);
		}

		else if (flag  == "help")
		{
			print_Ausage_A14();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if ((para_A14->InStr3).empty()  || (para_A14->InStr1).empty())
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	(para_A14->InStr3)= add_Asuffix ( (para_A14->InStr3) ) ;

	return 1 ;
}


void Talignment_format(  bam_hdr_t *header  ,bam1_t *aln , ogzstream  &  OUT , uint8_t * Base ,int shiftQ )
{


	string readID=bam_get_qname(aln);

	string ab="a";
	if  ( (((aln)->core.flag&64) != 0))
	{
		readID=readID + "/1" ;
	}
	else if ((((aln)->core.flag&128) != 0) )
	{
		readID=readID + "/2" ;
		ab="b";
	}

	uint8_t  * seqQ=bam_get_qual(aln);
	uint8_t   *seq=bam_get_seq(aln);
	string BaseSeq="";
	string BaseSeqQ="";
	for(int i=0; i < (aln->core).l_qseq ; i++)
	{
		BaseSeqQ=BaseSeqQ+char(seqQ[i]+shiftQ);
		BaseSeq=BaseSeq+char(Base[bam_seqi(seq, i)]);
	}

	string ZF="+";
	if  ((((aln)->core.flag&16) != 0))
	{
		ZF="-";
	}



	string mapchr=(header->target_name[(aln->core).tid]);
	int site=(aln->core).pos;

	int Hit=1;
	if  ((aln->core).qual==0)
	{
		Hit=10;
	}
	else if   ((aln->core).qual<30)
	{
		Hit=2;
	}


	string cigarStr="*";
	if (aln->core.n_cigar)
	{ // cigar
		uint32_t *cigar = bam_get_cigar(aln);
		cigarStr="";
		for (int i = 0; i < aln->core.n_cigar; ++i) 
		{
			   stringstream   sstrm ;
	           sstrm  <<  (bam_cigar_oplen(cigar[i])) ;
		  	  cigarStr=cigarStr+sstrm.str()+(char)bam_cigar_opchr(cigar[i]);
		}
	}
	OUT<<readID<<"\t"<<BaseSeq<<"\t"<<BaseSeqQ<<"\t"<<Hit<<"\t"<<ab<<"\t"<< (aln->core).l_qseq<<"\t"<<ZF<<"\t"<<mapchr<<"\t"<<(aln->core).pos<<"\t"<<cigarStr<<endl;

}


int Xam2Soap_main(int argc, char **argv)
//int main(int argc, char **argv)
{

	In3str1v * para_A14 = new In3str1v;
	if ( parse_Acmd_A14(argc, argv,para_A14)==0 )
	{
		delete  para_A14 ;
		return 1;
	}

	(para_A14->InInt)+=31;

	bam_hdr_t *header;
	bam1_t *aln = bam_init1();
	samFile *in = sam_open((para_A14->InStr1).c_str(), "r");
	header = sam_hdr_read(in);
	uint8_t Base[16] = {0,65,67,0,71,0,0,0,84,0,0,0,0,0,0,78};

	ogzstream OUT(para_A14->InStr3.c_str());
	if (!OUT.good())
	{
		cerr<<"Can't open OutFile "<<para_A14->InStr3<<endl;
		return 1;
	}



	while  ( (sam_read1(in, header, aln) >= 0)  )
	{
		if ((aln->core).tid < 0)
		{
			continue ;
		}
		Talignment_format(  header , aln ,   OUT , Base,para_A14->InInt);
	}


	OUT.close();
	sam_close(in);
	bam_destroy1(aln);
	bam_hdr_destroy(header);


	delete para_A14 ;
	return 0;
}
#endif // xam2soap_H_
///////// swimming in the sky and flying in the sea ////////////

/* format the sam text to the soap text*/


