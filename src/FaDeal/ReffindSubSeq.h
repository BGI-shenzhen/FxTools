#ifndef RefSubSeq_h_
#define RefSubSeq_h_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "../ALL/comm.h"
#include "../ALL/DataClass.h"
#include <cstdlib>
#include "../ALL/gzstream.h"
#include <algorithm>
#include <zlib.h>
#include <stdio.h>
#include "../ALL/kseq.h"

using namespace std;

int  print_AusageAt29()
{
	cout <<""
		"\n"
		"\tUsage: findSubSeq  -i <in.fa> -f <seq.fa> \n"
		"\n"
		"\t\t-i    <str>   input FASTA\n"
		"\t\t-o    <str>   output file, default [STDOUT]\n"
		"\t\t-f    <str>   input subsequences, should be formatted as FASTA\n"
		"\n"
		"\t\t-h            show more details for help\n" 
		"\n";
	return 1;
}

void More_HelpFA_5 ()
{
	cout<<"\n\
\n\
\t\t1.  findSubSeq  -i <in.fa> -f <seq.fa>\n\
\t\tthis would give the locations of subsequences found on the input FASTA. \n\
\t\tThe subsequence file stores all sequences user would like to locate on the input FASTA and should be formatted as FASTA. For example, if user would like to locate two subsequences (AATT and CCGG) on input FASTA, the subsequence file could be formatted below:\n\
\t\t>seq1\n\
\t\tAATT\n\
\t\t>seq2\n\
\t\tCCGG\n\
\t\tIf seq1 locates on chromosomeX from base site x1 to x2, and seq2 locates on chromosomeY from base site y1 to y2, then the output would be shown as:\n\
\t\tchromosomeX  x1  x2  seq1\n\
\t\tchromosomeY  y1  y2  seq2\n\
\t\tIf any sequence in the subsequence was not found in the input FASTA, it would not show up in the output.\n\
\n\
\n";
}


int parse_AcmdAt29(int argc, char **argv , In3str1v * paraAt29 )
{
	if (argc <=1 ) {print_AusageAt29();return 0;}

	for(int i = 1; i < argc  ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "AllSeq"  ||  flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg( flag ) ; return 0;}
			i++;
			paraAt29->InStr1=argv[i];
		}
		else if (flag  ==  "SubSeq"   ||  flag  == "f")
		{
			if(i + 1 == argc) {LogLackArg( flag ) ;return 0;}
			i++;
			paraAt29->InStr3=argv[i];
		}
		else if (flag  ==  "OutPut"  ||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg( flag ) ;return 0;}
			i++;
			paraAt29->InStr2=argv[i];
		}
		else if (flag  == "help"  ||  flag  == "h")
		{
			More_HelpFA_5();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((paraAt29->InStr1).empty()  || (paraAt29->InStr3).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	return 1 ;
}




//int main(int argc, char *argv[])
int FA_FSubSeq_main(int argc, char *argv[])
{
	In3str1v * paraAt29 = new In3str1v;
	if( parse_AcmdAt29(argc, argv, paraAt29 )==0)
	{
		delete  paraAt29 ;
		return 1;
	}

	gzFile fpV2;
	kseq_t *seqV2;
	fpV2 = gzopen((paraAt29->InStr3).c_str(), "r");
	seqV2=kseq_init(fpV2);
	int l;
	string subseq="";
	string subname="";

	map <string ,string> ChrSeq ;

	while ((l = kseq_read(seqV2)) >= 0)
	{
		subseq=(seqV2->seq.s) ;
		transform(subseq.begin(),subseq.end(),subseq.begin(),::toupper) ;
		subname=(seqV2->name.s);
		ChrSeq[subname]=subseq;
	}
	kseq_destroy(seqV2);
	gzclose(fpV2);

	if (ChrSeq.empty())
	{
		cerr<<"subSeq empty"<<endl;
		delete paraAt29 ;
		return 1;
	}
	int subseq_len=subseq.length();
	gzFile fp;
	kseq_t *seq;
	fp = gzopen((paraAt29->InStr1).c_str(), "r");
	seq = kseq_init(fp);

	std::ostream * OUT = &cout;
	std::ofstream fout;
	if (!(paraAt29->InStr2).empty())
	{
		fout.open((paraAt29->InStr2).c_str());
		OUT = &fout;
	}

	while ((l = kseq_read(seq)) >= 0) 
	{
		string Ref=(seq->seq.s) ;
		transform(Ref.begin(),Ref.end(),Ref.begin(),::toupper) ;

		map <string ,string > ::iterator itL = ChrSeq.begin();
		for ( itL=ChrSeq.begin(); itL!=ChrSeq.end(); itL++)
		{
			subseq=itL->second;
			subname=itL->first;
			size_t N_tail_position = 0 ;
			size_t N_head_position=Ref.find(subseq, 0 );
			while( N_head_position != string::npos)
			{
				N_tail_position=N_head_position+subseq_len;
				(*OUT)<<(seq->name.s)<<"\t"<<N_head_position+1<<"\t"<<N_tail_position<<"\t"<<subname<<endl;
				N_head_position=Ref.find(subseq,N_tail_position);
			}
		}
	}

	kseq_destroy(seq);
	gzclose(fp);
	delete paraAt29 ;
	return 0;

}

#endif // SubSeq
///////// swimming in the sky and flying in the sea ////////////
