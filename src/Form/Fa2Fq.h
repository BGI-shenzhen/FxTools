#ifndef Fa2Fq_H_
#define Fa2Fq_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "../ALL/comm.h"
#include <cstdlib>
#include <gzstream.h>
#include <algorithm>
#include <zlib.h>
#include <stdio.h>
#include "../ALL/kseq.h"
#include "../ALL/DataClass.h"

using namespace std;

int  print_Fa2Fq()
{
	cout <<""
		"\n"
		"\tUsage: Fa2Fq -i <in.fa> -o <out.fq>\n"
		"\n"
		"\t\t-i     <str>   Input Fasta File\n"
		"\t\t-o     <str>   Output Fastq file\n"
		"\n"
		"\t\t-h             show this help\n" 
		"\n";
	return 1;
}


int usageFa2Fq(int argc, char **argv , In3str1v * para )
{
	if (argc <=2 ) {print_Fa2Fq();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InFa" || flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut" || flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para->InStr2=argv[i];
		}
		else if (flag  == "help" || flag  == "h")
		{
			print_Fa2Fq();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para->InStr1).empty()  ||  (para->InStr2).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	(para->InStr2)=add_Asuffix(para->InStr2) ;
	return 1 ;
}

int Form_Fa2Fq_main(int argc, char *argv[])
	//int main(int argc, char *argv[])
{
	In3str1v * para = new In3str1v;
	if( usageFa2Fq(argc, argv, para )==0)
	{
		delete  para ;
		return 1;
	}


	ogzstream OUT((para->InStr2).c_str()) ;
	if (!OUT.good())
	{
		cerr<<"Can't open OutFile "<<(para->InStr2)<<endl;
		return 1;
	}

	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen((para->InStr1).c_str(), "r");
	seq = kseq_init(fp);

	while ((l = kseq_read(seq)) >= 0)
	{
		string Ref=(seq->seq.s);
		llong seq_length=(seq->seq.l);
		string ID=seq->name.s;
		if (seq->comment.l)
		{
			string tmp=seq->comment.s;
			tmp.erase(remove_if(tmp.begin(),tmp.end(),bind2nd(equal_to<char>(),' ')),tmp.end());
			tmp.erase(remove_if(tmp.begin(),tmp.end(),bind2nd(equal_to<char>(),'\t')),tmp.end());
			ID=ID+"#"+tmp+"/1";
		}
		else
		{
//			ID=ID+"#A/1" ;
		}

		OUT<<"@"<<ID<<"\n"<<Ref<<"\n+"<<"\n";
		for (int i=0 ; i<seq_length ; i++)
		{
			OUT<<"h";
		}
		OUT<<"\n";
	}
	OUT.close();

	kseq_destroy(seq);
	gzclose(fp);

	delete para ;
	return 0;
}
#endif /// Fa2Fq_H_ 
///////// swimming in the sky and flying in the sea ////////////

