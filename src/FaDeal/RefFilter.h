#ifndef RefFilter_H_
#define RefFilter_H_

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

//KSEQ_AINIT(gzFile, gzread)

using namespace std;
typedef long long  llong ;

int  usage_FA03()
{
	cout <<""
		"\n"
		"\tUsage: filter  -i <in.fa> -o <out.fa>\n"
		"\n"
		"\t\t-i     <str>   Input fa for filter\n"
		"\t\t-o     <str>   Output file\n"
		"\n"
		"\t\t-l     <int>   The MinLength for Seq to pass[1000]\n"
		"\t\t-n   <float>   The max ratio of miss N Length[0.5]\n"
		"\n"
		"\t\t-h             show this help\n" 
		"\n";
	return 1;
}

int parse_Acmd_FA03(int argc, char **argv , ParaClass * para_FA03 )
{
	if (argc <=2 ) {usage_FA03();return 0;}


	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InPut" || flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FA03->InPut1=argv[i];
		}
		else if (flag  ==  "MinLen"  || flag  == "l")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FA03->InInt1=atoi(argv[i]);
		}
		else if (flag  ==  "NRatio" || flag  == "n")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FA03->Infloat1=atof(argv[i]);
		}
		else if (flag  ==  "OutPut"  || flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FA03->OutPut1=argv[i];
		}
		else if (flag  == "help"  || flag  == "h")
		{
			usage_FA03();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para_FA03->InPut1).empty() || (para_FA03->OutPut1).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	(para_FA03->OutPut1)=add_Asuffix((para_FA03->OutPut1)) ;
	return 1 ;
}

//int main(int argc, char *argv[])
int FA_Filter_main(int argc, char *argv[])
{
	ParaClass * para_FA03 = new ParaClass ;
	para_FA03->InInt1=1000;
	if( parse_Acmd_FA03(argc, argv, para_FA03 )==0)
	{
		delete  para_FA03 ;
		return 1;
	}
	int linecut= FaCutLine ((para_FA03->InPut1));

	ogzstream  OUT ((para_FA03->OutPut1).c_str());
	if(!OUT.good())
	{
		cerr << "open OUT File error: "<<(para_FA03->OutPut1)<<endl;
		delete  para_FA03 ;  return 1;
	}


	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen((para_FA03->InPut1).c_str(), "r");
	seq = kseq_init(fp);

	while ((l = kseq_read(seq)) >= 0) 
	{
		string Ref=(seq->seq.s);
		llong seq_length=(seq->seq.l);
		if (seq_length<para_FA03->InInt1)
		{
			continue ;
		}
		llong N_length=0;
		size_t N_tail_position = 0 ;
		size_t N_head_position=Ref.find_first_of( "Nn", 0 );
		while( N_head_position != string::npos)
		{
			N_tail_position=Ref.find_first_of("AaCcTtGg",N_head_position+1);
			if (N_tail_position!= string::npos )
			{
				N_length+=(N_tail_position-N_head_position);
			}
			else
			{
				N_length+=(seq_length-N_head_position);
				break ;
			}
			N_head_position=Ref.find_first_of("Nn",N_tail_position);
		}
		if (((N_length*1.0)/seq_length)<(para_FA03->Infloat1))
		{
			string ID=seq->name.s;
			if (seq->comment.l)
			{
				string tmp=seq->comment.s;
				ID=ID+"\t"+tmp;
			}
			Display( Ref,  ID , OUT ,linecut);
		}
	} 
	OUT.close();

	kseq_destroy(seq);
	gzclose(fp);
	delete para_FA03 ;
	return 0;

}
#endif
///////// swimming in the sky and flying in the sea ////////////
