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
#include <zlib.h>
#include <stdio.h>
#include "../ALL/kseq.h"
#include "../ALL/DataClass.h"

//KSEQ_AINIT(gzFile, gzread)

using namespace std;

int  print_Headerusage()
{
	cout <<""
		"\n"
		"\tUsage: dict  -i <in.fa> \n"
		"\n"
		"\t\t-i    <str>   Input Ref.fa to new Ref.dict Samheader\n"
		"\n"
		"\t\t-o    <str>   Output file,otherwise[STDOUT]\n"
		"\n"
		"\t\t-h            show this help\n" 
		"\n";
	return 1;
}


int parse_HeaderCmd(int argc, char **argv , In3str1v * para_A1 )
{
	if (argc <=2 ) {print_Headerusage();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InPut" ||  flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_A1->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut"   ||  flag  == "o" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_A1->InStr2=argv[i];
		}
		else if (flag  == "help" ||  flag  == "h")
		{
			print_Headerusage();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para_A1->InStr1).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	return 1 ;
}




int FA_SamHeader_main(int argc, char *argv[])
{
	In3str1v * para_A1 = new In3str1v;
	if( parse_HeaderCmd(argc, argv, para_A1 )==0)
	{
		delete  para_A1 ;
		return 1;
	}
	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen((para_A1->InStr1).c_str(), "r");
	seq = kseq_init(fp);
	string chr_Ainfo=(para_A1->InStr1)+".dict";
	if ((para_A1->InStr2).empty())
	{
		if ( access(chr_Ainfo.c_str(), 0) == 0 )
		{
			while ((l = kseq_read(seq)) >= 0) 
			{
				cout<<"@SQ\tSN:"<<(seq->name.s)<<"\tLN:"<<(seq->seq.l)<<endl;
			}
		}
		else
		{
			ofstream  OUT (chr_Ainfo.c_str());
			if(!OUT.good())
			{
				cerr << "open OUT File error: "<<chr_Ainfo<<endl;
			} 
			while ((l = kseq_read(seq)) >= 0)
			{
				cout<<"@SQ\tSN:"<<(seq->name.s)<<"\tLN:"<<(seq->seq.l)<<endl;
				OUT<<"@SQ\tSN:"<<(seq->name.s)<<"\tLN:"<<(seq->seq.l)<<endl;
			}
			OUT.close();
		}
	}
	else
	{
		if ( chr_Ainfo == (para_A1->InStr2) || ( access(chr_Ainfo.c_str(), 0) == 0 ) )
		{
			ofstream  OUT ((para_A1->InStr2).c_str());
			if(!OUT.good())
			{
				cerr << "open OUT File error: "<<(para_A1->InStr2)<<endl;
				delete para_A1 ; return 1 ;
			}
			while ((l = kseq_read(seq)) >= 0) 
			{
				OUT<<"@SQ\tSN:"<<(seq->name.s)<<"\tLN:"<<(seq->seq.l)<<endl;
			}
			OUT.close();
		}
		else
		{
			ofstream  OUT2 ((para_A1->InStr2).c_str());
			if(!OUT2.good())
			{
				cerr << "open OUT File error: "<<(para_A1->InStr2)<<endl;
				delete para_A1 ; return 1 ;
			}
			ofstream  OUT (chr_Ainfo.c_str());
			if(!OUT.good())
			{
				cerr << "open OUT File error: "<<chr_Ainfo<<endl;
			}

			while ((l = kseq_read(seq)) >= 0) 
			{
				if (OUT.good())
				{
					OUT<<"@SQ\tSN:"<<(seq->name.s)<<"\tLN:"<<(seq->seq.l)<<endl;
				}
				OUT2<<"@SQ\tSN:"<<(seq->name.s)<<"\tLN:"<<(seq->seq.l)<<endl;
			}
			OUT.close();
			OUT2.close();
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	delete para_A1 ;
	return 0;

}

///////// swimming in the sky and flying in the sea ////////////
