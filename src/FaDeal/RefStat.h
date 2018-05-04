#ifndef RefStat_H_
#define RefStat_H_

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
#include <gzstream.h>
#include <zlib.h>
#include <stdio.h>
#include "../ALL/kseq.h"
#include "RefN50.h"


using namespace std;
typedef long long  llong ;

int  print_Ausage_A1()
{
	cout <<""
		"\n"
		"\tUsage: stat  -i <in.fa> \n"
		"\n"
		"\t\t-i         <str>   Input fa for stat the length\n"
		"\n"
		"\t\t-o         <str>   Output file,otherwise[STDOUT]\n"
		"\n"
		"\t\t-s                 To Stat N50 of fa\n"
		"\t\t-c         <int>   Pre-N50Stat Cutoff length < X scaffold[100]\n"
		"\n"
		"\t\t-h                 show more detailed help infomation\n"
		"\n";
	return 1;
}

int  print_Ausage_A1_detail()
{
 	cout<<"\n\
	stat  -InPut  <in.fa> -OutPut <out.fa>\n  \
		##      out.fa.gz would showd as :\n\
           #chr len N_len eff_len GC_count Rep_count\n\
 \n\
    	stat  -InPut  <in.fa> -OutPut <out.fa> -N50Stat\n\
           ##      out.fa would showd as:\n\
                    #scafold report#\n\
                    #contig report#\n\
                  #length distribution#\n\
\n";
 return 1;

}


int parse_Acmd_A1(int argc, char **argv , In3str1v * para_A1 )
{
	if (argc <=1 ) {print_Ausage_A1();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InPut" ||  flag  == "i" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_A1->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut" ||   flag  == "o" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_A1->InStr2=argv[i];
		}
		else if (flag  == "N50Stat" ||   flag  == "s"  )
		{
			para_A1->TF=false ;
		}
		else if (flag  == "CutLen"  ||   flag  == "c" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_A1->InInt=atoi(argv[i]);
		}

		else if (flag  == "help"  ||  flag  == "h")
		{
			print_Ausage_A1_detail();return 0;
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




int FA_stat_main(int argc, char *argv[])
{
	In3str1v * para_A1 = new In3str1v;
	para_A1->InInt=100;
	if( parse_Acmd_A1(argc, argv, para_A1 )==0)
	{
		delete  para_A1 ;
		return 1;
	}

	if (!( para_A1->TF) )
	{
		FaStatN50 ( para_A1->InStr2,  para_A1->InStr1 , para_A1->InInt );
		delete  para_A1 ;
		return 0;
	}

	igzstream  INID (( para_A1->InStr1).c_str(),ifstream::in);
	if(!INID.good())
	{
		cerr << "open InputFile error: "<<(para_A1->InStr1)<<endl;
		delete para_A1 ; return 1;
	}
	string chr_Ainfo=(para_A1->InStr1)+".chrlist";

	if ((para_A1->InStr2).empty())
	{
		if ( access(chr_Ainfo.c_str(), 0) == 0 )
		{
			cout<<"#chr\tLen\tN_Len\tEff_Len\tGC_Count\tRep_Count\n";
			string seqAA;
			getline(INID, seqAA ,'>');
			while(!INID.eof())
			{
				string line,ID;
				getline(INID, line,'\n');
				if (line.empty())
				{
					continue;
				}
				getline(INID, seqAA ,'>');
				istringstream isone (line,istringstream::in);
				isone>>ID;
				llong Ascii[256] = {0};
				llong ALength=seqAA.length();
				stat_str_base(seqAA , Ascii,ALength);
				llong Length=ALength-Ascii['\n'];
				llong N_Aleng= Ascii['N']+Ascii['n'] ;				
				llong eff_Acout=Length-N_Aleng ;
				llong GC_Aleng=Ascii['C']+Ascii['G'] + Ascii['c']+Ascii['g'] ;
				llong samllleter= Ascii['c']+Ascii['g'] + Ascii['t']+Ascii['a'];
				cout<<ID<<"\t"<<Length<<"\t"<<N_Aleng<<"\t"<<eff_Acout<<"\t"<<GC_Aleng<<"\t"<<samllleter<<endl;
			}
		}
		else
		{


			cout<<"#chr\tLen\tN_Len\tEff_Len\tGC_Count\tRep_Count\n";
			ofstream  OUT (chr_Ainfo.c_str());
			if(!OUT.good())
			{
				cerr << "open OUT File error: "<<chr_Ainfo<<endl;
			} 
			OUT<<"#chr\tLen\tN_Len\tEff_Len\tGC_Count\tRep_Count\n";

			string seqAA;
			getline(INID, seqAA ,'>');
			while(!INID.eof())
			{
				string line,ID;
				getline(INID, line,'\n');
				if (line.empty())
				{
					continue;
				}
				getline(INID, seqAA ,'>');
				istringstream isone (line,istringstream::in);
				isone>>ID;
				llong Ascii[256] = {0};

				llong ALength=seqAA.length();
				stat_str_base(seqAA , Ascii,ALength );
				llong Length=ALength-Ascii['\n'];
				llong N_Aleng= Ascii['N']+Ascii['n'] ;				
				llong eff_Acout=Length-N_Aleng ;
				llong GC_Aleng=Ascii['C']+Ascii['G'] + Ascii['c']+Ascii['g'] ;
				llong samllleter= Ascii['c']+Ascii['g'] + Ascii['t']+Ascii['a'];

				if (OUT.good())
				{
					OUT<<ID<<"\t"<<Length<<"\t"<<N_Aleng<<"\t"<<eff_Acout<<"\t"<<GC_Aleng<<"\t"<<samllleter<<endl;
				}
				cout<<ID<<"\t"<<Length<<"\t"<<N_Aleng<<"\t"<<eff_Acout<<"\t"<<GC_Aleng<<"\t"<<samllleter<<endl;
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
			OUT<<"#chr\tLen\tN_Len\tEff_Len\tGC_Count\tRep_Count\n";


			string seqAA;
			getline(INID, seqAA ,'>');
			while(!INID.eof())
			{
				string line,ID;
				getline(INID, line,'\n');
				if (line.empty())
				{
					continue;
				}
				getline(INID, seqAA ,'>');
				istringstream isone (line,istringstream::in);
				isone>>ID;

				llong Ascii[256] = {0};
				llong ALength=seqAA.length();
				stat_str_base(seqAA , Ascii,ALength );
				llong Length=ALength-Ascii['\n'];
				llong N_Aleng= Ascii['N']+Ascii['n'] ;				
				llong eff_Acout=Length-N_Aleng ;
				llong GC_Aleng=Ascii['C']+Ascii['G'] + Ascii['c']+Ascii['g'] ;
				llong samllleter= Ascii['c']+Ascii['g'] + Ascii['t']+Ascii['a'];
				OUT<<ID<<"\t"<<Length<<"\t"<<N_Aleng<<"\t"<<eff_Acout<<"\t"<<GC_Aleng<<"\t"<<samllleter<<endl;
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
			OUT2<<"#chr\tLen\tN_Len\tEff_Len\tGC_Count\tRep_Count\n";
			ofstream  OUT (chr_Ainfo.c_str());
			if(!OUT.good())
			{
				cerr << "open OUT File error: "<<chr_Ainfo<<endl;
			}
			OUT<<"#chr\tLen\tN_Len\tEff_Len\tGC_Count\tRep_Count\n";

			string seqAA;
			getline(INID, seqAA ,'>');
			while(!INID.eof())
			{
				string line,ID;
				getline(INID, line,'\n');
				if (line.empty())
				{
					continue;
				}
				getline(INID, seqAA ,'>');
				istringstream isone (line,istringstream::in);
				isone>>ID;

				llong Ascii[256] = {0};
				llong ALength=seqAA.length();
				stat_str_base(seqAA , Ascii,ALength );
				llong Length=ALength-Ascii['\n'];
				llong N_Aleng= Ascii['N']+Ascii['n'] ;				
				llong eff_Acout=Length-N_Aleng ;
				llong GC_Aleng=Ascii['C']+Ascii['G'] + Ascii['c']+Ascii['g'] ;
				llong samllleter= Ascii['c']+Ascii['g'] + Ascii['t']+Ascii['a'];
				if (OUT.good())
				{
					OUT<<ID<<"\t"<<Length<<"\t"<<N_Aleng<<"\t"<<eff_Acout<<"\t"<<GC_Aleng<<"\t"<<samllleter<<endl;
				}
				OUT2<<ID<<"\t"<<Length<<"\t"<<N_Aleng<<"\t"<<eff_Acout<<"\t"<<GC_Aleng<<"\t"<<samllleter<<endl;
			}
			OUT.close();
			OUT2.close();

		}
	}

	INID.close();
	delete para_A1 ;
	return 0;

}
#endif
///////// swimming in the sky and flying in the sea ////////////

