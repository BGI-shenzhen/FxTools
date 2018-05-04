#ifndef RefRand_H_
#define RefRand_H_

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


using namespace std;
typedef long long  llong ;

int  print_Ausage_rand()
{
	cout <<""
		"\n"
		"\tUsage: rand  -i <in.fa> \n"
		"\n"
		"\t\t-i    <str>    Input fa for stat the length\n"
		"\n"
		"\t\t-o    <str>    Output file,otherwise[STDOUT]\n"
		"\n"
		"\t\t-p    <float>  read rand out proportion [0.1]\n"
		"\t\t-s    <int>    rand seed,default time\n"
		"\n"
		"\t\t-h             show this help infomation\n"
		"\n";
	return 1;
}


int parse_Acmd_Rand(int argc, char **argv , ParaClass * para_A1 )
{
	if (argc <=1 ) {print_Ausage_rand();return 0;}

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
		else if (flag  ==  "OutPut"  ||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_A1->InStr2=argv[i];
		}
		else if (flag  == "proportion" ||  flag  == "p")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_A1->Infloat1=atof(argv[i]);
		}
		else if (flag  ==  "seed"  ||  flag  == "s")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para_A1->InInt1=atoi(argv[i]);
		}
		else if (flag  == "help"  ||  flag  == "h")
		{
			print_Ausage_rand();return 0;
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




//int main(int argc, char *argv[])
int FA_Rand_main(int argc, char *argv[])
{
	ParaClass * para_A1 = new ParaClass;
	para_A1->Infloat1=0.1;
	if(parse_Acmd_Rand(argc, argv, para_A1 )==0)
	{
		delete  para_A1 ;
		return 1;
	}

	igzstream  INID (( para_A1->InStr1).c_str(),ifstream::in);
	if(!INID.good())
	{
		cerr << "open InputFile error: "<<(para_A1->InStr1)<<endl;
		delete para_A1 ; return 1;
	}

	string seqAA;
	getline(INID, seqAA ,'>');

	if (para_A1->InInt1==0)
	{
		srand((unsigned)time(NULL));
	}
	else
	{
		srand(para_A1->InInt1);
	}

	int Y=int(para_A1->Infloat1*100);
	int X=0;
	
	if (!(para_A1->InStr2).empty())
	{
		para_A1->InStr2=add_Asuffix(para_A1->InStr2);
		gzFile OUTGZ;
		OUTGZ = gzopen ((para_A1->InStr2).c_str(), "wb");
		char EE[2] ; EE[0]='>';
		char FF[2];  FF[0]='\n';
		while(!INID.eof())
		{
			string line;
			getline(INID, line,'\n');
			if (line.empty())
			{
				continue;			
			}
			X=(rand()%100);
			getline(INID, seqAA ,'>');
			if (X>Y)  { continue; }
//			gzprintf(OUTGZ,">%s\n%s",line.c_str(),seqAA.c_str());
			gzwrite (OUTGZ, EE,1);
			gzwrite (OUTGZ,line.c_str(),line.length());
			gzwrite (OUTGZ, FF,1);
			gzwrite (OUTGZ,seqAA.c_str(),seqAA.length());
		}
		gzclose(OUTGZ);
	}
	else
	{
		while(!INID.eof())
		{
			string line;
			getline(INID, line,'\n');
			if (line.empty())
			{
				continue;
			}
			X=(rand()%100);
			getline(INID, seqAA ,'>');
			if (X>Y)  { continue  ; }
			fprintf(stdout,">%s\n%s",line.c_str(),seqAA.c_str());
		}
	}


	INID.close();
	delete para_A1 ;
	return 0;

}
#endif
///////// swimming in the sky and flying in the sea ////////////

