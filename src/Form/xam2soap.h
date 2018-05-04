#ifndef xam2soap_H_
#define xam2soap_H_

#include <iostream>
#include <fstream>
#include <gzstream.h>
#include "./xam2soap/TSamCtrl.h"
#include "./xam2soap/Ttools.h"
#include "../ALL/comm.h"
#include "../ALL/DataClass.h"
using namespace std;


int  print_Ausage_A14()
{
	cout <<""
		"\n"
		"\tUsage: Xam2Soap  -b <in.bam>  -o <Out.soap>\n"
		"\n"
		"\t\t-s     <str>   Input Sam Format file\n"
		"\t\t-b     <str>   Input Bam Format file\n"
		"\t\t-o     <str>   Output the SoapFormat file\n"
		"\n"
		"\t\t-Q     <int>   SeqQ shift Trans:(-31/+31/0)[0]\n"
		"\t\t               Sanger<->Illumina:! <-> @ \n"
		"\t\t-h             show this help\n" 
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

		if (flag  == "InSam"  ||   flag  == "s")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_A14->InStr1=argv[i];
		}
		else if (flag  == "InBam"  ||   flag  == "b")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_A14->InStr2=argv[i];
			para_A14->TF=false;
		}
		else if (flag  ==  "OutPut" ||   flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_A14->InStr3=argv[i];
		}
		else if (flag  ==  "QShift" ||   flag  == "Q" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_A14->InInt=atoi(argv[i]);
		}

		else if (flag  == "help" ||   flag  == "h")
		{
			print_Ausage_A14();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if ((para_A14->InStr3).empty()||((para_A14->InStr1).empty()&&(para_A14->InStr2).empty()))
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	(para_A14->InStr3)= add_Asuffix ( (para_A14->InStr3) ) ;

	return 1 ;
}

int Xam2Soap_main(int argc, char **argv)
	//int main(int argc, char **argv)
{

	In3str1v * para_A14 = new In3str1v;
	if (parse_Acmd_A14(argc, argv,para_A14  )==0)
	{
		delete  para_A14 ;
		return 1;
	}
	ogzstream OUT(para_A14->InStr3.c_str()) ;
	if (!OUT.good())
	{
		cerr<<"Can't open OutFile "<<para_A14->InStr3<<endl;
		return 1;
	}

	char in_mode[5] ={ 0 };
	string sampath ;

	if ((para_A14->TF))
	{
		in_mode[0]='r';
		sampath=(para_A14->InStr1);        
	}
	else
	{
		in_mode[0]='r';in_mode[1]='b';
		sampath=(para_A14->InStr2);
	}

	TSamCtrl samhandle;

	samhandle.open(sampath.c_str(),in_mode);
	string line;
	if ((para_A14->InInt)==0)
	{
		while(samhandle.readline(line)!=-1)
		{
			line = Talignment_format(line);
			if(line.compare(NOUSE_ALIGNMENT)==0) continue;
			OUT<<line<<endl;
		}
	}
	else
	{
		while(samhandle.readline(line)!=-1)
		{
			line = Talignment_format(line);
			if(line.compare(NOUSE_ALIGNMENT)==0) continue;
			vector<string> inf;
			split(line,inf," \t");
			size_t length=inf[2].length();
			for (size_t ii=0 ; ii<length ; ii++)
			{
				(inf[2])[ii]=((inf[2])[ii])+(para_A14->InInt) ;
			}
			length=inf.size();
			OUT<<inf[0];
			for (size_t ii=1;ii<length;ii++)
			{
				OUT<<"\t"<<inf[ii];
			}
			OUT<<endl;
		}
	}

	OUT.close();
	delete para_A14 ;
	return 0;
}
#endif // xam2soap_H_
///////// swimming in the sky and flying in the sea ////////////
