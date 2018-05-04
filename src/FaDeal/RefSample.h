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

int  print_AusageFA02()
{
	cout <<""
		"\n"
		"\tUsage: getSN  -i <in.fa>  -s 5-10\n"
		"\tUsage: getSN  -i <in.fa>  -s 5,3,8,10\n"
		"\n"
		"\t\t-i      <str>   InPut fa for geting  Sub range Seq\n"
		"\t\t-s      <str>   Get Seq by specified order range\n"
		"\t\t                Only give one num will only extract one Seq\n"
		"\n"
		"\t\t-o      <str>   OutPut Fasta file or [STDOUT]\n" 
		"\n"
		"\t\t-h              show this help\n" 
		"\n";
	return 1;
}


int parse_AcmdFA02(int argc, char **argv , In3str1v * paraFA02 )
{
	if (argc <=2 ) {print_AusageFA02();return 0;}

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
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA02->InStr1=argv[i];
		}
		else if (flag  ==  "Sample" ||  flag  == "s")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA02->InStr3=argv[i];
		}
		else if (flag  ==  "OutPut"||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA02->InStr2=argv[i];
		}
		else if (flag  == "help" ||  flag  == "h" )
		{
			print_AusageFA02();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((paraFA02->InStr1).empty() || (paraFA02->InStr3).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	if (!(paraFA02->InStr2).empty() )
	{
		(paraFA02->InStr2)=add_Asuffix(paraFA02->InStr2) ;
	}

	return 1 ;
}

int FA_Order_main(int argc, char *argv[])
//int main(int argc, char *argv[])
{
	In3str1v * paraFA02 = new In3str1v;
	if( parse_AcmdFA02(argc, argv, paraFA02 )==0)
	{
		delete  paraFA02 ;
		return 1 ;
	}

	int  Start=0; int End=0;
	map <int , bool > IDHash ;


	if  ((paraFA02->InStr3).find(',')==string::npos)
	{
		vector<string> Temp;
		split((paraFA02->InStr3),Temp,"-");
		Start=atoi(Temp[0].c_str());
		End=atoi(Temp[Temp.size()-1].c_str());
		if (Start>End)
		{
			cerr<<Start<<"biger than"<<endl;
			delete  paraFA02 ; return 1;
		}

		for (int jj=Start ; jj<=End  ;jj++)
		{		 
			IDHash.insert(map <int,bool>:: value_type(jj,true));
		}
	}
	else
	{
		vector<string> Temp;
		split((paraFA02->InStr3),Temp,",");
		int EE=Temp.size() ;
		for (int jj=0; jj<EE; jj++)
		{
			int BB=atoi(Temp[jj].c_str());
			IDHash.insert(map <int,bool>:: value_type(BB,true));
		}
		map <int, bool> :: iterator itID ;
		itID=IDHash.begin() ;
		Start=itID->first;
		End=itID->first;   itID++;
		for(  ; itID!=IDHash.end(); itID++)
		{
			if (Start>(itID->first))
			{
				Start=(itID->first);
			}
			if (End<(itID->first))
			{
				End=(itID->first);
			}
		}
	}


	igzstream IN ((paraFA02->InStr1).c_str(),ifstream::in);
	if(!IN.good())
	{
		cerr << "open InputFile error: "<<(paraFA02->InStr1)<<endl;
		delete  paraFA02 ; return 1;
	}
	if ((paraFA02->InStr2).empty())
	{
		string seqAA ;
		getline(IN, seqAA, '>');
		int A=0;
		while(!IN.eof())
		{
			string chr_line ;
			A++;
			getline(IN, chr_line , '\n') ;
			getline(IN, seqAA , '>') ;
			if (A>=Start)
			{
				if  (IDHash.find(A)==IDHash.end())
				{
					continue ;
				}
				cout<<">"<<chr_line<<"\n"<<seqAA;
			}
			if (A>=End)
			{
				break ;
			}
		}
	}
	else
	{

		ogzstream  OUT  ((paraFA02->InStr2).c_str());
		if (OUT.fail())
		{
			cerr << "open OUTFile error: "<<(paraFA02->InStr2)<<endl;
			delete  paraFA02 ; return 1;
		}

		string seqAA ;
		getline(IN, seqAA, '>');
		int A=0;
		while(!IN.eof())
		{
			string chr_line ;
			A++;
			getline(IN, chr_line , '\n') ;
			getline(IN, seqAA , '>') ;
			if (A>=Start)
			{
				if  (IDHash.find(A)==IDHash.end())
				{
					continue ;
				}
				OUT<<">"<<chr_line<<"\n"<<seqAA;
			}
			if (A>=End)
			{
				break ;
			}
		}

		OUT.close();
	}
	delete paraFA02 ;
	return 0;

}

///////// swimming in the sky and flying in the sea ////////////

