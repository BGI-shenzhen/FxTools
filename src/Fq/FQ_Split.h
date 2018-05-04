
#ifndef FQ_Split_H_
#define FQ_Split_H_ 

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <gzstream.h>
#include "../ALL/comm.h"
#include "FQ_Filter.h"
#include <sys/stat.h>

using namespace std;

int  print_usage_Sp03 ()
{
	cout <<""
		"\n"
		"Usage:splitFq -i A_1.fq A_2.fq \n"
		"\n"
		"\t\t-i    <str>   File name of InFq1/InFq2 Input\n"
		"\n"
		"\t\t-o    <str>   Output Dir [pwd]\n"
		"\t\t-n    <int>   Max Fq Read Number for per Fq[10000000]\n"
		"\t\t-h            show this help\n" 
		"\n";
	return 1;
}

int parse_cmd_Sp03(int argc, char **argv , Para_A24 * para)
{
	if (argc <=4  ) {print_usage_Sp03();return 0;}

	int err_flag = 0;
	for(int i = 1; i < argc || err_flag; i++)
	{
		if(argv[i][0] != '-' )
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i];
		flag=replace_all(flag,"-","");

		if (flag  == "InFq1" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->InFq1=argv[i];
		}
		else if (flag  ==  "InFq2" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->InFq2=argv[i];
		}
		else if (flag  ==  "i" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->InFq1=argv[i];
			if ( (i+1) < argc  && (argv[i+1][0]!='-'))
			{
				i++;
				para->InFq2=argv[i];
			}
		}
		else if (flag  ==  "MaxNum" ||  flag  ==  "n")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->adapter1=argv[i];
		}
		else if (flag  ==  "OutDir"  || flag  == "o")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ;return 0;}
			i++;
			para->adapter2=argv[i];
		}
		else if (flag  == "help"  || flag  == "h" )
		{
			print_usage_Sp03();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if ( (para->InFq1).empty() )
	{
		cerr<< "InPut OutPut Lack argument for the must"<<endl;
		return 0;
	}
	return 1 ;
}

///////// swimming in the sky and flying in the sea ////////////
int FQ_Split_main(int argc, char **argv)
//int main(int argc, char **argv)
{
	Para_A24 * para = new Para_A24 ;        
	(para->adapter1)="10000000";
	(para->adapter2)="./";
	int Flag_para=parse_cmd_Sp03(argc, argv ,para ) ;
	if ( Flag_para ==0)
	{
		delete  para ; 
		return 1;
	}

	long long read_Fisrt=-1 ;
	int Count=1;
	int MaxNum= atoi((para->adapter1).c_str());
	mkdir(((para->adapter2)).c_str() , 0755 ) ;

	if (!(para->InFq2).empty())
	{
		igzstream IN_1 ((para->InFq1).c_str(),ifstream::in); // ifstream  + gz
		igzstream IN_2 ((para->InFq2).c_str(),ifstream::in); // ifstream  + gz
		if(!IN_1.good())
		{
			cerr << "open IN File error: "<<para->InFq1<<endl;
			return 1;
		}
		if(!IN_2.good())
		{
			cerr << "open IN File error: "<<para->InFq2<<endl;
			return 1;
		}

		para->OutFq1= para->adapter2+"/OUT1_1.fq.gz";
		para->OutFq2= para->adapter2+"/OUT1_2.fq.gz";

		ogzstream OUT_1 ((para->OutFq1).c_str());
		ogzstream OUT_2 ((para->OutFq2).c_str());

		if(!OUT_1.good())
		{
			cerr << "open OUT File error: "<<para->OutFq1<<endl;
			return 1;
		} 
		if(!OUT_2.good())
		{
			cerr << "open OUT File error: "<<para->OutFq2<<endl;
			return 1;
		}

		string ID_1 ,seq_1,temp_1,Quly_1 ;
		string ID_2 ,seq_2,temp_2,Quly_2 ;


		while(!IN_1.eof())
		{
			getline(IN_1,ID_1);
			getline(IN_1,seq_1);
			getline(IN_1,temp_1);
			getline(IN_1,Quly_1);
			getline(IN_2,ID_2);
			getline(IN_2,seq_2);
			getline(IN_2,temp_2);
			getline(IN_2,Quly_2);
			if (ID_1.length()<=0)  { continue  ; }

			read_Fisrt++;
			if ((int(read_Fisrt/MaxNum)+1)<=Count)
			{
				OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";
				OUT_2<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n";
			}
			else
			{
				Count++;
				OUT_2.close();
				OUT_1.close();
				para->OutFq1= (para->adapter2)+"/OUT"+Int2Str(Count)+"_1.fq.gz";
				para->OutFq2= (para->adapter2)+"/OUT"+Int2Str(Count)+"_2.fq.gz";

				OUT_1.open((para->OutFq1).c_str());
				OUT_2.open((para->OutFq2).c_str());
				OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";
				OUT_2<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n";

			}
		}

		read_Fisrt++;
		IN_2.close();
		OUT_2.close();
		IN_1.close();
		OUT_1.close();
		cout<<"All Read Pair Number is "<<read_Fisrt<<"\n and total split to File Number is " <<Count<<" *2 "<<endl;

	}


	else
	{

		igzstream IN_1 ((para->InFq1).c_str(),ifstream::in); // ifstream  + gz 
		if(!IN_1.good())
		{
			cerr << "open IN File error: "<<para->InFq1<<endl;
			return 1;
		}

		para->OutFq1= para->adapter2+"/OUT1.fq.gz";
		ogzstream OUT_1 ((para->OutFq1).c_str());

		if(!OUT_1.good())
		{
			cerr << "open OUT File error: "<<para->OutFq1<<endl;
			return 1;
		} 

		string ID_1 ,seq_1,temp_1,Quly_1 ;


		while(!IN_1.eof())
		{
			getline(IN_1,ID_1);
			getline(IN_1,seq_1);
			getline(IN_1,temp_1);
			getline(IN_1,Quly_1);
			if (ID_1.length()<=0)  { continue  ; }

			read_Fisrt++;
			if ((int(read_Fisrt/MaxNum)+1)<=Count)
			{
				OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";
			}
			else
			{
				Count++;
				OUT_1.close();
				para->OutFq1= (para->adapter2)+"/OUT"+Int2Str(Count)+".fq.gz";

				OUT_1.open((para->OutFq1).c_str());
				OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";

			}
		}


		read_Fisrt++;

		IN_1.close();
		OUT_1.close();

		cout<<"All Read Pair Number is "<<read_Fisrt<<"\n and total split to File Number is " <<Count<<" *1 "<<endl;


	}


	delete para ;
	return 0;
}

#endif  

///////// swimming in the sky and flying in the sea ////////////
