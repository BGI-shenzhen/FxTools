
#ifndef FQ_cutIndex_H_
#define FQ_cutIndex_H_ 


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
#include "../ALL/DataClass.h"

using namespace std;

int  print_usage_B01 ()
{
	cout <<""
		"\n"
		"Usage:cut  -i In.fq  -o Out.fq\n"
		"\n"
		"\t\t-i    <str>   File name of InFq Input\n"
		"\n"
		"\t\t-o    <str>   OutFile Fq.gz,otherwise[STDOUT]\n"
		"\t\t-s    <int>   The Start position to cut X bp [5]\n"
		"\t\t-e    <int>   The End site of cut[Rleng]\n"
		"\t\t-h                show this help\n"
		"\n";
	return 1;
}
//
int parse_cmd_B01(int argc, char **argv , ParaClass * para )
{
	if (argc <=4  ) {print_usage_B01();return 0;}

	int err_flag = 0;
	for(int i = 1; i < argc || err_flag; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InFq"  || flag  == "i")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->InPut1=argv[i];
		}
		else if (flag  ==  "OutFq"  || flag  == "o")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->OutPut1=argv[i];
		}
		else if (flag  ==  "StartCut" || flag  == "s")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->InInt1=atoi(argv[i]);
		}
		else if (flag  ==  "EndCut" || flag  == "e")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->InInt2=atoi(argv[i]);
		}
		else if ( flag  == "h")
		{
			print_usage_B01();
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if ( (para->InPut1).empty() )
	{
		cerr<< "-InFq -OutFq lack argument for the must"<<endl;
		return 0;
	}
	para->InInt1--;
	return 1;
}




///programme entry
///////// swimming in the sky and flying in the sea ////////////
//int main(int argc, char **argv)
int FQ_IndexCut_main(int argc, char **argv)
{
	ParaClass * para = new ParaClass;
	para->InInt2=-1;
	para->InInt1=5;
	int Flag_para=parse_cmd_B01(argc, argv ,para ) ;
	if ( Flag_para ==0)
	{
		delete  para ; 
		return 1;
	}

	igzstream INFQ ((para->InPut1).c_str(),ifstream::in); // ifstream  + gz

	if(!INFQ.good())
	{
		cerr << "open IN File error: "<<para->InPut1<<endl;
		return 1;
	}


	int maxReadLeng=0 ;
	int BB=0;
	while( (!INFQ.eof()) && (BB< 16888 ))
	{
		string  line ;
		getline(INFQ,line);
		getline(INFQ,line);
		int tmp=line.length();
		if (tmp> maxReadLeng)
		{
			maxReadLeng=tmp;
		}
		getline(INFQ,line);
		getline(INFQ,line);
		BB++;
	}
	INFQ.close();

	int cutLength=maxReadLeng-(para->InInt1);
	if  ( (para->InInt1)  > maxReadLeng)
	{
		cerr << "Start site shoud smaller ReadLength"<<endl;
		delete  para ; 
		return 1;
	}
	else if ((para->InInt2)<1)
	{
		
	}
	else if ((para->InInt2) > maxReadLeng )
	{
		cerr<<"warning: cut End site : "<<(para->InInt2)<<"\tbiger than the Read length "<<maxReadLeng<<"\n";
	}
	else if ( (para->InInt2) < (para->InInt1))
	{
		cerr<<"End cut site fail ,End site should biger Start site"<<endl;
	}
	else
	{
		cutLength=(para->InInt2) -(para->InInt1);
	}


	igzstream IN_1 ((para->InPut1).c_str(),ifstream::in); // ifstream  + gz 
	if(!IN_1.good())
	{
		cerr << "open IN File error: "<<para->InPut1<<endl;
		return 1;
	}


	string ID_1 ,seq_1,temp_1,Quly_1 ;


	if ((para->OutPut1).empty())
	{
		while(!IN_1.eof())
		{
			getline(IN_1,ID_1);
			if (ID_1.empty())
			{
				continue ;
			}
			getline(IN_1,seq_1);
			getline(IN_1,temp_1);
			getline(IN_1,Quly_1);
			string A , B ;
			A = seq_1.substr(para->InInt1,cutLength); 
			B = Quly_1.substr(para->InInt1,cutLength);
			cout<<ID_1<<"\n"<<A<<"\n"<<temp_1<<"\n"<<B<<"\n";
		}
	}
	else
	{
		(para->OutPut1)=add_Asuffix(para->OutPut1);
		ogzstream OUT_1 ((para->OutPut1).c_str());
		if(!OUT_1.good())
		{
			cerr << "open OUT File error: "<<para->OutPut1<<endl;
			return 1;
		} 


		while(!IN_1.eof())
		{
			getline(IN_1,ID_1);
			if (ID_1.empty())
			{
				continue ;
			}
			getline(IN_1,seq_1);
			getline(IN_1,temp_1);
			getline(IN_1,Quly_1);
			string A , B ;
			A = seq_1.substr(para->InInt1,cutLength); 
			B = Quly_1.substr(para->InInt1,cutLength);
			OUT_1<<ID_1<<"\n"<<A<<"\n"<<temp_1<<"\n"<<B<<"\n";
		}
		OUT_1.close();

	}	
	delete para ;
	return 0 ;
}

#endif  // FQ_cutIndex_H_


///////// swimming in the sky and flying in the sea ////////////
