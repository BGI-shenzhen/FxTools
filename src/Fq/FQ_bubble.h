
#ifndef FQ_Bubble_H_
#define FQ_Bubble_H_ 

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

using namespace std;

int  print_usage_C01 ()
{
	cout <<""
		"\n"
		"Usage:bubble -i A1.fq  A2.fq -o  B1.fq  B2.fq -b  8,9 \n"
		"\n"
		"\t\t-i    <str>   File name of InFq1/InFq2 Input\n"
		"\t\t-o    <str>   File name of OutFq1/OutFq2.gz output\n"
		"\t\t-a    <str>   The Site of Bubble along the FqRead1 such as [3,8]\n"
		"\t\t-b    <str>   The Site of Bubble along the FqRead2 such as [1]\n"
		"\t\t-h            show this help\n" 
		"\n";
	return 1;
}
//
int parse_cmd_C01(int argc, char **argv , Para_A24 * para)
{
	if (argc <=4  ) {print_usage_C01();return 0;}

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

		if (flag  == "InFq1" || flag  == "i")
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
		else if (flag  ==  "InFq2" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->InFq2=argv[i];
		}
		else if (flag  ==  "OutFq1" || flag  == "o")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->OutFq1=argv[i];
			if ( (i+1) < argc  && (argv[i+1][0]!='-'))
			{
				i++;
				para->OutFq2=argv[i];
			}
		}
		else if (flag  ==  "OutFq2" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ;return 0;}
			i++;
			para->OutFq2=argv[i];
		}
		else if (flag  ==  "BubbleSite1" ||  flag  == "a")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->adapter1=argv[i];
		}
		else if (flag  ==  "BubbleSite2"  ||  flag  == "b" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ;return 0;}
			i++;
			para->adapter2=argv[i];
		}
		else if (flag  == "help" || flag  == "h")
		{
			print_usage_C01();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if ((para->OutFq1).empty() || (para->InFq1).empty() || (para->OutFq2).empty() || (para->InFq2).empty())
	{
		cerr<< "InPut OutPut Lack argument for the must"<<endl;
		return 0;
	}
	if (  (para->adapter1)== "NA " &&  (para->adapter2) == "NA" )
	{
		cerr<< "BubbleSite Lack argument for the must"<<endl;
		return 0;
	}
	(para->OutFq1)=add_Asuffix(para->OutFq1);
	(para->OutFq2)=add_Asuffix(para->OutFq2);
	return 1 ;
}

////////////////////////////

///////// swimming in the sky and flying in the sea ////////////
int FQ_bubble_main(int argc, char **argv)
	//int main(int argc, char **argv)
{
	Para_A24 * para = new Para_A24 ;        

	int Flag_para=parse_cmd_C01(argc, argv ,para ) ;
	if ( Flag_para ==0)
	{
		delete  para ; 
		return 1;
	}

	igzstream IN_1 ((para->InFq1).c_str(),ifstream::in); // ifstream  + gz 
	igzstream IN_2 ((para->InFq2).c_str(),ifstream::in); // ifstream  + gz 
	ogzstream OUT_1 ((para->OutFq1).c_str());
	ogzstream OUT_2 ((para->OutFq2).c_str());

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

	long long read_Fisrt=0 , read_Second=0 ,Base_Fisrt=0 ,Base_Second =0;


	if ( (para->adapter1)!="NA"   &&  (para->adapter2)=="NA" )
	{

		vector <string> inf ;
		vector <int> ID1 ;
		split( (para->adapter1) ,inf,",") ;
		for (unsigned int ii=0 ; ii<inf.size(); ii++)
		{
			ID1.push_back(atoi(inf[ii].c_str())-1);
		}

		int lengthID1=ID1.size();
		inf.clear();

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
			string::size_type i_end=seq_1.size()-1;
			string::size_type s_start=0 ,s_end =i_end ;
			int IF_RM=0;
			for (unsigned int ii=0 ; ii<lengthID1 ; ii++)
			{
				if ( seq_1[ID1[ii]]=='N'  )
				{
					IF_RM=1;
					break  ;
				}
			}

			if (IF_RM == 0 )
			{
				read_Second++;
				OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";
				OUT_2<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n";
			}
		}
	}
	else if ( (para->adapter1)=="NA"   &&  (para->adapter2)!="NA" )
	{

		vector <string> inf ;
		vector <int> ID2 ;
		split( (para->adapter2) ,inf,",") ;

		for (unsigned int ii=0 ; ii<inf.size(); ii++)
		{
			ID2.push_back(atoi(inf[ii].c_str())-1);
		}

		int lengthID2=ID2.size();
		inf.clear();

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
			string::size_type i_end=seq_1.size()-1;
			string::size_type s_start=0 ,s_end =i_end ;
			int IF_RM=0;
			for (unsigned int ii=0 ; ii<lengthID2 ; ii++)
			{
				if ( seq_2[ID2[ii]]=='N'  )
				{
					IF_RM=1;
					break  ;
				}
			}

			if (IF_RM == 0 )
			{
				read_Second++;
				OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";
				OUT_2<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n";
			}
		}
	}
	else
	{

		vector <string> inf ;
		vector <int> ID1 ;
		split( (para->adapter1) ,inf,",") ;
		for (unsigned int ii=0 ; ii<inf.size(); ii++)
		{
			ID1.push_back(atoi(inf[ii].c_str())-1);
		}

		int lengthID1=ID1.size();
		inf.clear();

		vector <int> ID2 ;
		split( (para->adapter2) ,inf,",") ;

		for (unsigned int ii=0 ; ii<inf.size(); ii++)
		{
			ID2.push_back(atoi(inf[ii].c_str())-1);
		}

		int lengthID2=ID2.size();
		inf.clear();

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
			string::size_type i_end=seq_1.size()-1;
			string::size_type s_start=0 ,s_end =i_end ;
			int IF_RM=0;
			for (unsigned int ii=0 ; ii<lengthID1 ; ii++)
			{
				if ( seq_1[ID1[ii]]=='N'  )
				{
					IF_RM=1;
					break  ;
				}
			}

			for (unsigned int ii=0 ; ii<lengthID2 ; ii++)
			{
				if ( seq_2[ID2[ii]]=='N'  )
				{
					IF_RM=1;
					break  ;
				}
			}


			if (IF_RM == 0 )
			{
				read_Second++;
				OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";
				OUT_2<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n";
			}
		}

	}



	IN_2.close();
	OUT_2.close();
	IN_1.close();
	OUT_1.close();

	read_Fisrt=(read_Fisrt)*2 ;
	read_Second=(read_Second)*2;
	cout<<"##Final Stat OUT##"<<endl;
	cout<<"#Original Read\t"<<read_Fisrt<<"\t#After Filter Read\t"<<read_Second<<endl;

	delete para ;
	return 0 ;
}

#endif  // FQ_Bubble_H_


///////// swimming in the sky and flying in the sea ////////////
