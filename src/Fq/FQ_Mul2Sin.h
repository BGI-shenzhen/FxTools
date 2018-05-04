#ifndef FQ_Mul2Sin_H_
#define FQ_Mul2Sin_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include "../ALL/comm.h"
#include "../ALL/gzstream.C"
#include "../ALL/kseq.h"
#include "../ALL/DataClass.h"
#include <gzstream.h>

using namespace std;
typedef long long llong ;



int  print_U_Mul2Sin()
{
	cout <<""
		"\n"
		"\tUsage: valid -i <in.fq>  -o <out.fq>\n"
		"\n"
		"\t\t-i    <str>     Input muti line FASTQ file\n"
		"\t\t-o    <str>     Output single line FASTQ file\n"
		"\n"
		"\t\t-h              show this help\n" 
		"\n";
	return 1;
}


int parse_Mul2Sin_A30(int argc, char **argv,In3str1v *  para_A30 )
{
	if (argc <=2 ) {print_U_Mul2Sin();return 0 ;}
	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-' )
		{
			cerr << "command option error! please check." << endl;
			return 0 ;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InFq" ||  flag  == "i" )
		{
			if(i + 1 == argc) { LogLackArg(flag); return 0 ;}
			i++;
			para_A30->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut" || flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0 ;}
			i++;
			para_A30->InStr3=argv[i];
		}
		else if (flag  == "help" || flag  == "h")
		{
			print_U_Mul2Sin();return 0 ;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0 ;
		}
	}
	if  (para_A30->InStr1.empty()  ||  para_A30->InStr3.empty())
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0 ;
	}
	para_A30->InStr3=add_Asuffix(para_A30->InStr3);
	return 1 ;
}


int FQ_Mul2Sin_main( int argc,char *argv[] )
//int main (int argc,char *argv[])
{
	In3str1v *  para_A30 =new In3str1v ;
	if( parse_Mul2Sin_A30( argc, argv,para_A30) ==0 )
	{
		delete  para_A30 ;
		return 1 ;
	}
	
	igzstream IN_A ((para_A30->InStr1).c_str(),ifstream::in); // ifstream  + gz
	ogzstream OUT ((para_A30->InStr3).c_str());


	string ID_1,ID_2;
	getline(IN_A,ID_2);
	while(!IN_A.eof())
	{
		string seq_1,temp_1,Quly_1;
		ID_1=ID_2;
		getline(IN_A,temp_1);
		while((temp_1[0]!='+') &&  (!IN_A.eof()))
		{
			seq_1=seq_1+temp_1;
			getline(IN_A,temp_1);
		}
		getline(IN_A,Quly_1);
		getline(IN_A,ID_2);
		string::size_type SeqLen=seq_1.length();
		
		while( ( (ID_2[0]!='@')  ||  (Quly_1.length()<SeqLen) )  && (!IN_A.eof())  )
		{
			Quly_1=Quly_1+ID_2;
			getline(IN_A,ID_2);
		}		

		OUT<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";
	}

	if (ID_1!=ID_2 && (!ID_2.empty()))
	{
		cerr<<"Sequence(s) with Error at Empty seq :"<<ID_2<<endl;
	}
	IN_A.close();
	OUT.close();
	delete para_A30 ;
	return 0 ;
}

///////// swimming in the sky and flying in the sea ////////////
#endif 
