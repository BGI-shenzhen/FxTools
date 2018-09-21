
#ifndef FQ_ChangQ_H_
#define FQ_ChangQ_H_ 


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <cstdlib>

#include "../ALL/gzstream.h"
#include "../ALL/comm.h"
#include "../ALL/DataClass.h"

using namespace std;

int  print_usageChangQ ()
{
	cout <<""
		"\n"
		"Usage: changeQ  -i <In.fq>  -o <out.fq>\n"
		"\n"
		"\t\t-i     <str>   input FASTQ\n"
		"\t\t-o     <str>   output file\n"
		"\n"
		"\t\t-s   <int>     1:Sanger2Solexa 2: Solexa2Sanger (-31) [1]\n"
		"\t\t               3:ASCII33-->ASCII64[+31] with ResetID\n"
		"\t\t               4:ASCII33-->ASCII64[+31] with ResetID & MaxQ:h\n"
		"\t\t-h             show more details for help\n" 
		"\n";
	return 1;
}
void More_HelpFQ14()
{
	cout<<"\n\
\n\
\t\t1.  changeQ  -i <in.fq> -o AAA -s 1\n\
\t\tthis will convert the quality score of the input FASTQ from Sanger to Solexa and output to a compressed file named AAA in current directory.\n\
\n\
\n";
}

//
int parse_cmdChangQ(int argc, char **argv , In3str1v  * para )
{
	if (argc <=1  ) {print_usageChangQ();return 0;}

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

		if (flag  == "InFq" || flag  == "i")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->InStr1=argv[i];
		}
		else if (flag  ==  "OutFq" || flag  == "o")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->InStr2=argv[i];
		}
		else if (flag  ==  "Model"  ||  flag  ==  "s")
		{
			if(i + 1 == argc) {  LogLackArg( flag ) ; return 0;}
			i++;
			int A=atoi(argv[i]);
			if  ( A<0 ||  A>4 )
			{
				cerr<<"Model  Wrong"<<endl; return 0;
			}
			para->InInt=A;
		}
		else if (flag  == "help"   ||  flag  ==  "h")
		{
			More_HelpFQ14();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0 ;
		}
	}
	if  ((para->InStr2).empty() || (para->InStr1).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	(para->InStr2)=add_Asuffix(para->InStr2);
	return 1 ;
}

////////////////////////////

inline void trimA (  string & Quli )
{
	string::size_type i_end=Quli.size();
	for(string::size_type ix=0 ; ix<i_end ; ix++)
	{
		Quli[ix]=Quli[ix]+31;
	}

}

inline void trimB (  string & Quli )
{
	string::size_type i_end=Quli.size();
	for(string::size_type ix=0 ; ix<i_end ; ix++)
	{
		Quli[ix]=Quli[ix]-31;
	}
}

inline void ID4_too_2 (  string & ID )
{
	vector<string> Temp;
	vector<string> inf;
	split(ID,Temp," \t");
	split(Temp[1],inf,":");
	ID=Temp[0]+"#"+inf[inf.size()-1]+"/"+inf[0];
}

inline void trimD (  string & Quli )
{
	string::size_type i_end=Quli.size();
	for(string::size_type ix=0 ; ix<i_end ; ix++)
	{
		Quli[ix]=Quli[ix]+31;
		if ( Quli[ix] > 'h')
		{
			Quli[ix]=  'h' ;
		}
	}
}


//programme entry
///////// swimming in the sky and flying in the sea ////////////
//int main(int argc, char **argv)
int FQ_ChangQ_main(int argc, char **argv)
{
	In3str1v * para = new In3str1v;
	para->InInt=1;
	if (parse_cmdChangQ(argc, argv ,para )==0)
	{
		delete  para ; 
		return 1;
	}
	igzstream IN_1 ((para->InStr1).c_str(),ifstream::in); // ifstream  + gz 
	ogzstream OUT_1 ((para->InStr2).c_str());

	if(!IN_1.good())
	{
		cerr << "open IN File error: "<<para->InStr1<<endl;
		return 1;
	}
	if(!OUT_1.good())
	{
		cerr << "open OUT File error: "<<para->InStr2<<endl;
		return 1;
	} 
	if ( ( para->InInt) ==1 )
	{
		while(!IN_1.eof())
		{
			string ID_1 ,seq_1,temp_1,Quly_1 ;
			getline(IN_1,ID_1);
			getline(IN_1,seq_1);
			getline(IN_1,temp_1);
			getline(IN_1,Quly_1);
			trimA (Quly_1 );
			if (ID_1.length()<=0)  { continue  ; }
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";
		}
	}
	else if  ( ( para->InInt) ==3 )
	{
		while(!IN_1.eof())
		{
			string ID_1 ,seq_1,temp_1,Quly_1 ;
			getline(IN_1,ID_1);
			getline(IN_1,seq_1);
			getline(IN_1,temp_1);
			getline(IN_1,Quly_1);
			trimA (Quly_1 );
			if (ID_1.length()<=0)  { continue  ; }
			ID4_too_2(ID_1);
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";
		}
	}
	else if    ( ( para->InInt) ==2 )
	{
		while(!IN_1.eof())
		{
			string ID_1 ,seq_1,temp_1,Quly_1 ;
			getline(IN_1,ID_1);
			getline(IN_1,seq_1);
			getline(IN_1,temp_1);
			getline(IN_1,Quly_1);
			trimB (Quly_1 );
			if (ID_1.length()<=0)  { continue  ; }
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";
		}
	}
	else if    ( ( para->InInt) ==4 )
	{
		while(!IN_1.eof())
		{
			string ID_1 ,seq_1,temp_1,Quly_1 ;
			getline(IN_1,ID_1);
			getline(IN_1,seq_1);
			getline(IN_1,temp_1);
			getline(IN_1,Quly_1);
			trimD (Quly_1);
			if (ID_1.length()<=0)  { continue  ; }
			ID4_too_2(ID_1);
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";
		}
	}

	IN_1.close();
	OUT_1.close();
	return 0 ;
}

#endif  // FQ_Filter_H_

///////// swimming in the sky and flying in the sea ////////////
