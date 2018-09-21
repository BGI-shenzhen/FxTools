#ifndef FQ_RmDup_H_
#define FQ_RmDup_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <semaphore.h> 
#include <string>
#include <iomanip>
#include <math.h>
#include <zlib.h>
#include "FQ_Filter.h"

#define MAX_L 512
#define MAX_Q 128

using namespace std;
//  #define _A 0 #define _C 1    #define _G 2   #define _T 3   #define _N 4

int  print_usage_rmDup()
{
	cout <<""
		"\n"
		"Usage:rmdup  -i <1.fq> <2.fq> -o <out.1> <out.2>\n"
		"\n"
		"\t\t-i       <str>   input PE FASTQ. For SE FASTQ, just omit <2.fq>\n"
		"\t\t-o       <str>   output for PE FASTQ. For SE FASTQ, just omit <out.2>, default [STDOUT]\n"
		"\n"
		"\t\t-m       <int>   remove duplicates block by block. Regard X reads as a block to limit memory usage, default [10M]\n"
		"\n"
		"\t\t-h               show this help\n" 
		"\n";
	return 1;
}

void MoreHelpRd()
{
	cout<<"\n\
\t\t1.	rmdup -i <in.fq> <2.fq> -o <out.1> <out.2>\n\
\t\t  	rmdup -I <in.fq>  -o <out.fq>\n\
\t\t  this will check the first 10M reads for duplicates and then the next 10M, and so on, which means rmdup is not able to find the duplicate reads not in the same block. User could develop the result by either repeating this operation for multiple times or expanding the value of block (default is 10M).\n\
\n\
\n";
}

int parse_cmd_rmDup( int argc, char **argv, ParaClass * PaAA )
{
	if (argc <=1) {print_usage_rmDup();return 0;}

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

		if (flag  == "InFq1" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			PaAA->InPut1=argv[i];
		}
		else if (flag  ==  "i" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			PaAA->InPut1=argv[i];
				if ( (i+1) < argc  && (argv[i+1][0]!='-'))
			{
				i++;
				PaAA->InPut2=argv[i];
			}
		}
		else if (flag  ==  "o" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			PaAA->OutPut1=argv[i];
				if ( (i+1) < argc  && (argv[i+1][0]!='-'))
			{
				i++;
				PaAA->OutPut2=argv[i];
			}
		}

		else if (flag  ==  "InFq2" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			PaAA->InPut2=argv[i];
		}
		else if (flag  ==  "OutFq1" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			PaAA->OutPut1=argv[i];
		}

		else if (flag  == "OutFq2")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			PaAA->OutPut2=argv[i];
		}
		else if (flag  ==  "ReadBk" ||  flag=="m")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			string AA=argv[i];
			char A=argv[i][strlen(argv[i])-1];
			int II=1;
			if (A == 'k' ||  A == 'K') 
			{
			II=1024;
			}
			else if (A == 'M' ||  A == 'm')
			{
			II=1024*1024;
			}
			else if  (A == 'G' ||  A == 'g')
			{
			II=1024*1024*1024;
			}
			if  (II==1)
			{
			PaAA->InInt1=atoi(argv[i]);
			}
			else
			{
			  argv[i][strlen(argv[i])-1]='\0';
			  PaAA->InInt1=atoi(argv[i])*II;
			}
		}
		else if (flag  == "help" ||  flag=="h")
		{
			MoreHelpRd();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if ( (PaAA->InPut1).empty() )
	{
		cerr<< "-InFq1  lack argument for the must"<<endl;
		return 0;
	}

	if (!(PaAA->InPut2).empty())
	{
		if ( (PaAA->OutPut1).empty()  ||  (PaAA->OutPut2).empty()  )
		{
			cerr<< "-OutFq1   -OutFq2  must be setted with -InFq2  together"<<endl;
		}
		return 0;
	}

	if (!(PaAA->OutPut1).empty())
	{
		(PaAA->OutPut1)=add_Asuffix(PaAA->OutPut1);
	}

	if (!(PaAA->OutPut2).empty())
	{
		(PaAA->OutPut2)=add_Asuffix(PaAA->OutPut2);
	}

	return 1 ;
}


int FqRmDupSE ( ParaClass * PaAA  )
{
	igzstream IN_A ((PaAA->InPut1).c_str(),ifstream::in); // ifstream  + gz
	if(!IN_A.good())
	{
		cerr << "open IN File error: "<<(PaAA->InPut1)<<endl;
		return 1;
	}

	int TotalRead=0;
	int DupRead=0;

	int Tarry=(PaAA->InInt1)*4;
	string *FqData =new string[Tarry];
	cerr<<"Regard "<<(PaAA->InInt1)<<" read as a Blocks to remove Reads "<<endl;

	if ((PaAA->OutPut1).empty())
	{
		while(!IN_A.eof())
		{
			int ii =0 ;
			map <string ,int> SeqIndex ;
			map <string ,int> ::iterator it ;

			for ( ; ii< Tarry &&  (!IN_A.eof() ); ii+=4)
			{
				getline(IN_A,FqData[ii]);
				getline(IN_A,FqData[ii+1]);
				getline(IN_A,FqData[ii+2]);
				getline(IN_A,FqData[ii+3]);
				if  (FqData[ii+1].length()==0)	{		 continue ;		}
				it=SeqIndex.find(FqData[ii+1]);
				if (it!=SeqIndex.end())
				{
					DupRead++;
					it->second=ii;
				}
				else
				{
					SeqIndex.insert(map <string ,int > ::value_type(FqData[ii+1],ii));
				}
				TotalRead++;
			}
			for (it=SeqIndex.begin(); it!=SeqIndex.end();)
			{
				int II=it->second;
				fprintf(stdout,"%s\n%s\n%s\n%s\n",FqData[II].c_str(),FqData[II+1].c_str(),FqData[II+2].c_str(),FqData[II+3].c_str());
				SeqIndex.erase(it++);
			}

		}
		cerr<<"Total Read Number: "<<TotalRead<<"\n"<<"Remove Read Number "<<DupRead<<endl;
	}
	else
	{

		gzFile OUT_A;  		OUT_A = gzopen ((PaAA->OutPut1).c_str(), "wb");
		while(!(IN_A).eof())
		{

			int ii =0 ;
			map <string ,int> SeqIndex ;
			map <string ,int> ::iterator it ;

			for ( ; ii<Tarry &&  (!IN_A.eof() ); ii+=4)
			{
				getline(IN_A,FqData[ii]);
				getline(IN_A,FqData[ii+1]);
				getline(IN_A,FqData[ii+2]);
				getline(IN_A,FqData[ii+3]);
				if  (FqData[ii+1].length()==0)	{		 continue ;		}
				it=SeqIndex.find(FqData[ii+1]);
				if (it!=SeqIndex.end())
				{
					DupRead++;
					it->second=ii;
				}
				else
				{
					SeqIndex.insert(map <string ,int > ::value_type(FqData[ii+1],ii));
				}
				TotalRead++;
			}
			for (it=SeqIndex.begin(); it!=SeqIndex.end(); )
			{
				int II=it->second;
				gzprintf(OUT_A,"%s\n%s\n%s\n%s\n",FqData[II].c_str(),FqData[II+1].c_str(),FqData[II+2].c_str(),FqData[II+3].c_str());
				SeqIndex.erase(it++);
			}

		}
		gzclose(OUT_A);
		cerr<<"Total Read Number: "<<TotalRead<<"\n"<<"Remove Read Number "<<DupRead<<endl;
	}
	IN_A.close();
	delete [] FqData ;
	return 1;
}


int FqRmDupPE ( ParaClass * PaAA  )
{

	igzstream IN_A ((PaAA->InPut1).c_str(),ifstream::in); // ifstream  + gz
	igzstream IN_B ((PaAA->InPut2).c_str(),ifstream::in); // ifstream  + gz
	gzFile OUT_A;
	OUT_A = gzopen ((PaAA->OutPut1).c_str(), "wb");
	gzFile OUT_B;
	OUT_B = gzopen ((PaAA->OutPut2).c_str(), "wb");

	if(!IN_A.good())
	{
		cerr << "open IN File error: "<<(PaAA->InPut1)<<endl;
		return 1;
	}
	if(!IN_B.good())
	{
		cerr << "open IN File error: "<<(PaAA->InPut2)<<endl;
		return 1;
	}

	int TotalRead=0;
	int DupRead=0;

	int Tarry=(PaAA->InInt1)*4;
	string *FqDataA =new string[Tarry];
	string *FqDataB =new string[Tarry];
	cerr<<"Regard "<<(PaAA->InInt1)<<" read PE as a Blocks to remove Reads "<<endl;

	
	while(!IN_A.eof())
	{
			int ii =0 ;
			map <string ,int> SeqIndex ;
			map <string ,int> ::iterator it ;
			for ( ; ii<Tarry &&  (!IN_A.eof() ); ii+=4)
			{
				getline(IN_A,FqDataA[ii]);
				getline(IN_A,FqDataA[ii+1]);
				getline(IN_A,FqDataA[ii+2]);
				getline(IN_A,FqDataA[ii+3]);
				getline(IN_B,FqDataB[ii]);
				getline(IN_B,FqDataB[ii+1]);
				getline(IN_B,FqDataB[ii+2]);
				getline(IN_B,FqDataB[ii+3]);
				if  (FqDataA[ii+1].length()==0)	{ continue ;	}
				string CatSeq=FqDataA[ii+1]+FqDataB[ii+1];
				it=SeqIndex.find(CatSeq);
				if (it!=SeqIndex.end())
				{
					DupRead++;
					it->second=ii;
				}
				else
				{
					SeqIndex.insert(map <string ,int > ::value_type(CatSeq,ii));
				}
				TotalRead++;
			}
			for (it=SeqIndex.begin(); it!=SeqIndex.end();)
			{
				int II=it->second;
				gzprintf(OUT_A,"%s\n%s\n%s\n%s\n",FqDataA[II].c_str(),FqDataA[II+1].c_str(),FqDataA[II+2].c_str(),FqDataA[II+3].c_str());
				gzprintf(OUT_B,"%s\n%s\n%s\n%s\n",FqDataB[II].c_str(),FqDataB[II+1].c_str(),FqDataB[II+2].c_str(),FqDataB[II+3].c_str());
				SeqIndex.erase(it++);
			}

	}

	cout<<"Total Read PE Number: "<<TotalRead<<"\n"<<"Remove Read PE Number "<<DupRead<<endl;
	gzclose(OUT_A);
	gzclose(OUT_B);
	delete [] FqDataA ;
	delete [] FqDataB ;
	IN_B.close();
	IN_A.close();
	return 1;
}

int FQ_RmDup_main(int argc, char **argv)
//int main (int argc, char *argv[ ])
{

	ParaClass * PaAA_A24 = new  ParaClass ;
	PaAA_A24->InInt1=1024*1024*10;
	if( parse_cmd_rmDup(argc, argv, PaAA_A24)==0)
	{
		delete  PaAA_A24 ;
		return 1;
	}

	if  ((PaAA_A24->InPut2).empty())
	{
		FqRmDupSE( PaAA_A24  );
	}
	else
	{
		FqRmDupPE( PaAA_A24  ) ;
	}

	delete PaAA_A24 ;
	return 0;
}

#endif  // FQ_RmDup_






