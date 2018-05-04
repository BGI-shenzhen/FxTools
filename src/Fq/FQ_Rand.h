#ifndef FQ_Rand_H_
#define FQ_Rand_H_

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

int  print_usage_rand()
{
	cout <<""
		"\n"
		"Usage:rand  -i A1.fq A2.fq -o B1.fq B2.fq [options]\n"
		"\n"
		"\t\t-i        <str>   File name of InFq1/InFq2 Input\n"
		"\t\t-o        <str>   OUT Fq File1/File2,otherwise[STDOUT]\n"
		"\n"
		"\t\t-p      <float>   read rand out proportion [0.1]\n"
		"\t\t-s        <int>   rand seed,default time\n"
		"\n"
		"\t\t-h                show this help\n" 
		"\n";
	return 1;
}



int parse_cmd_rand( int argc, char **argv, ParaClass * PaAA )
{
	if (argc <=2) {print_usage_rand();return 0;}

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
		else if (flag  ==  "proportion"  ||  flag  ==  "p" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			PaAA->Infloat1=atof(argv[i]);
		}
		else if (flag  ==  "seed"   ||  flag  ==  "s" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			PaAA->InInt1=atoi(argv[i]);
		}
		else if (flag  == "help"  ||  flag  ==  "h")
		{
			print_usage_rand();return 0;
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


int FqRandSE ( ParaClass * PaAA  )
{	

	short int Y=(short int)(PaAA->Infloat1*100);
	short int X ;
	char buf[8196];
	FILE *fp ;
	string path=(PaAA->InPut1);
	if  (GzipTrue(path))
	{
		path="gzip -cd "+ (PaAA->InPut1);
	}
	else
	{
		path  = "cat "+(PaAA->InPut1) ;
	}
	sprintf(buf, path.c_str());
	fp = popen(buf, "r");
	if( fp == NULL )
	{
		printf("popen() error!/n");
		return -1;
	}


	if ((PaAA->OutPut1).empty())
	{
		while(fgets(buf, sizeof(buf), fp))
		{
			X=(rand()%100);
			if (X>Y)  
			{
				fgets(buf, sizeof(buf), fp);
				fgets(buf, sizeof(buf), fp);
				fgets(buf, sizeof(buf), fp);
				continue;
			}
			else
			{
				fprintf(stdout,"%s",buf);
				fgets(buf, sizeof(buf), fp); fprintf(stdout,"%s",buf);
				fgets(buf, sizeof(buf), fp); fprintf(stdout,"%s",buf);
				fgets(buf, sizeof(buf), fp); fprintf(stdout,"%s",buf);
			}
		}
	}
	else
	{
		gzFile OUT_A;
		OUT_A = gzopen ((PaAA->OutPut1).c_str(), "wb");
		while(fgets(buf, sizeof(buf), fp))
		{
			X=(rand()%100);
			if (X>Y)  
			{
				fgets(buf, sizeof(buf), fp);
				fgets(buf, sizeof(buf), fp);
				fgets(buf, sizeof(buf), fp);
				continue;
			}
			else
			{
				gzprintf(OUT_A,"%s",buf);
				fgets(buf, sizeof(buf), fp); gzprintf(OUT_A,"%s",buf);
				fgets(buf, sizeof(buf), fp); gzprintf(OUT_A,"%s",buf);
				fgets(buf, sizeof(buf), fp); gzprintf(OUT_A,"%s",buf);
			}
		}
		gzclose(OUT_A);
	}
	pclose(fp);
	return 1;
}


int FqRandPE ( ParaClass * PaAA  )
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

	string seq_1,temp_1,Quly_1;
	string ID_2,seq_2,temp_2,Quly_2;
	int Y=int(PaAA->Infloat1*100)-1;
	int X=0;
	while(!IN_A.eof())
	{
		string ID_1;
		getline(IN_A,ID_1);
		getline(IN_A,seq_1);
		getline(IN_A,temp_1);
		getline(IN_A,Quly_1);
		getline(IN_B,ID_2);
		getline(IN_B,seq_2);
		getline(IN_B,temp_2);
		getline(IN_B,Quly_2);

		X=(rand()%100);
		if (X>Y)  { continue  ; }
		if (ID_1.length()<=0)  { continue  ; }
		gzprintf(OUT_A,"%s\n%s\n%s\n%s\n",ID_1.c_str(),seq_1.c_str(),temp_1.c_str(),Quly_1.c_str());
		gzprintf(OUT_B,"%s\n%s\n%s\n%s\n",ID_2.c_str(),seq_2.c_str(),temp_2.c_str(),Quly_2.c_str());
	}

	gzclose(OUT_A);
	gzclose(OUT_B);
	IN_B.close();
	IN_A.close();
}

int FQ_Rand_main(int argc, char **argv)
//int main (int argc, char *argv[ ])
{

	ParaClass * PaAA_A24 = new  ParaClass ;
	PaAA_A24->Infloat1=0.1;
	if( parse_cmd_rand(argc, argv, PaAA_A24)==0)
	{
		delete  PaAA_A24 ;
		return 1;
	}

	if (PaAA_A24->InInt1==0)
	{
		srand((unsigned)time(NULL));
	}
	else
	{
		srand(PaAA_A24->InInt1);
	}

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);



	if  ((PaAA_A24->InPut2).empty())
	{
		FqRandSE( PaAA_A24  );
	}
	else
	{
		FqRandPE( PaAA_A24  ) ;
	}

	delete PaAA_A24 ;
	return 0;
}

#endif  // FQ_Rand_






