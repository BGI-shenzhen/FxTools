#ifndef FQ_rmAdpter_H_
#define FQ_rmAdpter_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <math.h>
#include <zlib.h>
#include "FQ_Filter.h"

#define MAX_L 512
#define MAX_Q 128
using namespace std;

//  #define _A 0 #define _C 1    #define _G 2   #define _T 3   #define _N 4

int  print_usage_rmAdpter()
{
	cout <<""
		"\n"
		"Usage:rmAdapter -i A_1.fq A_1.fq -o B_1.fq B_2.fq [options]\n"
		"\n"
		"\t\t-i       <str>   File name of InFq1/InFq2 Input\n"
		"\t\t-o       <str>   OUT Fq File1/File2 name\n"
		"\n"
		"\t\t-a     <str>   Input adapters1 ListFile created by fqcheck\n"
		"\t\t               or seq [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]\n"
		"\t\t-b     <str>   Input adapters2 ListFile created by fqcheck\n"
		"\t\t               or seq [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA]\n"
		"\n"
		"\t\t-h             show this help\n" 
		"\n";
	return 1;
}



int parse_cmd_rmAdpter( int argc, char **argv, ParaClass * PaAA )
{
	if (argc <=4 ) {print_usage_rmAdpter();return 0;}

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
		else if (flag  ==  "AdaList1" ||  flag  ==  "a"	)
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			PaAA->InStr1=argv[i];
		}
		else if (flag  ==  "AdaList2" ||  flag  ==  "b" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			PaAA->InStr2=argv[i];
		}
		else if (flag  == "help")
		{
			print_usage_rmAdpter();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if ( (PaAA->InPut1).empty() || (PaAA->InPut2).empty() )
	{
		cerr<< "-InFq1   -InFq2 lack argument for the must"<<endl;
		return 0;
	}
	if ( (PaAA->OutPut2).empty() || (PaAA->OutPut1).empty() )
	{
		cerr<< "-OutFq1   -OutFq2 lack argument for the must"<<endl;
		return 0;
	}
	if ((access(PaAA->InStr1.c_str(), 0) == 0 )  ||  (access(PaAA->InStr2.c_str(), 0) == 0))
	{
		PaAA->TF=false;		
	}

	if ((PaAA->InStr1).empty()) { PaAA->InStr1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";}
	if ((PaAA->InStr2).empty()) { PaAA->InStr2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA";}

	(PaAA->OutPut1)=add_Asuffix(PaAA->OutPut1);
	(PaAA->OutPut2)=add_Asuffix(PaAA->OutPut2);

	return 1 ;
}



int RunFqRmAdpter( ParaClass * PaAA  )
{

	string adptid_AA, adpt_se_AA;
	string adptid_BB, adpt_se_BB;
	adpt_se_AA=PaAA->InStr1;
	adpt_se_BB=PaAA->InStr2;
	//	getline(ADF,adpt_id);	adpt_id.erase(0,1);
	//	getline(ADF,adpt_se);
	int adptlen_AA = adpt_se_AA.length();
	int adptlen_BB = adpt_se_BB.length();

	ogzstream OUT_A (PaAA->OutPut1.c_str());
	ogzstream OUT_B (PaAA->OutPut2.c_str());


	cout<<"adapter_id\tpolluted_reads\tempty_reads\tadapter_sequence"<<endl;


	//input file
	char read_id_AA[MAX_L], read_se_AA[MAX_L], qual_id_AA[MAX_L], qual_se_AA[MAX_L];
	int read_len_AA;

	gzFile read_fs_AA;
	read_fs_AA = gzopen (PaAA->InPut1.c_str(), "rb");
	if (read_fs_AA==NULL)
	{
		cerr<<"Cann't open files\t"<<PaAA->InPut1<<endl;
		return 1;
	}

	gzgets(read_fs_AA, read_id_AA, MAX_L);
	gzgets(read_fs_AA, read_se_AA, MAX_L);
	string read_AA = read_se_AA;
	read_len_AA = read_AA.length()-1;
	gzrewind(read_fs_AA);
	long read_num=0;
	long adpt_num=0;
	long empt_num=0;

	int Find_Len_AA=read_len_AA-9;

	char read_id_BB[MAX_L], read_se_BB[MAX_L], qual_id_BB[MAX_L], qual_se_BB[MAX_L];
	int read_len_BB;

	gzFile read_fs_BB;
	read_fs_BB = gzopen (PaAA->InPut2.c_str(), "rb");
	if (read_fs_BB==NULL)
	{
		cerr<<"Cann't open files\t"<<PaAA->InPut2<<endl;
		return 1;
	}

	gzgets(read_fs_BB, read_id_BB, MAX_L);
	gzgets(read_fs_BB, read_se_BB, MAX_L);
	string read_BB = read_se_BB;
	read_len_BB = read_BB.length()-1;
	gzrewind(read_fs_BB);

	int Find_Len_BB=read_len_BB-9;



	while (gzgets (read_fs_AA, read_id_AA, MAX_L)) 
	{
		gzgets (read_fs_AA, read_se_AA, MAX_L);
		gzgets (read_fs_AA, qual_id_AA, MAX_L);
		gzgets (read_fs_AA, qual_se_AA, MAX_L);

		gzgets (read_fs_BB, read_id_BB, MAX_L) ;
		gzgets (read_fs_BB, read_se_BB, MAX_L);
		gzgets (read_fs_BB, qual_id_BB, MAX_L);
		gzgets (read_fs_BB, qual_se_BB, MAX_L);


		read_num++;
		bool find=true;


		int a1 = adptlen_AA-15;
		int r1 = 0;
		int len;
		int mis;

		for ( ;a1>0 ; a1--)
		{
			int len1 = adptlen_AA - a1;
			len=(len1<read_len_AA)? len1:read_len_AA;
			mis = 0;
			int map[MAX_L];
			map[0]=0;
			for (int c=0; c<len; c++)
			{
				if (adpt_se_AA[a1+c]==read_se_AA[c]) {map[mis]++;}
				else {mis++;map[mis]=0;}
			}			

			int max_map=map[0];
			for (int c=1; c<=mis; c++)
			{
				if (map[c]>max_map) {max_map=map[c];}
			}

			if ( (max_map>14)    || (len>(mis*5)) )
			{
				find = false ;
				break;
			}
		}

		if (find)
		{
			for (r1=0; r1<Find_Len_AA; r1++)
			{
				int len2 = read_len_AA - r1;
				len=(adptlen_AA<len2)? adptlen_AA:len2;
				mis = 0;
				int map[MAX_L];
				map[0]=0;
				for (int c=0; c<len; c++)
				{
					if (adpt_se_AA[c]==read_se_AA[r1+c]) {map[mis]++;}
					else {mis++;map[mis]=0;}
				}

				int max_map=map[0];
				for (int c=1; c<=mis; c++)
				{
					if (map[c]>max_map) {max_map=map[c];}
				}
				if ( (max_map>14)   ||  (len>(mis*5))  )
				{
					find = false;
					break;
				}
			}

		}




		if  (!find)
		{	adpt_num++;
			if (r1<=3){empt_num++;}
			continue ;
		}
		else
		{



			a1 = adptlen_BB-15;
			r1 = 0;
			len=0;
			mis=0;

			for ( ;a1>0 ; a1--)
			{
				int len1 = adptlen_BB - a1;
				len=(len1<read_len_BB)? len1:read_len_BB;
				mis = 0;
				int map[MAX_L];
				map[0]=0;
				for (int c=0; c<len; c++)
				{
					if (adpt_se_BB[a1+c]==read_se_BB[c]) {map[mis]++;}
					else {mis++;map[mis]=0;}
				}

				int max_map=map[0];
				for (int c=1; c<=mis; c++)
				{
					if (map[c]>max_map) {max_map=map[c];}
				}

				if ( (max_map>14)    || (len>(mis*5)) )
				{
					find = false ;
					break;
				}
			}

			if (find)
			{
				for (r1=0; r1<Find_Len_BB; r1++)
				{
					int len2 = read_len_BB - r1;
					len=(adptlen_BB<len2)? adptlen_BB:len2;
					mis = 0;
					int map[MAX_L];
					map[0]=0;
					for (int c=0; c<len; c++)
					{
						if (adpt_se_BB[c]==read_se_BB[r1+c]) {map[mis]++;}
						else {mis++;map[mis]=0;}
					}

					int max_map=map[0];
					for (int c=1; c<=mis; c++)
					{
						if (map[c]>max_map) {max_map=map[c];}
					}
					if ( (max_map>14)   ||  (len>(mis*5))  )
					{
						find = false;
						break;
					}
				}

			}



			if  (!find)
			{	adpt_num++;
				if (r1<=3){empt_num++;}
				continue ;
			}
			else
			{
				OUT_A<<read_id_AA<<read_se_AA<<qual_id_AA<<qual_se_AA<<"\n";
				OUT_B<<read_id_BB<<read_se_BB<<qual_id_BB<<qual_se_BB<<"\n";
			}



		}

	}
	gzclose(read_fs_AA);
	gzclose(read_fs_BB);
	OUT_A.close();
	OUT_B.close();


	float adpt_percent = ((float) adpt_num)/read_num*100;
	float empt_percent = ((float) empt_num)/read_num*100;
	cout<<"Total adpt Num:"<<adpt_num<<" ("<<setprecision(3)<<fixed<<adpt_percent<<"%)"<<'\t'<<empt_num<<" ("<<empt_percent<<"%)"<<endl;
	cout<<"\ntotal_reads: "<<read_num<<endl;



	return 0;
}



int RunFqRmAdpterList( ParaClass * PaAA  )
{
	igzstream IN_1 ((PaAA->InStr1).c_str(),ifstream::in); 
	igzstream IN_2 ((PaAA->InStr2).c_str(),ifstream::in);
	map <string ,bool> ListA_1; 
	map <string ,bool> ListA_2;
	string readID;
	string temp; getline(IN_1,temp); getline(IN_2,temp);
	while(!IN_1.eof())
	{
		string line ;
		getline(IN_1,line);
		if (line.length()<=0)  { continue  ; }
		istringstream isone (line,istringstream::in);
		isone>>readID;
		ListA_1.insert(map <string ,bool>  :: value_type(readID,true));
	}
	IN_1.close();

	while(!IN_2.eof())
	{
		string line ;
		getline(IN_2,line);
		if (line.length()<=0)  { continue  ; }
		istringstream isone (line,istringstream::in);
		isone>>readID;
		ListA_2.insert(map <string ,bool>  :: value_type(readID,true));
	}
	IN_2.close();

	igzstream IN_A ((PaAA->InPut1).c_str(),ifstream::in); // ifstream  + gz
	igzstream IN_B ((PaAA->InPut2).c_str(),ifstream::in); // ifstream  + gz
	ogzstream OUT_A ((PaAA->OutPut1).c_str());
	ogzstream OUT_B ((PaAA->OutPut2).c_str());

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
	if(!OUT_A.good())
	{
		cerr << "open OUT File error: "<<(PaAA->OutPut1)<<endl;
		return 1;
	}
	if(!OUT_B.good())
	{
		cerr << "open OUT File error: "<<(PaAA->OutPut2)<<endl;
		return 1;
	}


	while(!IN_A.eof())
	{
		string ID_1,seq_1,temp_1,Quly_1;
		string ID_2,seq_2,temp_2,Quly_2;
		getline(IN_A,ID_1);
		getline(IN_A,seq_1);
		getline(IN_A,temp_1);
		getline(IN_A,Quly_1);
		getline(IN_B,ID_2);
		getline(IN_B,seq_2);
		getline(IN_B,temp_2);
		getline(IN_B,Quly_2);
		if (ID_1.length()<=0)  { continue  ; }
		string AA=ID_1.substr(1);
		string BB=ID_2.substr(1);
		map <string, bool> :: iterator itInt1= ListA_1.find(AA);
		map <string, bool> :: iterator itInt2= ListA_2.find(BB);
		if  (itInt1== ListA_1.end() ||  itInt2 ==ListA_2.end())
		{
			OUT_A<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";
			OUT_B<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n";
		}
	}

	OUT_B.close();
	OUT_A.close();
	IN_A.close();
	IN_B.close();


}

int FQ_RmAdpter_main(int argc, char **argv)
//int main (int argc, char *argv[ ])
{

	ParaClass * PaAA_A24 = new  ParaClass ;

	if( parse_cmd_rmAdpter(argc, argv, PaAA_A24)==0)
	{
		delete  PaAA_A24 ;
		return 1;
	}


	if  (PaAA_A24->TF)
	{
		RunFqRmAdpter( PaAA_A24  );
	}
	else
	{
		RunFqRmAdpterList( PaAA_A24  ) ;
	}
	delete PaAA_A24 ;
	return 0;
}

#endif  // FQ_check_H_






