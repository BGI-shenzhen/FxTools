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
#include "FQ_Filter.h"
#include <sys/stat.h>

//KSEQ_AINIT(gzFile, gzread)

using namespace std;
typedef long long  llong ;

int  print_Ausage_RG()
{
	cout <<""
		"\n"
		"\tUsage: splitpool -InFq1 <In1.fq> -InFq2 <In2.fq> -Index <sample.info> -Flag <ferment.seq>\n"
		"\n"
		"\t\t-InFq1    <str>   Input #_1.Fq to split RAD(GBS)\n"
		"\t\t-InFq2    <str>   Input #_2.Fq to split RAD(GBS)\n"
		"\t\t-Index    <str>   Input File with (sample seq)\n"
		"\t\t-Flag     <str>   Input File with Flag(ferment) seq\n"
		"\n"
		"\t\t-OutDir   <str>   Output Dir for Split Files[PWD]\n"
		"\t\t-MisMatch         Allow one misMatch on the sample seq\n"
		"\n"
		"\t\t-help             show this help\n"
		"\n";
	return 1;
}

int parse_Acmd_RG(int argc, char **argv , Para_A24 * para_A24 )
{
	if (argc <=2 ) {print_Ausage_RG();return 0;}

	for(int i = 1; i < argc ; i++)
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
			if(i + 1 == argc) {LogLackArg( flag ) ; return 0;}
			i++;
			para_A24->InFq1=argv[i];
		}
		else if (flag  == "InFq2" )
		{
			if(i + 1 == argc) {LogLackArg( flag ) ; return 0;}
			i++;
			para_A24->InFq2=argv[i];
		}
		else if (flag  == "Index" )
		{
			if(i + 1 == argc) {LogLackArg( flag ) ; return 0;}
			i++;
			para_A24->adapter1=argv[i];
		}
		else if (flag  == "Flag" )
		{
			if(i + 1 == argc) {LogLackArg( flag ) ; return 0;}
			i++;
			para_A24->adapter2=argv[i];
		}
		else if (flag  ==  "OutDir")
		{
			if(i + 1 == argc) {LogLackArg( flag ) ; return 0;}
			i++;
			para_A24->OutFq1=argv[i];
			para_A24->OutFq1=(para_A24->OutFq1)+"/";
		}
		else if (flag  ==  "MisMatch")
		{
			para_A24->indexcut=1;
		}
		else if (flag  == "help")
		{
			print_Ausage_RG();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para_A24->InFq1).empty() ||  (para_A24->InFq2).empty() ||  (para_A24->adapter1).empty()  ||  (para_A24->adapter2).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	if ((para_A24->OutFq1).empty())
	{
		(para_A24->OutFq1)="./" ;
	}
	return 1 ;
}


int FQ_Split_Pooing_main(int argc, char *argv[])
	//	int main(int argc, char *argv[])
{
	Para_A24 * para_A24 = new Para_A24;

	if( parse_Acmd_RG(argc, argv, para_A24 )==0)
	{
		delete  para_A24 ;
		return 1;
	}

	mkdir((para_A24->OutFq1).c_str() , 0755 );
	igzstream  INA ((para_A24->InFq1).c_str(),ifstream::in);

	if(!INA.good())
	{
		cerr << "open InputFile error: "<<(para_A24->InFq1)<<endl;
		delete  para_A24 ; return 1;
	}

	igzstream  INB ((para_A24->InFq2).c_str(),ifstream::in);
	if(!INB.good())
	{
		cerr << "open InputFile error: "<<(para_A24->InFq2)<<endl;
		delete  para_A24 ; return 1;
	}

	igzstream INDEX ((para_A24->adapter1).c_str(),ifstream::in);
	if(!INDEX.good())
	{
		cerr << "open InputFile error: "<<(para_A24->adapter1)<<endl;
		delete  para_A24 ; return 1;
	}

	igzstream MEI ((para_A24->adapter2).c_str(),ifstream::in);
	if(!(MEI.good()))
	{
		cerr << "open InputFile error: "<<(para_A24->adapter2)<<endl;
		delete  para_A24 ; return 1;
	}

	int MEI_len=0;
	map <string,bool> HashMEI ;
	while(!MEI.eof())
	{
		string  line ;
		getline(MEI,line);
		if (line.length()<=0)  { continue  ; }
		else if (MEI_len<line.length())
		{
			MEI_len=line.length();
		}
		HashMEI[line]=true;
	}
	MEI.close();

	int sample_count=0;
	map <string,bool> CheckSampleCount;
	map <string,bool> :: iterator map_it ;
	while(!INDEX.eof())
	{
		string  line ,sample_ID,Seq;
		getline(INDEX,line);
		if (line.length()<=0)  { continue  ;}
		sample_count++;

		istringstream isone (line,istringstream::in);
		isone>>sample_ID>>Seq;
		map_it=CheckSampleCount.find(sample_ID);
		if (map_it!=CheckSampleCount.end())
		{
			cerr<<"Warming : sample double in this INDEX Files. Sample ID: "<<sample_ID<<"; please renamed it diff\n";
			delete para_A24 ;
			return 1;
		}
		else
		{
			CheckSampleCount.insert(map <string,int>::value_type(sample_ID,true));
		}
	}
	INDEX.close();INDEX.clear();


	ogzstream *OUT1 = new ogzstream[sample_count] ;
	ogzstream *OUT2 = new ogzstream[sample_count] ;
	sample_count=0;
	int SeqMaxLen=0;
	map <string,int> Seq2Sample ;

	igzstream INDEX2 ((para_A24->adapter1).c_str(),ifstream::in);	
	if(!INDEX2.good())
	{
		cerr << "open InputFile error: "<<(para_A24->adapter1)<<endl;
		delete  para_A24 ; return 1;
	}

	if ((para_A24->indexcut)==1)
	{
		char Base[5];
		Base[0]='A';  Base[1]='C';  Base[2]='T';  Base[3]='G';  Base[4]='N';
		map <string,int> MapTemp ;
		map <string,int> Seq2Sample_Rm ;
		while(!INDEX2.eof())
		{
			string  line ,sample_ID,Seq;
			getline(INDEX2,line);
			if (line.length()<=0)  { continue; }
			istringstream isone (line,istringstream::in);
			isone>>sample_ID>>Seq;
			map <string, int> :: iterator it=MapTemp.find(Seq);
			if (it == MapTemp.end())
			{
				MapTemp.insert(map <string, int> :: value_type(Seq,sample_count));
			}
			else
			{
				cerr<<"Warning: the seq "<<Seq<<"\thad two sample: "<<sample_ID<<" and the "<<MapTemp[Seq]<<endl;
				cerr<<"plese check it"<<endl;
			}

			int seq_len=Seq.length() ;
			if (SeqMaxLen<seq_len)
			{
				SeqMaxLen=seq_len;
			}
			for (int i=0 ; i<seq_len ; i++)
			{
				for (int j=0  ;j<5;  j++)
				{
					string new_seq=Seq ;
					new_seq[i]=Base[j];
					map <string, int> :: iterator it=Seq2Sample.find(new_seq);
					if (it == Seq2Sample.end())
					{
						Seq2Sample.insert(map <string, int> :: value_type(new_seq,sample_count));
					}
					else
					{
						it=MapTemp.find(new_seq);
						if (it == MapTemp.end())
						{
							cerr<<"Warning: the seq "<<new_seq<<"\thad two sample: "<<sample_ID<<" and the "<<Seq2Sample[new_seq]<<",so we classify them to the nuknow"<<endl;
							it=Seq2Sample.find(new_seq);
							Seq2Sample.erase(it);
							Seq2Sample_Rm.insert(map <string, int> :: value_type(new_seq,sample_count));
						}
					}
				}
			}

			string OUT_name1=(para_A24->OutFq1)+sample_ID+"_1.fq.gz";
			OUT1[sample_count].open(OUT_name1.c_str());	
			string OUT_name2=(para_A24->OutFq1)+sample_ID+"_2.fq.gz";
			OUT2[sample_count].open(OUT_name2.c_str());
			sample_count++;		
		}

		map <string, int> :: iterator map_it_rm=Seq2Sample_Rm.begin();
		while(map_it_rm!=Seq2Sample_Rm.end())
		{
			string new_seq_rm=map_it_rm->first ;
			map <string, int> :: iterator map_it_Check=Seq2Sample.find(new_seq_rm);
			if (map_it_Check!=Seq2Sample.end())
			{
				Seq2Sample.erase(map_it_Check);
			}
			map_it_rm++;
		}
	}
	else
	{
		while(!INDEX2.eof())
		{
			string  line ,sample_ID,Seq;
			getline(INDEX2,line);
			if (line.length()<=0)  { continue; }
			istringstream isone (line,istringstream::in);
			isone>>sample_ID>>Seq;
			int seq_len=Seq.length() ;
			if (SeqMaxLen<seq_len)
			{
				SeqMaxLen=seq_len;
			}

			map <string, int> :: iterator it=Seq2Sample.find(Seq);
			if (it == Seq2Sample.end())
			{
				Seq2Sample.insert(map <string, int> :: value_type(Seq,sample_count));
			}
			else
			{
				cerr<<"Warning: the seq "<<Seq<<"\thad two sample: "<<sample_ID<<" and the "<<Seq2Sample[Seq]<<endl;
				cerr<<"plese check it"<<endl;
			}
			string OUT_name1=(para_A24->OutFq1)+sample_ID+"_1.fq.gz";
			OUT1[sample_count].open(OUT_name1.c_str());	
			string OUT_name2=(para_A24->OutFq1)+sample_ID+"_2.fq.gz";
			OUT2[sample_count].open(OUT_name2.c_str());
			sample_count++;
		}
	}
	INDEX2.close();INDEX2.clear();



	string OUT_name1=(para_A24->OutFq1)+"UnKnow_1.fq.gz";
	string OUT_name2=(para_A24->OutFq1)+"UnKnow_2.fq.gz";
	ogzstream OUT_U1 (OUT_name1.c_str());
	ogzstream OUT_U2 (OUT_name2.c_str());

	int End_find=SeqMaxLen+MEI_len;
	map<string,bool>::iterator innerit ;
	map <string, int> :: iterator it ;
	string  ID_2, Seq_2, Flag_2,Q_2 ,  Seq_1, Flag_1,Q_1;
	bool A=false ;	size_t find_position=0;
	int ID=1;
	while(!INA.eof())
	{
		string  ID_1 ;
		getline(INA,ID_1);
		if (ID_1.length()<=0)  { continue ;}
		getline(INA,Seq_1);
		getline(INA,Flag_1);
		getline(INA,Q_1);
		getline(INB,ID_2);
		getline(INB,Seq_2);
		getline(INB,Flag_2);
		getline(INB,Q_2);
		A=false ; find_position=0;
		size_t    Flag_position=0;
		for (  innerit=HashMEI.begin();  innerit!=HashMEI.end();  innerit++)
		{
			string subMEIseq=innerit->first;
			Flag_position=Seq_1.find(subMEIseq,3);
			if ( (Flag_position != string::npos) && (Flag_position <End_find ) )
			{				
				A=true;
				if (find_position==0)
				{
					find_position=Flag_position;
				}
				else if (find_position>Flag_position)
				{
					find_position=Flag_position;
				}
			}
		}
		if (A)
		{
			string BB=Seq_1.substr(0,find_position) ; 
			it=Seq2Sample.find(BB);
			if (it!=Seq2Sample.end())
			{
				ID=it->second;
				Q_1=Q_1.substr(find_position); Seq_1=Seq_1.substr(find_position);
				OUT1[ID]<<ID_1<<"\n"<<Seq_1<<"\n"<<Flag_1<<"\n"<<Q_1<<endl;
				OUT2[ID]<<ID_2<<"\n"<<Seq_2<<"\n"<<Flag_2<<"\n"<<Q_2<<endl;
			}
			else
			{
				OUT_U1<<ID_1<<"\n"<<Seq_1<<"\n"<<Flag_1<<"\n"<<Q_1<<endl;
				OUT_U2<<ID_2<<"\n"<<Seq_2<<"\n"<<Flag_2<<"\n"<<Q_2<<endl;
			}
		}
		else
		{
			OUT_U1<<ID_1<<"\n"<<Seq_1<<"\n"<<Flag_1<<"\n"<<Q_1<<endl;
			OUT_U2<<ID_2<<"\n"<<Seq_2<<"\n"<<Flag_2<<"\n"<<Q_2<<endl;
		}
	}
	INA.close();
	INB.close();
	for (int i=0 ; i<sample_count ; i++)
	{
		OUT1[i].close();
		OUT2[i].close();
	}
	OUT_U1.close(); OUT_U2.close();

	delete para_A24 ;
	return 0;
}

///////// swimming in the sky and flying in the sea ////////////

