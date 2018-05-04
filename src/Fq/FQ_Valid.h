#ifndef FQ_Valid_H_
#define FQ_Valid_H_

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



int  print_U_Valid()
{
	cout <<""
		"\n"
		"\tUsage: valid -i <in.fq>  \n"
		"\n"
		"\t\t-i    <str>     Input FASTQ file\n"
		"\n"
		"\t\t-o    <str>     Output a new valid FASTQ by removing the reads with error\n"
		"\t\t-w    <str>     ID of reads with error, default [STDOUT]\n" 
		"\n"
		"\t\t-h              show this help\n" 
		"\n";
	return 1;
}


int parse_Vaild_A30(int argc, char **argv,In3str1v *  para_A30 )
{
	if (argc <=2 ) {print_U_Valid();return 0 ;}
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
			print_U_Valid();return 0 ;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0 ;
		}
	}
	if  (para_A30->InStr1.empty()  )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0 ;
	}
	if (para_A30->InStr3.empty() )
	{
		para_A30->InStr3="./.tmpValid.fq.gz";
	}
	else
	{
		para_A30->InStr3=add_Asuffix(para_A30->InStr3);
	}
	return 1 ;
}


int FQ_Valid_main( int argc,char *argv[] )
//int main (int argc,char *argv[])
{
	In3str1v *  para_A30 =new In3str1v ;
	if( parse_Vaild_A30( argc, argv,para_A30) ==0 )
	{
		delete  para_A30 ;
		return 1 ;
	}

	igzstream IN_A ((para_A30->InStr1).c_str(),ifstream::in); // ifstream  + gz
	ogzstream OUT ((para_A30->InStr3).c_str());
	bool OUTTF=true;
	bool OUTTFALL=true;
	bool OUTTFALLV2=true;

	string ID_1,ID_2;
	getline(IN_A,ID_2);

	if (!(para_A30->InStr2).empty())
	{
		ofstream Info ((para_A30->InStr2).c_str());
		while(!IN_A.eof())
		{
			OUTTF=true;
			string seq_1,temp_1,Quly_1;
			ID_1=ID_2;
			if (ID_1[0]!='@')
			{
				Info<<"wrong ID nearby at "<<ID_1<<endl;
				OUTTF=false ;
			}

			getline(IN_A,temp_1);
			while((temp_1[0]!='+') &&  (!IN_A.eof()))
			{
				seq_1=seq_1+temp_1;
				getline(IN_A,temp_1);
			}

			if (temp_1.length()>1)
			{
				temp_1[0]='@';
				if (ID_1!=temp_1)
				{
					Info<<"Sequence(s) with Error with Diff Read ID :"<<ID_1<<"\t\t"<<temp_1<<endl;
					OUTTF=false ;
				}
				temp_1="+";
			}

			getline(IN_A,Quly_1);
			getline(IN_A,ID_2);
			string::size_type SeqLen=seq_1.length();

			while( ((ID_2[0]!='@')  ||  (Quly_1.length()<SeqLen) )  && (!IN_A.eof())   )
			{
				Quly_1=Quly_1+ID_2;
				getline(IN_A,ID_2);
			}


			string::size_type  pos(0);
			pos=0;
			if ((pos=seq_1.find(" "))!=string::npos )
			{
				Info<<"Sequence(s) with space ascii, Read ID :"<<ID_1<<endl;
				OUTTF=false ;
			}

			pos=0;
			if ((pos=seq_1.find("\t"))!=string::npos )
			{
				Info<<"Sequence(s) with Tab ascii, Read ID :"<<ID_1<<endl;
				OUTTF=false ;
			}

			pos=0;
			if ((pos=Quly_1.find(" "))!=string::npos )
			{
				Info<<"Sequence(s) with space ascii, Read ID :"<<ID_1<<endl;
				OUTTF=false ;
			}
			pos=0;
			if ((pos=Quly_1.find("\t"))!=string::npos )
			{
				Info<<"Sequence(s) with Tab ascii, Read ID :"<<ID_1<<endl;
				OUTTF=false ;
			}


			string::size_type s_Aend=Quly_1.length();

			if (SeqLen!=s_Aend)
			{
				Info<<"Sequence(s) with Error at seq and quality length Diff :"<<ID_1<<endl;			
				OUTTF=false ;
			}
			if (s_Aend<1)
			{
				Info<<"Sequence(s) with Error at Empty quality :"<<ID_1<<endl;
				OUTTF=false;
			}
			if  (SeqLen<1)
			{
				Info<<"Sequence(s) with Error at Empty seq :"<<ID_1<<endl;
				OUTTF=false;
			}

			string::size_type ix=ID_1.find("PHRED scores");


			if ( ix==string::npos )
			{
				ix=ID_1.find("Solexa scores") ;
				if (ix!=string::npos )
				{
					ix++;
				}
			}
			if ( ix!=string::npos ) // PHRED
			{
				OUTTFALL=false;
				string Tmp=ID_1.substr(ix+13);
				vector<string> inf;
				split(Tmp,inf," \t");
				vector <int > Num ;
				if (inf[0]== "from")
				{
					int Stat=atoi(inf[1].c_str());
					int End=atoi(inf[3].c_str());
					if  (Stat<End)
					{
						for (int jj=Stat ; jj<=End ; jj++)
						{
							Num.push_back(jj);
						}
					}
					else
					{
						for (int jj=Stat ; jj>=End; jj--)
						{
							Num.push_back(jj);
						}
					}
				}
				else
				{
					for (int ii=1 ; ii<inf.size(); ii++)
					{
						if( (ix=inf[ii].find(','))!=string::npos)
						{
							inf[ii].replace(ix,1,"");
							int Stat=atoi(inf[ii].c_str());
							Num.push_back(Stat);
						}
						else
						{
							int Stat=atoi(inf[ii].c_str());
							Num.push_back(Stat);
							break;
						}
					}
				}

				int CCbin=Num.size();
				int Tbintt=0;
				for(ix=0 ; ix<s_Aend ; ix++)
				{
					Tbintt=int(ix % CCbin);
					int cc=Num[Tbintt]+33;
					if (cc>73)
					{
						cc=73;
					}
					else if (cc<33)
					{
						cc=33;
					}
					Quly_1[ix]=cc;
				}
				istringstream isone (ID_1,istringstream::in);
				isone>>ID_1;
				//			ID_1=Tmp;
			}
			for(ix=0 ; ix<s_Aend ; ix++)
			{
				if  (Quly_1[ix]<33 ||  Quly_1[ix]>115)
				{
					Info<<"Sequence(s) with Error at too Low/Hight quality ASCII : "<<ID_1<<"\t at site "<<ix<<endl;
					OUTTF=false ;
				}
			}
			if (OUTTF)
			{
				OUT<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";
			}
			else
			{
				OUTTFALLV2=false;
			}
		}

		if (ID_1!=ID_2 && (!ID_2.empty()))
		{
			Info<<"Sequence(s) with Error at Empty seq :"<<ID_2<<endl;
			OUTTF=false;
			OUTTFALLV2=false;
		}	

		if (OUTTFALLV2)
		{
			if (OUTTFALL)
			{
				Info<<"\t\t\tVALID FASTQ File"<<endl;
				cout<<"\t\t\tVALID FASTQ File"<<endl;
			}
			else
			{
				if ((para_A30->InStr3) ==  "./.tmpValid.fq.gz" ) 
				{
					Info<<"\t\t\tVALID FASTQ File. Note:File is not in normal format, please check your file or add -o to generate a new normal FASTQ"<<endl;
					cout<<"\t\t\tVALID FASTQ File. Note:File is not in normal format, please check your file or add -o to generate a new normal FASTQ"<<endl;
				}
				else
				{
					Info<<"\t\t\tVALID FASTQ File, Note:File is not in normal format, we had change it to the normal Fq "<<(para_A30->InStr3)<<endl;
					cout<<"\t\t\tVALID FASTQ File, Note:File is not in normal format, we had change it to the normal Fq "<<(para_A30->InStr3)<<endl;
				}
			}
		}
		else
		{
			if ((para_A30->InStr3) ==  "./.tmpValid.fq.gz" )
			{
				cout<<"\t\t\tINVALID FASTQ File ; please check your file or add -o to generate a new normal FASTQ"<<endl;
				Info<<"\t\t\tINVALID FASTQ File ; please check your file or add -o to generate a new normal FASTQ"<<endl;
			}
			else
			{
				cout<<"\t\t\tINVALID FASTQ File"<<endl;
				Info<<"\t\t\tINVALID FASTQ File"<<endl;
			}
		}

		Info.close();
	}
	else
	{

		if ((para_A30->InStr3) ==  "./.tmpValid.fq.gz" )
		{


			while(!IN_A.eof())
			{
				OUTTF=true;
				string seq_1,temp_1,Quly_1;
				ID_1=ID_2;
				if (ID_1[0]!='@')
				{
					OUTTF=false ;
				}

				getline(IN_A,temp_1);
				while((temp_1[0]!='+') &&  (!IN_A.eof()))
				{
					seq_1=seq_1+temp_1;
					getline(IN_A,temp_1);
				}

				if (temp_1.length()>1)
				{
					temp_1[0]='@';
					if (ID_1!=temp_1)
					{
						OUTTF=false ;
					}
					temp_1="+";
				}

				getline(IN_A,Quly_1);
				getline(IN_A,ID_2);
				string::size_type SeqLen=seq_1.length();

				while( ((ID_2[0]!='@')  ||  (Quly_1.length()<SeqLen) )  && (!IN_A.eof())   )
				{
					Quly_1=Quly_1+ID_2;
					getline(IN_A,ID_2);
				}

				string::size_type  pos(0);
				pos=0;
				if ((pos=seq_1.find(" "))!=string::npos )
				{
					OUTTF=false ;
				}

				pos=0;
				if ((pos=seq_1.find("\t"))!=string::npos )
				{
					OUTTF=false ;
				}

				pos=0;
				if ((pos=Quly_1.find(" "))!=string::npos )
				{
					OUTTF=false ;
				}
				pos=0;
				if ((pos=Quly_1.find("\t"))!=string::npos )
				{
					OUTTF=false ;
				}


				string::size_type s_Aend=Quly_1.length();

				if (SeqLen!=s_Aend)
				{
					OUTTF=false ;
				}
				if (s_Aend<1)
				{
					OUTTF=false;
				}
				if  (SeqLen<1)
				{
					OUTTF=false;
				}

				string::size_type ix=ID_1.find("PHRED scores");


				if ( ix==string::npos )
				{
					ix=ID_1.find("Solexa scores") ;
					if (ix!=string::npos )
					{
						ix++;
					}
				}
				if ( ix!=string::npos ) // PHRED
				{
					OUTTFALL=false;
					string Tmp=ID_1.substr(ix+13);
					vector<string> inf;
					split(Tmp,inf," \t");
					vector <int > Num ;
					if (inf[0]== "from")
					{
						int Stat=atoi(inf[1].c_str());
						int End=atoi(inf[3].c_str());
						if  (Stat<End)
						{
							for (int jj=Stat ; jj<=End ; jj++)
							{
								Num.push_back(jj);
							}
						}
						else
						{
							for (int jj=Stat ; jj>=End; jj--)
							{
								Num.push_back(jj);
							}
						}
					}
					else
					{
						for (int ii=1 ; ii<inf.size(); ii++)
						{
							if( (ix=inf[ii].find(','))!=string::npos)
							{
								inf[ii].replace(ix,1,"");
								int Stat=atoi(inf[ii].c_str());
								Num.push_back(Stat);
							}
							else
							{
								int Stat=atoi(inf[ii].c_str());
								Num.push_back(Stat);
								break;
							}
						}
					}

					int CCbin=Num.size();
					int Tbintt=0;
					for(ix=0 ; ix<s_Aend ; ix++)
					{
						Tbintt=int(ix % CCbin);
						int cc=Num[Tbintt]+33;
						if (cc>73)
						{
							cc=73;
						}
						else if (cc<33)
						{
							cc=33;
						}
						Quly_1[ix]=cc;
					}
					istringstream isone (ID_1,istringstream::in);
					isone>>ID_1;
					//			ID_1=Tmp;
				}
				for(ix=0 ; ix<s_Aend ; ix++)
				{
					if  (Quly_1[ix]<33 ||  Quly_1[ix]>115)
					{
						OUTTF=false ;
					}
				}
				if (OUTTF)
				{

				}
				else
				{
					OUTTFALLV2=false;
					break;
				}
			}

			if (ID_1!=ID_2 && (!ID_2.empty()))
			{
				OUTTF=false;
				OUTTFALLV2=false;
			}

			if (OUTTFALLV2)
			{
				if  (OUTTFALL)
				{
					cout<<"\t\t\tVALID FASTQ File"<<endl;
				}
				else
				{
					cout<<"\t\t\tVALID FASTQ File. Note:File is not in normal format, please check your file or add -o to generate a new normal FASTQ"<<endl;
				}

			}
			else
			{
				cout<<"\t\t\tINVALID FASTQ File; you may add Para -o to generate the new valid fq"<<endl;
			}











		}
		else
		{

			while(!IN_A.eof())
			{
				OUTTF=true;
				string seq_1,temp_1,Quly_1;
				ID_1=ID_2;
				if (ID_1[0]!='@')
				{
					OUTTF=false ;
				}

				getline(IN_A,temp_1);
				while((temp_1[0]!='+') &&  (!IN_A.eof()))
				{
					seq_1=seq_1+temp_1;
					getline(IN_A,temp_1);
				}

				if (temp_1.length()>1)
				{
					temp_1[0]='@';
					if (ID_1!=temp_1)
					{
						OUTTF=false ;
					}
					temp_1="+";
				}

				getline(IN_A,Quly_1);
				getline(IN_A,ID_2);
				string::size_type SeqLen=seq_1.length();

				while( ((ID_2[0]!='@')  ||  (Quly_1.length()<SeqLen) )  && (!IN_A.eof())   )
				{
					Quly_1=Quly_1+ID_2;
					getline(IN_A,ID_2);
				}

				string::size_type  pos(0);

				pos=0;
				if ((pos=seq_1.find(" "))!=string::npos )
				{
					OUTTF=false ;
				}

				pos=0;
				if ((pos=seq_1.find("\t"))!=string::npos )
				{
					OUTTF=false ;
				}

				pos=0;
				if ((pos=Quly_1.find(" "))!=string::npos )
				{
					OUTTF=false ;
				}
				pos=0;
				if ((pos=Quly_1.find("\t"))!=string::npos )
				{
					OUTTF=false ;
				}


				string::size_type s_Aend=Quly_1.length();

				if (SeqLen!=s_Aend)
				{
					OUTTF=false ;
				}
				if (s_Aend<1)
				{
					OUTTF=false;
				}
				if  (SeqLen<1)
				{
					OUTTF=false;
				}

				string::size_type ix=ID_1.find("PHRED scores");


				if ( ix==string::npos )
				{
					ix=ID_1.find("Solexa scores") ;
					if (ix!=string::npos )
					{
						ix++;
					}
				}
				if ( ix!=string::npos ) // PHRED
				{
					OUTTFALL=false;
					string Tmp=ID_1.substr(ix+13);
					vector<string> inf;
					split(Tmp,inf," \t");
					vector <int > Num ;
					if (inf[0]== "from")
					{
						int Stat=atoi(inf[1].c_str());
						int End=atoi(inf[3].c_str());
						if  (Stat<End)
						{
							for (int jj=Stat ; jj<=End ; jj++)
							{
								Num.push_back(jj);
							}
						}
						else
						{
							for (int jj=Stat ; jj>=End; jj--)
							{
								Num.push_back(jj);
							}
						}
					}
					else
					{
						for (int ii=1 ; ii<inf.size(); ii++)
						{
							if( (ix=inf[ii].find(','))!=string::npos)
							{
								inf[ii].replace(ix,1,"");
								int Stat=atoi(inf[ii].c_str());
								Num.push_back(Stat);
							}
							else
							{
								int Stat=atoi(inf[ii].c_str());
								Num.push_back(Stat);
								break;
							}
						}
					}

					int CCbin=Num.size();
					int Tbintt=0;
					for(ix=0 ; ix<s_Aend ; ix++)
					{
						Tbintt=int(ix % CCbin);
						int cc=Num[Tbintt]+33;
						if (cc>73)
						{
							cc=73;
						}
						else if (cc<33)
						{
							cc=33;
						}
						Quly_1[ix]=cc;
					}
					istringstream isone (ID_1,istringstream::in);
					isone>>ID_1;
					//			ID_1=Tmp;
				}
				for(ix=0 ; ix<s_Aend ; ix++)
				{
					if  (Quly_1[ix]<33 ||  Quly_1[ix]>115)
					{
						OUTTF=false ;
					}
				}
				if (OUTTF)
				{
					OUT<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";
				}
				else
				{
					OUTTFALLV2=false;
				}
			}

			if (ID_1!=ID_2 && (!ID_2.empty()))
			{
				OUTTF=false;
				OUTTFALLV2=false;
			}	

			if (OUTTFALLV2)
			{
				if (OUTTFALL)
				{
					cout<<"\t\t\tVALID FASTQ File"<<endl;
				}
				else
				{
					cout<<"\t\t\tVALID FASTQ File. Note:File is not in normal format, we had change it to the normal Fq "<<(para_A30->InStr3)<<endl;
				}
			}
			else
			{
				cout<<"\t\t\tINVALID FASTQ File"<<endl;
			}


		}









	}


	if ((para_A30->InStr3) ==  "./.tmpValid.fq.gz" )
	{
		string temRm=" rm -rf  ./.tmpValid.fq.gz ";
		std::system(temRm.c_str()) ;
	}



	IN_A.close();
	OUT.close();
	delete para_A30 ;
	return 0 ;
}

///////// swimming in the sky and flying in the sea ////////////
#endif 
