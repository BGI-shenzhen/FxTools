#ifndef VCF2Genotype_H_
#define VCF2Genotype_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "../ALL/comm.h"
#include "../ALL/DataClass.h"
#include <cstdlib>
#include <gzstream.h>
#include <zlib.h>
#include <stdio.h>

using namespace std;

int  VCF2GenotypeSNP_help()
{
	cout <<""
		"\n"
		"\tUsage: VCF2Genotype -InPut  <in.vcf>  -OutPut  <out.add_ref>\n"
		"\n"
		"\t\t-InPut   <str>   Input GATK VCF genotype File\n"
		"\t\t-OutPut  <str>   OutPut Genotype file\n"
		"\n"
		"\t\t-SubPop <str>    SubGroup SampleList of VCFFile [ALLsample]\n"
		"\t\t-NoRef           Do not pring OUT Ref base[with ref base]\n"
		"\t\t-WithHeader      OutPut file add the header infomation\n"
		"\t\t-Warn            warning the Indel site for no chang\n"
		"\t\t-help            show this help\n"
		"\n";
	return 1;
}

int VCF2Genotype_help01(int argc, char **argv , In3str1v * paraFA04 )
{
	if (argc <=2 ) {VCF2GenotypeSNP_help();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InPut")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
		}
		else if (flag  ==  "SubPop")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr3=argv[i];
		}
		else if (flag == "Warn")
		{
			paraFA04->TF=false;
		}
		else if (flag == "WithHeader")
		{
			paraFA04->InInt=1;
		}
		else if (flag == "NoRef")
		{
			paraFA04->TF2=false;
		}
		else if (flag == "help")
		{
			VCF2GenotypeSNP_help();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((paraFA04->InStr1).empty() || (paraFA04->InStr2).empty()  )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	(paraFA04->InStr2)=add_Asuffix(paraFA04->InStr2);
	return 1 ;
}


int VCF2GenotypeSNP_main(int argc, char *argv[])
//	int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	if ((VCF2Genotype_help01(argc, argv, paraFA04)==0))
	{
		delete paraFA04 ;
		return 1;
	}

	map <string ,char > SNP_Allele ;
	SNP_Allele["AC"]='M'; SNP_Allele["CA"]='M'; SNP_Allele["GT"]='K'; SNP_Allele["TG"]='K';
	SNP_Allele["CT"]='Y'; SNP_Allele["TC"]='Y'; SNP_Allele["AG"]='R'; SNP_Allele["GA"]='R';
	SNP_Allele["AT"]='W'; SNP_Allele["TA"]='W'; SNP_Allele["CG"]='S'; SNP_Allele["GC"]='S';
	SNP_Allele["AA"]='A'; SNP_Allele["TT"]='T'; SNP_Allele["CC"]='C'; SNP_Allele["GG"]='G';
	SNP_Allele["AN"]='A'; SNP_Allele["NA"]='A';
	SNP_Allele["CN"]='C'; SNP_Allele["NC"]='C';
	SNP_Allele["TN"]='T'; SNP_Allele["NT"]='T';
	SNP_Allele["GN"]='G'; SNP_Allele["NG"]='T';


	ogzstream OUT ((paraFA04->InStr2).c_str());

	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<(paraFA04->InStr2)<<endl;
		delete  paraFA04 ; return 1;
	}


	igzstream INVCFFile ((paraFA04->InStr1).c_str(),ifstream::in);
	if (INVCFFile.fail())
	{
		cerr << "open Genotype File IN File error: "<<(paraFA04->InStr1)<<endl;
		delete  paraFA04 ; return 1;
	}

	vector <int> SampleSite;
	int NumberSubGroup=0;

	while(!INVCFFile.eof())
	{
		string  line ;
		getline(INVCFFile,line);
		if (line.length()<=0 )  { continue  ; }
		else if (line[0] == '#' && line[1] == '#' )  { continue  ; }
		else if( line[0] == '#' && line[1] != '#')
		{
			vector<string> Vsample ;
			split(line,Vsample," \t");
			if  ( Vsample[0]  != "#CHROM")
			{
				continue  ;
			}
			int ASize=Vsample.size();

			if  ((paraFA04->InStr3).empty())
			{
				for (int jj=9 ; jj< ASize ;jj++)
				{
					SampleSite.push_back(jj);
				}
				NumberSubGroup=SampleSite.size();
			}
			else
			{

				igzstream SampleList ((paraFA04->InStr3).c_str(),ifstream::in);
				if (SampleList.fail())
				{
					cerr << "open Sub Group IN File error: "<<(paraFA04->InStr3)<<endl;
					return 1;
				}

				map <string ,int >  SubVetor;
				map <string ,int >  :: iterator it;

				while(!SampleList.eof())
				{
					string  line ;
					getline(SampleList,line);
					if (line.length()<=0 || line[0] == '#' )  { continue  ; }
					vector<string> inf ;
					split(line,inf," \t");
					int A=inf.size();
					for(int ii=0 ; ii<A ; ii++)
					{
						it=SubVetor.find(inf[ii]);
						if (it==SubVetor.end())
						{
							SubVetor.insert(map <string ,int> ::value_type(inf[ii],1));
						}
						else
						{
							(it->second)++;
						}
					}
				}
				SampleList.close();

				for(int ii=9 ; ii< ASize ; ii++)
				{
					it=SubVetor.find(Vsample[ii]);
					if (it!=SubVetor.end())
					{
						SampleSite.push_back(ii);
					}
				}

				NumberSubGroup=SampleSite.size();
				cout<<"the Number of subPop samples[found in VCF] is "<<NumberSubGroup<<endl;

			}

			if ((paraFA04->InInt)==1)
			{
				if (paraFA04->TF2)
				{
					OUT<<Vsample[0]<<"\t"<<Vsample[1]<<"\t"<<Vsample[3];
					for (int kk=0 ; kk< NumberSubGroup ; kk++)
					{
						OUT<<" "<<Vsample[SampleSite[kk]];
					}
					OUT<<endl;
				}
				else
				{
					OUT<<Vsample[0]<<"\t"<<Vsample[1]<<"\t"<<Vsample[SampleSite[0]];

					for (int kk=1 ; kk< NumberSubGroup ; kk++)
					{
						OUT<<" "<<Vsample[SampleSite[kk]];
					}
					OUT<<endl;

				}


			}
			break ;
		}
		else if ( line[0] != '#')
		{
			cerr<<"wrong Line : "<<line<<endl;
			cerr<<"Genotype Header same thing wrong, can find sample info before site info"<<endl;
			cerr<<"Genotype Header sample info Flag : [  #CHROM  ] "<<endl;
			return 1;
			break;
		}
	}






	if (paraFA04->TF)
	{
		while(!INVCFFile.eof())
		{
			string  line ;
			getline(INVCFFile,line);
			if (line.length()<=0 || line[0] == '#' )  { continue  ; }
			string GeneID  ;
			vector<string> inf ;
			split(line,inf," \t");
			int  Base_len=inf[3].length();
			vector<string> Alt ;
			split(inf[4],Alt,",");
			map <int,string> Num2Base ;
			Num2Base[0]=inf[3];
			for (int ii=0 ; ii<Alt.size() ;ii++)
			{
				if (Alt[ii].length()>Base_len)
				{
					Base_len=Alt[ii].length();
				}
				int j=ii+1;
				Num2Base[j]=Alt[ii];
			}
			if (Base_len>1)
			{
				continue ;
			}

			map <string,char> Num2Genotype ;
			map <int,string> ::iterator  key1 ;
			map <int,string> ::iterator  key2 ;
			for (key1=Num2Base.begin(); key1!=Num2Base.end() ; key1++)
			{

				string keyA=Int2Str(key1->first);
				char  ValueS=(key1->second)[0];
				Num2Genotype[keyA]=ValueS;
				for (key2=Num2Base.begin(); key2!=Num2Base.end() ; key2++)
				{
					string doub_A=(key1->second)+(key2->second);
					char Value=SNP_Allele[doub_A];
					string keyB=Int2Str(key2->first);
					string Final_key=keyA+"/"+keyB;
					Num2Genotype[Final_key]=Value;
					Final_key=keyB+"/"+keyA;
					Num2Genotype[Final_key]=Value;
					Final_key=keyA+"|"+keyB;
					Num2Genotype[Final_key]=Value;
					Final_key=keyB+"|"+keyA;
					Num2Genotype[Final_key]=Value;
				}
			}
			string temp="./." ;			Num2Genotype[temp]='-';
			temp=".|." ;			Num2Genotype[temp]='-';
			temp="." ;			Num2Genotype[temp]='-';

			int A=inf.size();
			if (paraFA04->TF2)
			{

				OUT<<inf[0]<<"\t"<<inf[1]<<"\t"<<inf[3];
				for (int kk=0 ; kk< NumberSubGroup ; kk++)
				{
					vector<string> Btmp  ; 
					split(inf[SampleSite[kk]], Btmp,":");
					char ctmp='N';
					ctmp=Num2Genotype[Btmp[0]];
					OUT<<" "<<ctmp;
				}
				OUT<<endl;
			}
			else
			{
				vector<string> BtmpA  ;
				split(inf[SampleSite[0]], BtmpA,":");
				char ctmpA='N';
				ctmpA=Num2Genotype[BtmpA[0]];

				OUT<<inf[0]<<"\t"<<inf[1]<<"\t"<<ctmpA;
				for (int kk=1 ; kk< NumberSubGroup ; kk++)
				{
					vector<string> Btmp  ; 
					split(inf[SampleSite[kk]], Btmp,":");
					char ctmp='N';
					ctmp=Num2Genotype[Btmp[0]];
					OUT<<" "<<ctmp;
				}
				OUT<<endl;
			}
		}
	}
	else
	{
		while(!INVCFFile.eof())
		{
			string  line ;
			getline(INVCFFile,line);
			if (line.length()<=0 || line[0] == '#' )  { continue  ; }
			string GeneID  ;
			vector<string> inf ;
			split(line,inf," \t");
			int  Base_len=inf[3].length();
			vector<string> Alt ;
			split(inf[4],Alt,",");
			map <int,string> Num2Base ;
			Num2Base[0]=inf[3];
			for (int ii=0 ; ii<Alt.size() ;ii++)
			{
				if (Alt[ii].length()>Base_len)
				{
					Base_len=Alt[ii].length();
				}
				int j=ii+1;
				Num2Base[j]=Alt[ii]; 
			}
			if (Base_len>1)
			{
				cerr<<"warning Indel site\t"<<inf[0]<<"\t"<<inf[1]<<endl;
				continue ;
			}

			map <string,char> Num2Genotype ;
			map <int,string> ::iterator  key1 ;
			map <int,string> ::iterator  key2 ;
			for (key1=Num2Base.begin(); key1!=Num2Base.end() ; key1++)
			{
				for (key2=Num2Base.begin(); key2!=Num2Base.end() ; key2++)
				{
					string doub_A=(key1->second)+(key2->second);
					char Value=SNP_Allele[doub_A];
					string keyA=Int2Str(key1->first);
					string keyB=Int2Str(key2->first);
					string Final_key=keyA+"/"+keyB;
					Num2Genotype[Final_key]=Value;
					Final_key=keyB+"/"+keyA;
					Num2Genotype[Final_key]=Value;
					Final_key=keyA+"|"+keyB;
					Num2Genotype[Final_key]=Value;
					Final_key=keyB+"|"+keyA;
					Num2Genotype[Final_key]=Value;
				}
			}
			string temp="./." ;			Num2Genotype[temp]='-';
			temp=".|." ;			Num2Genotype[temp]='-';

			int A=inf.size();

			if (paraFA04->TF2)
			{

				OUT<<inf[0]<<"\t"<<inf[1]<<"\t"<<inf[3];
				for (int kk=0 ; kk< NumberSubGroup ; kk++)
				{
					vector<string> Btmp  ; 
					split(inf[SampleSite[kk]], Btmp,":");
					char ctmp='N';
					ctmp=Num2Genotype[Btmp[0]];
					OUT<<" "<<ctmp;
				}
				OUT<<endl;
			}
			else
			{
				vector<string> BtmpA  ;
				split(inf[SampleSite[0]], BtmpA,":");
				char ctmpA='N';
				ctmpA=Num2Genotype[BtmpA[0]];

				OUT<<inf[0]<<"\t"<<inf[1]<<"\t"<<ctmpA;
				for (int kk=1 ; kk< NumberSubGroup ; kk++)
				{
					vector<string> Btmp  ; 
					split(inf[SampleSite[kk]], Btmp,":");
					char ctmp='N';
					ctmp=Num2Genotype[Btmp[0]];
					OUT<<" "<<ctmp;
				}
				OUT<<endl;
			}


		}
	}


	INVCFFile.close();
	OUT.close();
	delete paraFA04 ;
	return 0;
}
#endif // VCF2Genotype_H_  //
///////// swimming in the sky and flying in the sea ////////////

