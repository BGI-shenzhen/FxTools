#ifndef RefSubSeq_H_
#define RefSubSeq_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "../ALL/comm.h"
#include <cstdlib>
#include "../ALL/gzstream.h"
#include <zlib.h>
#include <stdio.h>
#include "../ALL/kseq.h"
#include "RefRegenerate.h"

using namespace std;
typedef long long  llong ;


int  print_AusageFA01()
{
	cout <<""
		"\n"
		"\tUsage: grep  -i <in.fa> -o <out.fa> -s chr:start:end\n"
		"\n"
		"\t\t-i       <str>   input FASTA\n"
		"\t\t-o       <str>   output file, default [STDOUT]\n"
		"\t\t-s       <str>   single region to extract, format as [chr:start:end]\n"
		"\t\t-m       <str>   file containing multiple regions to extract, format as [chr start end]\n"
		"\n"
		"\t\t-r               reverse of the output sequences\n"
		"\t\t-c               complement of the output sequences\n"
		"\n"
		"\t\t-h               show more details for help\n" 
		"\n";
	return 1;
}

void More_HelpFA_6 ()
{
	cout<<"\n\
\n\
\t\t1.  grep -i <in.fa> -o AAA -s seqX:x1:x2\n\
\t\tthis will extract the read which locates on seqX from base site x1 to x2 in FASTA and output the result in a compressed file named AAA.gz. \n\
\n\
\t\t2.  grep -i <in.fa> -o AAA -m multiple_region_file\n\
\t\tif user has more than one regions to extract from FASTA, please list all the regions in a file. In this example, we list all the regions in a file named multiple_region_file and it is shown as:\n\
\t\tseq1 x1 x2\n\
\t\tseq2 x3 x4\n\
\t\tBy doing this, reads which locating on seq1 from base site x1 to x2 and reads locating on seq2 from base site x3 to x4 in FASTA would be output to a compressed file named AAA.gz.\n\
\n\
\t\t3.  if user would like to get the reverse complement of the output sequences, just simply add the -r and -c option.\n\
\n\
\n";
}

int parse_AcmdFA01(int argc, char **argv , ParaClass * paraFA01 )
{
	if (argc <=1 ) {print_AusageFA01();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InPut" || flag  == "i")
		{
			if(i + 1 == argc) { LogLackArg(flag); return 0;}
			i++;
			paraFA01->InPut1=argv[i];
		}
		else if (flag  ==  "OutPut"  ||  flag  == "o")
		{
			if(i + 1 == argc) { LogLackArg(flag); return 0;}
			i++;
			paraFA01->OutPut1=argv[i];
		}
		else if (flag  ==  "ORegion" || flag  == "s" )
		{
			if(i + 1 == argc) { LogLackArg(flag); return 0;}
			i++;
			paraFA01->InStr1=argv[i];
		}
		else if (flag  ==  "MRegion"  ||  flag  == "m")
		{
			if(i + 1 == argc) { LogLackArg(flag); return 0;}
			i++;
			paraFA01->InStr2=argv[i];
		}
		else if (flag  ==  "Reverse"  ||  flag  == "r")
		{
			paraFA01->InInt1 = (paraFA01->InInt1) | 0x1;
		}
		else if (flag  ==  "Complement" || flag  == "c")
		{
			paraFA01->InInt1 = (paraFA01->InInt1) | 0x2;
		}
		else if (flag  == "help" || flag  == "h" )
		{
			More_HelpFA_6();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((paraFA01->InPut1).empty() || ( (paraFA01->InStr1).empty() &&   (paraFA01->InStr2).empty()   )  )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	if (!(paraFA01->OutPut1).empty())
	{
		(paraFA01->OutPut1)=add_Asuffix(paraFA01->OutPut1) ;
	}

	return 1 ;
}

int FA_SubSeq_main(int argc, char *argv[])
	//int main(int argc, char *argv[])
{
	ParaClass * paraFA01 = new ParaClass;
	if( parse_AcmdFA01(argc, argv, paraFA01 )==0)
	{
		delete  paraFA01 ;
		return 1 ;
	}

	map <string,vector <Region> > RegionList ;

	if (!(paraFA01->InStr1).empty())
	{
		vector<string> Temp;
		split((paraFA01->InStr1),Temp,":");
		Region tmp ;
		tmp.Start=atol(Temp[1].c_str());
		tmp.End=atol(Temp[2].c_str());
		vector <Region> R (1,tmp);
		RegionList.insert(map  <string,vector <Region> > :: value_type(Temp[0],R));
	}

	if (!(paraFA01->InStr2).empty())
	{
		igzstream  MRG ((paraFA01->InStr2).c_str(),ifstream::in);
		if(!MRG.good())
		{
			cerr << "open InputFile error: "<<(paraFA01->InStr2)<<endl;
			delete  paraFA01 ; return 1;
		}
		while(!MRG.eof())
		{
			string  line ;
			getline(MRG,line);
			if (line.length()<=0)  { continue  ; }
			//            istringstream isone (line,istringstream::in);
			string chr;
			Region tmp ;
			vector<string> Temp;
			split(line,Temp,": \t");
			tmp.Start=atol(Temp[1].c_str());
			tmp.End=atol(Temp[2].c_str());
			chr=Temp[0];
			//          isone>>chr>>tmp.Start>>tmp.End;
			map  <string,vector <Region> >  :: iterator it=RegionList.find(chr);
			if (it == RegionList.end() )
			{
				vector <Region> R (1,tmp);
				RegionList.insert(map  <string,vector <Region> > :: value_type(chr,R));
			}
			else
			{
				(it->second).push_back(tmp);
			}
		}
		MRG.close();
	}

	int linecut= FaCutLine ((paraFA01->InPut1));

	gzFile fp;
	kseq_t *seq;
	int l;

	fp = gzopen((paraFA01->InPut1).c_str(), "r");
	seq = kseq_init(fp);    





	map <char ,char >  Complement;
	Complement['A']='T'; Complement['C']='G';  Complement['T']='A'; Complement['G']='C'; Complement['N']='N';
	Complement['a']='t'; Complement['c']='g';  Complement['t']='a'; Complement['g']='c'; Complement['n']='n';
	Complement['M']='K'; Complement['R']='Y';  Complement['W']='W'; Complement['S']='S'; Complement['Y']='R';
	Complement['K']='M';
	// $$seq_p=~tr/AGCTagctMRWSYK/TCGAtcgaKYWSRM/;

	if (!(paraFA01->OutPut1).empty())
	{
		ogzstream OUT ((paraFA01->OutPut1).c_str());
		if(!OUT.good())
		{
			cerr << "open OUT File error: "<<(paraFA01->OutPut1)<<endl;
			delete  paraFA01 ; return  1;
		}

		while ((l = kseq_read(seq)) >= 0)
		{
			string chr=seq->name.s ;
			map <string,vector <Region> >  :: iterator it=RegionList.find(chr);
			if (it != RegionList.end())
			{
				int Vec_Size=(it->second).size();
				string Ref= seq->seq.s ;

				for(int ii=0 ; ii< Vec_Size ; ii++ )
				{
					llong START=((it->second)[ii]).Start-1;
					if (START<0)
					{
						cerr<<"Region OUT Range at : "<<chr<<"\t"<<(((it->second)[ii]).Start)<<"\t"<<((it->second)[ii]).End<<endl;
						continue ;
					}
					int Leng=((it->second)[ii]).End-START;
					string BB=Ref.substr(START,Leng) ;
					if ( (paraFA01->InInt1) & 0x1  )
					{
						reverse(BB.begin(), BB.end());
					}
					if  ( (paraFA01->InInt1) & 0x2  )
					{
						int leng=BB.size();
						for (int i=0 ; i<leng ; i++)
						{
							BB[i]=Complement[BB[i]];
						}
					}
					string ID=chr+"_"+Int2Str(((it->second)[ii]).Start)+"_"+Int2Str(((it->second)[ii]).End);
					Display(BB , ID ,  OUT ,linecut  );
				}
			}
		}
		OUT.close();
	}
	else
	{

		while ((l = kseq_read(seq)) >= 0)
		{
			string chr=seq->name.s ;
			map <string,vector <Region> >  :: iterator it=RegionList.find(chr);
			if (it != RegionList.end())
			{
				int Vec_Size=(it->second).size();
				string Ref= seq->seq.s ;

				for(int ii=0 ; ii< Vec_Size ; ii++ )
				{
					llong START=((it->second)[ii]).Start-1;
					if (START<0)
					{
						cerr<<"Region OUT Range at : "<<chr<<"\t"<<(((it->second)[ii]).Start)<<"\t"<<((it->second)[ii]).End<<endl;
						continue ;
					}
					int Leng=((it->second)[ii]).End-START;
					string BB=Ref.substr(START,Leng) ;
					if ( (paraFA01->InInt1) & 0x1  )
					{
						reverse(BB.begin(), BB.end());
					}
					if  ( (paraFA01->InInt1) & 0x2  )
					{
						int leng=BB.size();
						for (int i=0 ; i<leng ; i++)
						{
							BB[i]=Complement[BB[i]];
						}
					}
					string ID=chr+"_"+Int2Str(((it->second)[ii]).Start)+"_"+Int2Str(((it->second)[ii]).End);
					Display(BB , ID , linecut  );
				}
			}
		}



	}



	kseq_destroy(seq);
	gzclose(fp);
	delete paraFA01 ;
	return 0;

}

#endif // RefSubSeq_H_ //
///////// swimming in the sky and flying in the sea ////////////

