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
#include "../ALL/kseq.h"

//KSEQ_AINIT(gzFile, gzread)

using namespace std;

int  print_usage_gtf2gf()
{
	cout <<""
		"\n"
		"\tUsage: Gtf2Gff  -InGTF <in.gtf>   -OutGFF  <out.gff>\n"
		"\n"
		"\t\t-InGTF     <str>   Input the GTF file\n"
		"\t\t-OutGFF    <str>   Output the GFF file\n"
		"\n" 
		"\t\t-change            change the exon 2 CDS\n"
		"\t\t-help              show this help\n" 
		"\n"
		"\t\tFIXME: stop_codon add to last CDS (also change to pipe)\n"
		"\n";
	return 1;
}

int parse_cmd_gtf2gf(int argc, char **argv , In3str1v * para_A1)
{
	if (argc <=2 ) {print_usage_gtf2gf();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InGTF" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_A1->InStr1=argv[i];
		}
		else if (flag  ==  "OutGFF")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_A1->InStr2=argv[i];
		}
		else if (flag  == "chang" ||  flag  == "change" )
		{
			para_A1->TF=false ;
		}
		else if (flag  == "help")
		{
			print_usage_gtf2gf();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para_A1->InStr1).empty() ||  (para_A1->InStr2).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	(para_A1->InStr2)= add_Asuffix ( (para_A1->InStr2) ) ;

	return 1 ;
}


int GTF2GFF_main(int argc, char *argv[])
{
	In3str1v * para_A1 = new In3str1v;
	if( parse_cmd_gtf2gf(argc, argv, para_A1)==0)
	{
		delete  para_A1 ;
		return 1;
	}
	ogzstream OUT((para_A1->InStr2).c_str()) ;
	if (!OUT.good())
	{
		cerr<<"Can't open OutFile "<<(para_A1->InStr2)<<endl;
		return 1;
	}
	igzstream IN ((para_A1->InStr1).c_str(),ifstream::in);
	if(!IN.good())
	{
		cerr << "open Input  error: "<<(para_A1->InStr1)<<endl;
		return 1;
	}

	OUT<<"##gff-version 3"<<endl;
	while(!IN.eof())
	{
		string  line ;
		getline(IN,line);
		if (line.length()<=0)  { continue  ; }
		if (line[0] == '#')  { OUT<<line<<endl; continue  ; }
		vector<string> inf;
		split(line,inf,"\t");
		if   (!(para_A1->TF)) 
		{
			if ( (inf[2] == "stop_codon" || inf[2] == "start_codon")  ) { continue  ; }
			else if (inf[2] == "transcript" ) { inf[2]="mRNA" ;}
			else if (inf[2] == "5UTR" ) { inf[2]="five_prime_UTR" ;}
			else if (inf[2] == "3UTR" ) { inf[2]="three_prime_UTR" ;}
			else { }
			vector<string> BB;
			split(inf[8],BB,";");
			int Len=BB.size();
			string Info="";
			for(int  i=0 ; i<Len  ; i++ )
			{
				vector<string> CC;
				vector<string> DD;
				split(BB[i],CC," ");
				int LenPP=CC.size();
				for(int  j=0 ; j<LenPP ; j++)
				{
					if (!(CC[j].empty()))
					{
						DD.push_back(CC[j]);
					}
				}
				if (DD[0] == "transcript_id")
				{
					DD[0]="Parent";
					if ( inf[2]=="mRNA")
					{
						DD[0]="ID";
					}
				}
				else if (DD[0] == "exon_number"  &&  inf[2] == "exon" )
				{
					inf[2] ="CDS";
				}
				DD[1]=replace_all(DD[1],"\"","");
				Info=Info+DD[0]+"="+DD[1]+";";
			}
			for (int i=0 ; i<8 ; i++)
			{
				OUT<<inf[i]<<"\t";
			}
			OUT<<Info<<endl;
		}
		else
		{

			vector<string> BB;
			split(inf[8],BB,";");
			int Len=BB.size();
			string Info="";
			for(int  i=0 ; i<Len  ; i++ )
			{
				vector<string> CC;
				vector<string> DD;
				split(BB[i],CC," ");
				int LenPP=CC.size();
				for(int  j=0 ; j<LenPP ; j++)
				{
					if (!(CC[j].empty()))
					{
						DD.push_back(CC[j]);
					}
				}
				DD[1]=replace_all(DD[1],"\"","");
				Info=Info+DD[0]+"="+DD[1]+";";
			}
			for (int i=0 ; i<8 ; i++)
			{
				OUT<<inf[i]<<"\t";
			}
			OUT<<Info<<endl;
		}
	}

	IN.close();
	OUT.close();
	delete para_A1 ;
	return 0;
}

///////// swimming in the sky and flying in the sea ////////////
