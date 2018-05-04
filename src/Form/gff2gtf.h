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

int  print_usage_gff2gtf()
{
	cout <<""
		"\n"
		"\tUsage: Gff2Gtf  -InGFF  <in.gff>  -OutGTF <out.gtf>\n"
		"\n"
		"\t\t-InGFF     <str>   Input the GFF file\n"
		"\t\t-OutGTF    <str>   Output the GTF file\n"
		"\n"
		"\t\t-help              show this help\n" 
		"\n";
	return 1;
}

int parse_cmd_gff2gtf(int argc, char **argv , In3str1v * para_A1)
{
	if (argc <=2 ) {print_usage_gff2gtf();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InGFF" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_A1->InStr1=argv[i];
		}
		///*
		else if (flag  ==  "OutGTF")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_A1->InStr2=argv[i];
		}
		///*////
		else if (flag  == "help")
		{
			print_usage_gff2gtf();return 0;
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


int GFF2GTF_main(int argc, char *argv[])
{
	In3str1v * para_A1 = new In3str1v;
	if( parse_cmd_gff2gtf(argc, argv, para_A1)==0)
	{
		delete  para_A1 ;
		return 0;
	}
	ogzstream OUT((para_A1->InStr2).c_str()) ;
	if (!OUT.good())
	{
		cerr<<"Can't open OutFile "<<(para_A1->InStr2)<<endl;
		return 0;
	}
	igzstream IN ((para_A1->InStr1).c_str(),ifstream::in);
	if(!IN.good())
	{
		cerr << "open Input  error: "<<(para_A1->InStr1)<<endl;
		return 0;
	}
	int count=0 ;
	string gtfid="";
	string geneid="";

	while(!IN.eof())
	{
		string  line ;
		getline(IN,line);
		if (line.length()<=0)  { continue  ; }
		if (line[0] == '#')  { OUT<<line<<endl; continue  ; }

		vector<string> inf;
		split(line,inf," \t");
		if (inf.size()<3){ continue  ; }
		if (inf[2] == "mRNA" )
		{
			count=1; 
			vector<string> BB;
			split(inf[8],BB,";");
			vector<string> CC;
			split(BB[0],CC,"=");
			geneid=CC[1];
			gtfid=CC[1]+".1";
			OUT<<inf[0]<<"\tCufflinks\ttranscript\t"<<inf[3]<<"\t"<<inf[4]<<"\t0\t"<<inf[6]<<"\t"<<inf[7]<<"\tgene_id \""<<geneid<<"\"; transcript_id \" "<<gtfid<<"\";"<<endl;
		}
		else if  (inf[2] == "CDS" )
		{
			OUT<<inf[0]<<"\tCufflinks\texon\t"<<inf[3]<<"\t"<<inf[4]<<"\t0\t"<<inf[6]<<"\t"<<inf[7]<<"\tgene_id \""<<geneid<<"\"; transcript_id \" "<<gtfid<<"\"; exon_number \""<<count<<"\";"<<endl;
			count++;
		}
	}

	IN.close();
	OUT.close();
	delete para_A1 ;
	return 0;
}

///////// swimming in the sky and flying in the sea ////////////
