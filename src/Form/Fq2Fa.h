#ifndef Fq2Fa_H_
#define Fq2Fa_H_
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
#include <algorithm>
#include <zlib.h>
#include <stdio.h>
#include "../ALL/kseq.h"
#include "../ALL/DataClass.h"

using namespace std;

int  print_Fq2Fa()
{
	cout <<""
		"\n"
		"\tUsage: Fq2Fa -i <in.fq> -o <out.fa>\n"
		"\n"
		"\t\t-i      <str>   Input Fastq File\n"
		"\t\t-o      <str>   Output Fasta file\n"
		"\t\t                with Hight Quality score \"h\"\n"
		"\n"
		"\t\t-g              No Output with gzip\n"
		"\n"
		"\t\t-h              show this help\n" 
		"\n";
	return 1;
}


int usageFq2Fa(int argc, char **argv , In3str1v * para )
{
	if (argc <=2 ) {print_Fq2Fa();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InFq"  || flag  == "i" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut" || flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para->InStr2=argv[i];
		}
		else if (flag  == "help"  || flag  == "h")
		{
			print_Fq2Fa();return 0;
		}
		else if (flag  == "Nogzip" || flag  == "g" )
		{
			para->TF=false ;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para->InStr1).empty()  ||  (para->InStr2).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	(para->InStr2)=add_Asuffix(para->InStr2) ;
	return 1 ;
}

int Form_Fq2Fa_main(int argc, char *argv[])
{
	In3str1v * para = new In3str1v;
	if( usageFq2Fa(argc, argv, para )==0)
	{
		delete  para ;
		return 1;
	}

    std::ios::sync_with_stdio(false);
	std::cin.tie(0);


	igzstream INFQ ;
  //  char buf[64 * 1024];
	   
    //INFQ.rdbuf()->pubsetbuf(buf, sizeof buf);
	INFQ.open((para->InStr1).c_str(),ifstream::in);
	//igzstream INFQ ((para->InStr1).c_str(),ifstream::in);
	if(!INFQ.good())
	{
		cerr << "open InputFQ  error: "<<(para->InStr1)<<endl;
		return 1;
	}



	if  (para->TF)
	{
	/*///

		gzFile OUT;
		OUT = gzopen ((para->InStr2).c_str(), "wb");
		string seq;
		int IDlength=0;
		char BB[2]; BB[0]='\n'; BB[1]='\0';
		while(!INFQ.eof())
		{
			string  ID ;
			getline(INFQ,ID);
			IDlength=ID.length();
			if (IDlength<=0)  { continue  ; }
			ID[0]='>';
			getline(INFQ,seq);
			gzwrite(OUT,ID.c_str(),IDlength);	
			gzwrite(OUT,BB,1);	
			gzwrite(OUT,seq.c_str(),seq.length());	
			gzwrite(OUT,BB,1);	
			getline(INFQ,ID);
			getline(INFQ,seq);
		}

		gzclose(OUT);
		////*///

		//*///
		   ogzstream OUT((para->InStr2).c_str()) ;
		   if (!OUT.good())
		   {
		   cerr<<"Can't open OutFile "<<(para->InStr2)<<endl;
		   return 1;
		   }
					   
		   string seq;

		   while(!INFQ.eof())
		   {
		   string  ID ;
		   getline(INFQ,ID);
		   if (ID.length()<=0)  { continue  ; }
		   ID[0]='>';
		   getline(INFQ,seq);
		   OUT<<ID<<"\n"<<seq<<"\n";
		   getline(INFQ,ID);
		   getline(INFQ,seq);
		   }
		   OUT.close();
	////	   */
	}
	else
	{
		string OUTFileAA=(para->InStr2).substr(0,(para->InStr2).length()-3);
		ofstream OUT(OUTFileAA.c_str()) ;
		//const unsigned int length = 8192 ;
		//char buffer[length];
	   //OUT.write(buffer, length);
		if (!OUT.good())
		{
			cerr<<"Can't open OutFile "<<OUTFileAA<<endl;
			return 1;
		}
		string seq;
//		seq.reserve(64 * 1024);
		while(!INFQ.eof())
		{
			string  ID ;
			getline(INFQ,ID);
			if (ID.length()<=0)  { continue  ; }
			ID[0]='>';
			getline(INFQ,seq);
			OUT<<ID<<"\n"<<seq<<"\n";
			getline(INFQ,ID);
			getline(INFQ,seq);
		}
		OUT.close();


	}


	INFQ.close();
	delete para ;
	return 0;
}
#endif /// Fq2Fa_H_ 
///////// swimming in the sky and flying in the sea ////////////

