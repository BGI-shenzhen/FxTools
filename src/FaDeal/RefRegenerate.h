#ifndef RefRegenerate_H_
#define RefRegenerate_H_

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
#include "RefStat.h"
#include "RefFilter.h"

using namespace std;
typedef long long  llong ;

int  usage_FA05()
{
	cout <<""
		"\n"
		"\tUsage: regenerate -i  <in.fa> -o <out.fa> \n"
		"\n"
		"\t\t-i     <str>    Input fa for regenerate\n"
		"\t\t-o     <str>    Output file fa of new Ref\n"
		"\n"
		"\t\t-c     <str>    Original chr name for NoChang[Chr]\n"
		"\t\t-n     <str>    New Chr for merge Seq[NewChr]\n"
		"\t\t-b     <int>    Insert Num of N between two Scaf [150]\n"
		"\t\t-s     <int>    the Num for merge new Seq [20]\n"
		"\t\t-l     <int>    Len of each new merge Seq [NA]\n"
		"\t\t-g              OutPut file with No gz [NA]\n"
		"\n"
		"\t\t-h              show this help\n"
		"\n";
	return 1;
}


int parse_Acmd_FA05(int argc, char **argv , ParaClass * para_FA05 )
{
	if (argc <=2 ) {usage_FA05();return 0;}

	int tmp=0 ;
	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InPut" ||  flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_FA05->InPut1=argv[i];
		}
		else if (flag  ==  "OutPut" ||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FA05->OutPut1=argv[i];
		}       
		else if (flag  ==  "InsertN" || flag  == "b")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FA05->InInt1=atoi(argv[i]);
		}
		else if (flag  == "NumSeq" || flag  == "s")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FA05->InInt2=atoi(argv[i]);
			tmp++;
		}
		else if (flag  == "Length" || flag  == "l")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FA05->Inllong1=atol(argv[i]);
			tmp++;
		}
		else if (flag  == "OriChr" ||  flag  == "c")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FA05->InStr1=argv[i];
		}
		else if (flag  == "NewChr" ||  flag  == "n")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FA05->InStr2=argv[i];
		}
		else if (flag  == "NoOutgz" || flag  == "g")
		{
			para_FA05->TF=true;
		}
		else if (flag  == "help"||  flag  == "h")
		{
			usage_FA05();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para_FA05->InPut1).empty() || (para_FA05->OutPut1).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	if (tmp>1)
	{
		cerr<< "para -l & -s  should not together"<<endl ;
		return 0;
	}
	return 1 ;
}




//int main(int argc, char *argv[])
int FA_Regenerate_main(int argc, char *argv[])
{
	ParaClass * para_FA05 = new ParaClass ;
	para_FA05->InInt1=150;
	para_FA05->TF=false;
	para_FA05->InInt2=20;
	para_FA05->InStr1="Chr";
	para_FA05->InStr2="NewChr";
	if( parse_Acmd_FA05(argc, argv, para_FA05 )==0)
	{
		delete  para_FA05 ;
		return 1;    
	}

	int linecut= FaCutLine ((para_FA05->InPut1));

	string mergelist=(para_FA05->OutPut1)+".merlist";
	(para_FA05->OutPut1)=add_Asuffix((para_FA05->OutPut1));
	string chrlist=(para_FA05->OutPut1)+".chrlist";
	ofstream  MER (mergelist.c_str());
	ogzstream  OUT ((para_FA05->OutPut1).c_str());
	if(!OUT.good())
	{
		cerr << "open OUT File error: "<<(para_FA05->OutPut1)<<endl;
		return 1;
	}
	if (MER.fail())
	{
		cerr << "open OUT File error: "<<mergelist<<endl;
		return 1;
	}

	MER<<"##OldChr\tStart\tEnd\tNewChr\tStart\tEnd"<<endl;

	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen((para_FA05->InPut1).c_str(), "r");
	seq = kseq_init(fp);
	map <string ,llong > ChrLen ;
	map <string ,string> ChrSeq ;
	map <llong ,vector <string> > SortLeng ;
	llong SumLeng=0;

	while ((l = kseq_read(seq)) >= 0) 
	{
		string RefSeq=(seq->seq.s);
		string chr=(seq->name.s);
		llong seq_length = (seq->seq.l);
		if (chr.find(para_FA05->InStr1) == string::npos)
		{
			ChrLen[chr]=seq_length;
			ChrSeq[chr]=RefSeq;
			SumLeng+=(seq_length+(para_FA05->InInt1));
			map <llong ,vector <string> > ::iterator it=SortLeng.find(seq_length) ;
			if (it==SortLeng.end())
			{
				vector <string> tmp;
				tmp.push_back(chr);
				SortLeng.insert(map <llong,vector <string> > ::value_type(seq_length,tmp));
			}
			else
			{
				(it->second).push_back(chr);
			}
		}
		else
		{
			if (seq->comment.l)
			{
				string tmp=seq->comment.s;
				chr=chr+"\t"+tmp;
			}
			Display( RefSeq , chr  , OUT ,linecut);
			MER<<chr<<"\t1\t"<<seq_length<<"\t"<<chr<<"\t1\t"<<seq_length<<"\n";
		}
	}

	kseq_destroy(seq);
	gzclose(fp);
	vector <string>  SortChr ;
	map <llong ,vector <string> > ::iterator it=SortLeng.begin();
	for ( ; it!=SortLeng.end();it++)
	{
		int Size=(it->second).size();
		for (int ii=0; ii<Size ; ii++)
		{
			SortChr.push_back((it->second)[ii]);
		}
	}
	SortLeng.clear();

	llong Each_NewSeq_length=(SumLeng/(para_FA05->InInt2));
	if ((para_FA05->Inllong1)!=0)
	{
		Each_NewSeq_length=para_FA05->Inllong1;
	}

	string N; 
	for (int jj=0; jj<(para_FA05->InInt1) ;  jj++ )
	{
		N+="N";
	}


	int size_length=SortChr.size();
	llong Tmp=0;
	int New_chr_ID=1;
	string  ID=(para_FA05->InStr2)+Int2Str(New_chr_ID);

	for (int i=size_length-1; i>-1; i--)
	{
		if (Tmp==0)
		{
			Tmp=ChrLen[SortChr[i]];
			ID=(para_FA05->InStr2)+Int2Str(New_chr_ID);
			New_chr_ID++;
			OUT<<">"<<ID<<"\n" ;
			Display( ChrSeq[SortChr[i]] , OUT , linecut);
			MER<<SortChr[i]<<"\t1\t"<<Tmp<<"\t"<<ID<<"\t1\t"<<Tmp<<"\n";
			if (Tmp>Each_NewSeq_length)
			{
				Tmp=0;
			}			
		}
		else
		{
			Display( N  , OUT , linecut );
			Display( ChrSeq[SortChr[i]] , OUT , linecut);
			llong tmpLength=ChrLen[SortChr[i]];
			Tmp+=(para_FA05->InInt1);
			MER<<SortChr[i]<<"\t1\t"<<tmpLength<<"\t"<<ID<<"\t"<<(Tmp+1)<<"\t";
			Tmp+=tmpLength;
			MER<<Tmp<<"\n";
			if (Tmp>Each_NewSeq_length)
			{
				Tmp=0;
			}
		}
	}

	OUT.close();

	char * A = const_cast<char*>((para_FA05->OutPut1).c_str());
	string DD=(para_FA05->OutPut1)+".filter.gz";
	char * C = const_cast<char*>((DD).c_str());
	char * TmpFF[9]={ (char *)"filter", (char *)"-InPut", A , (char *)"-OutPut" , C , (char *)"-MinLen" ,(char *)"100", (char *)"-NRatio", (char *)"0.88" };
	FA_Filter_main( 9 , TmpFF );
	string  MV="mv  "+DD+" "+(para_FA05->OutPut1) ;
	std::system(MV.c_str()) ;


	char * B = const_cast<char*>((chrlist).c_str());
	char * ssTmp[5]={(char *)"stat", (char *)"-InPut", A ,(char *)"-OutPut" , B } ;
	FA_stat_main(5 , ssTmp ) ;

	if  (para_FA05->TF)
	{
		string Gzip="gzip -d  "+(para_FA05->OutPut1);
		std::system(Gzip.c_str());
	}

	delete para_FA05 ;
	return 0;
}

#endif
///////// swimming in the sky and flying in the sea ////////////
