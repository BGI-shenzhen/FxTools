
#ifndef RefN50Stat_H_
#define RefN50Stat_H_

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
#include <math.h>

using namespace std;
typedef long long  llong ;

void FaStatN50 ( string  InStr2 ,  string  InStr1 , int  InInt1 )
{

	std::ostream  *OUT = & cout;
	std::ofstream fout;
	if (!(InStr2.empty()))
	{
		fout.open(InStr2.c_str());
		OUT = &fout;
	}

	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen(InStr1.c_str(),"r");
	seq = kseq_init(fp);
	map <llong ,int > ChrLen ;
	map <int,int > BinDis ;
	llong  BaseCount = 0;
	int ReadCount=0;

	map <llong ,int > ChrLen_Conti ;
	map <int,int > BinDis_Conti ;
	llong  BaseCount_Conti = 0;
	int ReadCount_Conti=0;
	llong Ascii[256] = {0};

	map <int,int > BinDisTemp ;


	while ((l = kseq_read(seq)) >= 0) 
	{
		string RefSeq=(seq->seq.s);
		llong seq_length = (seq->seq.l);
		if  ( seq_length   <  InInt1 )
		{
			continue ;
		}

		stat_str_base(RefSeq , Ascii , seq_length);
		BaseCount+=seq_length ;
		ReadCount++;
		map <llong , int  > ::iterator it=ChrLen.find(seq_length);
		int AA=int(log10(seq_length) );
		BinDis[AA]++;
		BinDisTemp[AA]++;
		if (it==ChrLen.end())
		{
			ChrLen.insert(map <llong,int > ::value_type(seq_length,1));
		}
		else
		{
			(it->second)++;
		}

		//////////////     Conti  Stat //////////////
		vector<string> inf;
		split(RefSeq,inf,"Nn");
		int  a=inf.size();
		for(int ii=0 ; ii<a ; ii++)
		{

			long seq_length_Conti =inf[ii].length();
			BaseCount_Conti+=seq_length_Conti ;
			ReadCount_Conti++;
			map <llong , int  > ::iterator it=ChrLen_Conti.find(seq_length_Conti);
			int AA=int(log10(seq_length_Conti) );
			BinDis_Conti[AA]++;
			BinDisTemp[AA]++;
			if (it==ChrLen_Conti.end())
			{
				ChrLen_Conti.insert(map <llong,int > ::value_type(seq_length_Conti,1));
			}
			else
			{
				(it->second)++;
			}
		}
	}
	kseq_destroy(seq);
	gzclose(fp);




	(*OUT)<<"\n##########\tsummarize Scaffold report\t##########"<<endl;

	llong  BaseSum_A=Ascii['A']+Ascii['a'] ;
	llong  BaseSum_T=Ascii['T']+Ascii['t'] ;
	llong  BaseSum_C=Ascii['C']+Ascii['c'] ;
	llong  BaseSum_G=Ascii['G']+Ascii['g'] ;
	llong  BaseSum_N=Ascii['N']+Ascii['n'] ;

	map <llong , int  > ::iterator it=ChrLen.end();
	it--;
	llong maxReadLeng=it->first;
	llong MeanLenth=int(BaseCount/ReadCount);
	(*OUT)<<"#ScaNum: "<<ReadCount<<"\tBaseNum: "<<BaseCount<<"\n#MaxScaLen: "<<maxReadLeng<<"\tMeanScaLen: "<<MeanLenth<<endl;

	double pc=BaseSum_C*100.0/BaseCount ;
	double pg=BaseSum_G*100.0/BaseCount ;

	double pN=BaseSum_N*100.0/BaseCount ;
	double pA=BaseSum_A*100.0/BaseCount ;
	double pT=BaseSum_T*100.0/BaseCount ;

	double pGC=(BaseSum_C+BaseSum_G)*100.0/(BaseSum_C+BaseSum_G+BaseSum_A+BaseSum_T);

	(*OUT)<<"\n#GC%: "<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<pGC<<"%\tAT%: "<<100-pGC<<"%"<<endl;
	(*OUT)<<"#A BaseNum: "<<BaseSum_A<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<pA<<"%"<<endl;
	(*OUT)<<"#C BaseNum: "<<BaseSum_C<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<pc<<"%"<<endl;
	(*OUT)<<"#T BaseNum: "<<BaseSum_T<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<pT<<"%"<<endl;
	(*OUT)<<"#G BaseNum: "<<BaseSum_G<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<pg<<"%"<<endl;
	(*OUT)<<"#N BaseNum: "<<BaseSum_N<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<pN<<"%"<<endl;


	///////////////////////////////
	(*OUT)<<"\n##########\tsummarize ConTig report\t##########"<<endl;

	it=ChrLen_Conti.end();
	it--;
	llong maxReadLeng_Conti=it->first;
	llong MeanLenth_Conti=int(BaseCount_Conti/ReadCount_Conti);
	(*OUT)<<"#ContigNum: "<<ReadCount_Conti<<"\tBaseNum: "<<BaseCount_Conti<<"\n#MaxContigLen: "<<maxReadLeng_Conti<<"\tMeanContigLen: "<<MeanLenth_Conti<<endl;

	pc=BaseSum_C*100.0/BaseCount_Conti ;
	pg=BaseSum_G*100.0/BaseCount_Conti ;
	pA=BaseSum_A*100.0/BaseCount_Conti ;
	pT=BaseSum_T*100.0/BaseCount_Conti ;


	(*OUT)<<"\n#GC%: "<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<pc+pg<<"%\tAT%: "<<pA+pT<<"%"<<endl;
	(*OUT)<<"#A BaseNum: "<<BaseSum_A<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<pA<<"%"<<endl;
	(*OUT)<<"#C BaseNum: "<<BaseSum_C<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<pc<<"%"<<endl;
	(*OUT)<<"#T BaseNum: "<<BaseSum_T<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<pT<<"%"<<endl;
	(*OUT)<<"#G BaseNum: "<<BaseSum_G<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<pg<<"%"<<endl;





	map <int , int  > ::iterator BinIt=BinDisTemp.begin();
	llong	ScaCount_Temp=0;
	llong	ConTigCount_Temp=0;

	(*OUT)<<"\n##########\tsummarize Length distribution\t##########"<<endl;
	(*OUT)<<"\n#LenDis\tScafflod_Number\t>End:ScaPer%\tConTig_Number\t>End:ConTigPer%"<<endl;
	for ( ; BinIt!=BinDisTemp.end(); BinIt++)
	{
		int Start=int(pow(10.0,BinIt->first));
		int End=int(pow(10.0,(BinIt->first)+1));

		map <int , int  > ::iterator ItScaffold=BinDis.find(BinIt->first);
		double P_scaf=0.0;
		int  Number_scaf=0;

		if (ItScaffold!=BinDis.end())
		{
			Number_scaf=ItScaffold->second ;
		}

		ScaCount_Temp+=Number_scaf;
		P_scaf=(Number_scaf)*100.0/ReadCount ;
		double PP_scaf=100-(ScaCount_Temp)*100.0/ReadCount;


		map <int , int  > ::iterator ItConTig=BinDis_Conti.find(BinIt->first);
		double P_ConTig=0.0;
		int Number_ConTig=0;

		if (ItConTig!=BinDis_Conti.end())
		{
			Number_ConTig=ItConTig->second ;
		}

		ConTigCount_Temp+=Number_ConTig;
		P_ConTig=(Number_ConTig)*100.0/ReadCount_Conti ;
		double PP_ConTig=100-(ConTigCount_Temp)*100.0/ReadCount_Conti;


		(*OUT)<<Start<<"--"<<End<<" \t"<<Number_scaf<<"("<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<P_ConTig<<"%)\t>"<<End<<": "<<PP_scaf<<"%\t"<<Number_ConTig<<"("<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<P_ConTig<<"%)\t>"<<End<<": "<<PP_ConTig<<"%"<<endl;
	}



	(*OUT)<<"\n##########\tsummarize N50 distribution\t##########"<<endl;
	(*OUT)<<"\n#NX\tScaffold_Length(bp)\tScaffold_Number(%)\tConTig_Length(bp)\tConTig_Number(%)"<<endl;
	it=ChrLen.end();
	int NX=10;
	llong BaseCount_Temp=0;
	ScaCount_Temp=0 ;
	map <int,llong > ScaAAA;
	map <int,llong > ScaBBB;

	while(1)
	{
		it--;
		ScaCount_Temp+=it->second;
		BaseCount_Temp+=(it->first)*(it->second);
		if (BaseCount_Temp>BaseCount*NX/100)
		{
			ScaAAA[NX]=it->first;
			ScaBBB[NX]=ScaCount_Temp;
			NX=(int(BaseCount_Temp*10.0/BaseCount)+1)*10;
		}
		if (it==ChrLen.begin())
		{
			break;
		}
	}


	it=ChrLen_Conti.end();
	NX=10;
	BaseCount_Temp=0;
	ScaCount_Temp=0 ;
	int ScaAAA_TmpC=-1;
	int ScaBBB_TmpC=-1;
	while(1)
	{
		it--;
		ScaCount_Temp+=it->second;
		BaseCount_Temp+=(it->first)*(it->second);
		if (BaseCount_Temp>BaseCount_Conti*NX/100)
		{
			double P_Contig=ScaCount_Temp*100.0/ReadCount_Conti;
			if (ScaAAA.find(NX)!=ScaAAA.end())
			{
				ScaBBB_TmpC=ScaBBB[NX];
				ScaAAA_TmpC=ScaAAA[NX];

				double P_Sca=ScaBBB_TmpC*100.0/ReadCount ;
				(*OUT)<<"N"<<NX<<"\t"<<ScaAAA_TmpC<<"\t"<<ScaBBB_TmpC<<"("<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<P_Sca<<"%)\t"<<it->first<<"\t"<<ScaCount_Temp<<"("<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<P_Contig<<"%)"<<endl;
			}
			else
			{
				if  (ScaAAA_TmpC==-1)
				{
					(*OUT)<<"N"<<NX<<"\t-\t-(NA%)\t"<<it->first<<"\t"<<ScaCount_Temp<<"("<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<P_Contig<<"%)"<<endl;
				}
				else
				{
					double P_Sca=ScaBBB_TmpC*100.0/ReadCount ;
					(*OUT)<<"N"<<NX<<"\t"<<ScaAAA_TmpC<<"\t"<<ScaBBB_TmpC<<"("<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<P_Sca<<"%)\t"<<it->first<<"\t"<<ScaCount_Temp<<"("<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<P_Contig<<"%)"<<endl;
				}
			}
			NX=(int(BaseCount_Temp*10.0/BaseCount_Conti)+1)*10;
		}

		if (it==ChrLen_Conti.begin())
		{
			break;
		}

	}


}

#endif // RefN50Stat_H_

///////// swimming in the sky and flying in the sea ////////////
