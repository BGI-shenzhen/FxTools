#ifndef FQ_Stat_H_
#define FQ_Stat_H_


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
#include "../ALL/kseq.h"
#include <gzstream.h>
#include "FQ_Filter.h" 
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/bind.hpp>
#include <mutex>

using namespace std;
typedef long long llong ;

boost::mutex  mutexOUT ;
///////////////////


int  print_Ausage_A30()
{
	cout <<""
		"\n"
		"\tUsage: stat -i <1.fq> <2.fq>  -o <info.out>\n"
		"\n"
		"\t\t-i    <str>   [Repeat] Input Fq File for stat\n"
		"\t\t-l    <str>   Input Fq File List for stat\n"
		"\n"
		"\t\t-o    <str>   Out file of stat,otherwise[STDOUT]\n"
		"\t\t-p    <int>   Number of Thread [1]\n" 
		"\n"
		"\t\t-h            show this help\n" 
		"\n";
	return 1;
}


int parse_Acmd_A30(int argc, char **argv,In3str1v *  para_A30 )
{
	if (argc <=2 ) {print_Ausage_A30();return 0 ;}
	int file_count=0;
	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0 ;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InFq" ||  flag  == "i")
		{
			if(i + 1 == argc) { LogLackArg(flag); return 0 ;}
			i++;
			string A=argv[i];
			(para_A30->List).push_back(A);
			file_count++;
			bool RunT=true;
			while(RunT)
			{
				if ( (i+1) < argc  && (argv[i+1][0]!='-'))
				{
					i++;
					A=argv[i];
					(para_A30->List).push_back(A);
					file_count++;
				}
				else
				{
					RunT=false;
				}
			}
		}
		else if (flag  ==  "InFqList" || flag  == "l")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0 ;}
			i++;
			string A=argv[i];
			file_count+=(ReadList (A ,(para_A30->List)));
		}
		else if (flag  ==  "CPU" || flag  == "p")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0 ;}
			i++;
			para_A30->InInt=atoi(argv[i]);
		}
		else if (flag  ==  "OutStat" || flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0 ;}
			i++;
			para_A30->InStr2=argv[i];
		}
		else if (flag  == "help" || flag  == "h")
		{
			print_Ausage_A30();return 0 ;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0 ;
		}
	}
	if  ( file_count<1 )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0 ;
	}
	return 1 ;
}


void StatFq ( string &file , std::ostream * OUT ,  In3str1v  * para_A30 )
{
	string ext =file.substr(file.rfind('/') ==string::npos ? file.length() : file.rfind('/') + 1);
	if (ext.empty())
	{
		ext=file ;
	}
	igzstream INFQ (file.c_str(),ifstream::in);
	int UC[256]={4};
	UC['A']=0; UC['C']=1;UC['G']=3;UC['T']=2 ; UC['N']=4;
	UC['a']=0; UC['c']=1;UC['g']=3;UC['t']=2 ; UC['n']=4;

	int  minBaseQ=64 ;
	if (!(para_A30->TF))
	{
		minBaseQ=33;
	}
	if(!INFQ.good())
	{
		cerr << "open InputFQ  error: "<<file<<endl;
		return ;
	}

	int maxReadLeng=0 ;
	int BB=0;
	while( (!INFQ.eof()) && (BB< 16888 ))
	{
		string  line ;
		getline(INFQ,line);
		getline(INFQ,line);
		int tmp=line.length();
		if (tmp> maxReadLeng)
		{
			maxReadLeng=tmp;
		}
		getline(INFQ,line);
		getline(INFQ,line);
		BB++;
	}
	INFQ.close();
	/*
	   if ( maxReadLeng==90 || maxReadLeng==100 ||  maxReadLeng==150 || maxReadLeng==50 || maxReadLeng==1000 || maxReadLeng==500 || maxReadLeng==200 ||  maxReadLeng==75 ||  maxReadLeng==45 ||  maxReadLeng==44 ||   maxReadLeng==35  ||   maxReadLeng==125  )
	   {

	   }
	   else 
	   {
	   cerr<<"Warn:: In the head 16888 Reads, the maxReadLength is\t"<<maxReadLeng<<endl;
	   cerr<<"   Which is not general Readlength:35/45/50/75/90/100/125/150"<<endl;
	   }
	   */
	igzstream INFQ2 (file.c_str(),ifstream::in);
	if(!INFQ2.good())
	{
		cerr << "open InputFQ  error: "<<file<<endl;
		return ;
	}

	DA *BaseD =new DA[maxReadLeng];
	llong ReadCount=0 , BaseCount=0 ;
	llong QCount[10]={0};
	llong RCount[10]={0};
	llong BBaseCount[128]={0};

	while( !INFQ2.eof() )
	{
		string  ID , seq ,tmp ,Qseq ;
		getline(INFQ2,ID);
		if (ID.length()<=0)  { continue ; }
		getline(INFQ2,seq);  getline(INFQ2,tmp); getline(INFQ2,Qseq);
		int ReadLeng=seq.length();
		ReadCount++;  BaseCount+=ReadLeng ;
		int SumQ=0 ;
		if (ReadLeng==maxReadLeng)
		{
			for (int i=0 ; i<ReadLeng ; i++ )
			{
				BBaseCount[seq[i]]++;
				SumQ+=Qseq[i];
				int qbase=(Qseq[i]-minBaseQ)/10 ;
				QCount[qbase]++;
				int ey=UC[(seq[i])];
				(BaseD[i]).ADD(ey);
			}
		}
		else
		{
			for (int i=0 ; i<ReadLeng ; i++ )
			{
				BBaseCount[seq[i]]++;
				SumQ+=Qseq[i];
				int qbase=(Qseq[i]-minBaseQ)/10 ;
				QCount[qbase]++;
			}
		}
		int  Rq=((SumQ/ReadLeng)-minBaseQ)/10 ;
		RCount[Rq]++;
	}

	INFQ2.close();

	boost::mutex::scoped_lock lock(mutexOUT);

	(*OUT)<<"##"<<ext<<"##"<<endl ;
	(*OUT)<<"#ReadNum: "<<ReadCount<<"\tBaseNum: "<<BaseCount<<"\tReadLeng: "<<maxReadLeng<<endl;

	llong BaseSum_A=(BBaseCount['A']+BBaseCount['a']);
	llong BaseSum_T=(BBaseCount['T']+BBaseCount['t']);
	llong BaseSum_C=(BBaseCount['C']+BBaseCount['c']);
	llong BaseSum_G=(BBaseCount['G']+BBaseCount['g']);
	llong BaseSum_N=(BBaseCount['N']+BBaseCount['n']);

	double pc=BaseSum_C*100.0/BaseCount ;
	double pg=BaseSum_G*100.0/BaseCount ;

	double pN=BaseSum_N*100.0/BaseCount ;
	double pA=BaseSum_A*100.0/BaseCount ;
	double pT=BaseSum_T*100.0/BaseCount ;

	(*OUT)<<"#GC%: "<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<(pc+pg)<<"%\tAT%: "<<(pA+pT)<<"%"<<endl;

	(*OUT)<<"#A BaseNum: "<<BaseSum_A<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<pA<<"%"<<endl;
	(*OUT)<<"#C BaseNum: "<<BaseSum_C<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<pc<<"%"<<endl;
	(*OUT)<<"#T BaseNum: "<<BaseSum_T<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<pT<<"%"<<endl;
	(*OUT)<<"#G BaseNum: "<<BaseSum_G<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<pg<<"%"<<endl;
	(*OUT)<<"#N BaseNum: "<<BaseSum_N<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<pN<<"%"<<endl;

	(*OUT)<<endl;

	double sump=0;
	for (int iii=0 ;iii<10 ; iii++ )
	{
		if  ( QCount[iii]==0 )
		{
			continue ;
		}
		int Start=iii*10;
		int End=Start+10 ;
		double p=QCount[iii]*100.0/BaseCount ;
		sump+=p;
		(*OUT)<<"#BaseQ:"<<Start<<"--"<<End<<" : "<<QCount[iii]<<"("<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<p<<"%)"<<"\t>Q"<<End<<": "<<(100-sump)<<"%"<<endl;
	}
	(*OUT)<<endl;

	sump=0;
	for (int iii=0 ;iii<10 ; iii++ )
	{
		if  ( RCount[iii]==0 )
		{
			continue ;
		}
		int Start=iii*10;
		int End=Start+10 ;
		double p=(RCount[iii])*100.0/ReadCount ;
		sump+=p;
		double Big_p=(100-sump) ;
		if (Big_p<0) {Big_p=0 ;}
		(*OUT)<<"#ReadQ:"<<Start<<"--"<<End<<" : "<<RCount[iii]<<"("<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<p<<"%)"<<"\t>Q"<<End<<": "<<Big_p<<"%"<<endl;
	}

	(*OUT)<<"Posi\tA\tC\tT\tG\tN"<<endl;
	for (int ii=0 ; ii<maxReadLeng; ii++)
	{
		(*OUT)<<ii+1<<"\t"<<(BaseD[ii]).get(0)<<"\t"<<(BaseD[ii]).get(1)<<"\t"<<(BaseD[ii]).get(2)<<"\t"<<(BaseD[ii]).get(3)<<"\t"<<(BaseD[ii]).get(4)<<"\t"<<endl;
	}

	delete [] BaseD ;
	return ;
}


int FQ_Stat_main( int argc,char *argv[] )
//int main( int argc,char *argv[] )
{
	In3str1v *  para_A30 =new In3str1v ;
	para_A30->InInt=1 ;
	if( parse_Acmd_A30( argc, argv,para_A30) ==0 )
	{
		delete  para_A30 ;
		return 1 ;
	}
	int tempshiftQ=GetShiftQ( ((para_A30->List)[0]) );

	if (tempshiftQ !=64)
	{
		para_A30->TF=false ;
	}


	std::ostream * Info = &cout;
	std::ofstream fout;
	if (!(para_A30->InStr2).empty())
	{
		fout.open((para_A30->InStr2).c_str());
		Info = &fout;

	}



	int List_Acount=(para_A30->List).size();
	if  ((para_A30->InInt==1) || (List_Acount==1 ))
	{    
		for (int i=0; i<List_Acount ; i++)
		{
			string FQStat_ANow=(para_A30->List)[i];
			StatFq (FQStat_ANow , Info ,  para_A30 );
		}
	}
	else
	{
		int i=0;
		while(true)
		{
			boost::thread_group TRead ;

			for ( int j=0 ;j<(para_A30->InInt) ;j++ )
			{
				if (i<List_Acount)
				{
					string FQStat_ANow=(para_A30->List)[i];
					i++;
					TRead.create_thread(bind(StatFq,FQStat_ANow,boost::ref(Info),boost::ref(para_A30)));
				}
				else
				{
					break ;
				}
			}
			TRead.join_all();
			if (i>=List_Acount)
			{
				break ;
			}
		}
	}
	delete para_A30 ;

	//(*Info).close();
	return 0 ;
}

///////// swimming in the sky and flying in the sea ////////////
#endif 
