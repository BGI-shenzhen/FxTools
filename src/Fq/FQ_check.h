#ifndef FQ_CheckAdpter_H_
#define FQ_CheckAdpter_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <vector>
#include <iomanip>
#include <math.h>
#include <zlib.h>
#include <gzstream.h>
#include "FQ_Filter.h"

#define MAX_L 512
#define MAX_Q 128
using namespace std;

//  #define _A 0 #define _C 1    #define _G 2   #define _T 3   #define _N 4

int  print_usage_FqCheck()
{
	cout <<""
		"\n"
		"Usage:fqcheck  -i A1.fq  A2.fq -o 1.fqcheck 2.fqcheck  [options]\n"
		"\n"
		"\t\t-i     <str>   File name of InFq1/InFq2 Input\n"
		"\t\t-o     <str>   Prefix of OUT File1/File2 name\n"
		"\n"
		"\t\t-a     <str>   Input adapters1/adapters2 fa file\n"
		"\n"
		"\t\t-h             show this help [v2.05]\n"
		"\n";
	return 1;
}

int parse_cmd_FqCheck(int argc, char **argv , Para_A24 * para ,   Para_A24 * P2In  )
{
	if (argc <=4  ) {print_usage_FqCheck();return 0;}

	int err_flag = 0;
	for(int i = 1; i < argc || err_flag; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");
		if (flag  == "InFq1" || flag  == "r" ||  flag  == "i" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->InFq1=argv[i];
			if ( (i+1) < argc  && (argv[i+1][0]!='-'))
			{
				i++;
				P2In->InFq1=argv[i];
			}
		}
		else if (flag  ==  "Adapter1"   || flag  == "a")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->InFq2=argv[i];
			if ( (i+1) < argc  && (argv[i+1][0]!='-'))
			{
				i++;
				P2In->InFq2=argv[i];
			}
		}
		else if (flag  ==  "OutStat1" || flag  == "c"   || flag  == "o")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->OutFq1=argv[i];
			if  ( (i+1) < argc && (argv[i+1][0]!='-'))
			{
				i++;
				P2In->OutFq1=argv[i];
			}
		}
		else if (flag  == "InFq2" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->InFq1=argv[i];
		}
		else if (flag  ==  "Adapter2" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->InFq2=argv[i];
		}
		else if (flag  ==  "OutStat2" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->OutFq1=argv[i];
		}
		else if (flag  == "help"   ||  flag  == "h")
		{
			print_usage_FqCheck();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if ( (para->OutFq1).empty() || (para->InFq1).empty() )
	{
		cerr<< "-InFq1 -OutStat1 lack argument for the must"<<endl;
		return 0;
	}

	if  ((P2In->OutFq1).empty() &&  (P2In->InFq1).empty())
	{
		return 1 ;
	}
	else if ( (P2In->OutFq1).empty() || (P2In->InFq1).empty() )
	{
		cerr<< "-InFq2 -OutStat2 must set togther"<<endl;
		return 0;
	}
	else
	{
		return  2 ;
	}
}



int RunSingleFq( string InFq, string OutStat,  Para_A24 * para  )
{

	static	long int base[5];//记录每种碱基的总数
	static	long int qual[MAX_Q];//每个质量值的总频率
	static	long int cycle_qual[MAX_L][MAX_Q];//每个cycle的质量分布
	static	long int cycle_base[MAX_L][5];//每个cycle的碱基分布
	static	double   q_err[MAX_Q];//记录质量为q的碱基的错误率
	static	double   err = 0.;
	static	long int q20 = 0;
	static	long int q30 = 0;
	int read_num = 0 ;

	int UC[256]={4};
	UC['A']=0; UC['C']=1;UC['G']=2;UC['T']=3 ; UC['N']=4;
	UC['a']=0; UC['c']=1;UC['g']=2;UC['t']=3 ; UC['n']=4;

	read_num = 0 ;
	q20 = 0;
	q30 = 0;
	err = 0.0;
	base[0]=0 ; base[1]=0 ; base[2]=0 ; base[3]=0 ; base[4]=0 ;

	for (int i=0 ; i<MAX_Q ; i++)
	{
		qual[i]=0;
		q_err[i]=0;
		for (int j=0 ; j<MAX_L ; j++)
		{
			cycle_qual[j][i]=0;
		}
	}

	for (int j=0 ; j<MAX_L ; j++)
	{
		cycle_base[j][0]=0;
		cycle_base[j][1]=0;
		cycle_base[j][2]=0;
		cycle_base[j][3]=0;
		cycle_base[j][4]=0;
	}

	//计算所有可能的质量值的err
	for(int q=0;q<MAX_Q;q++) {q_err[q]=1.0/(pow((double)10.0,(double)(q*0.1)));}
	{
		q_err[2]=0.121821900;
		q_err[3]=0.501187233;
		q_err[4]=0.398107170;
		q_err[5]=0.316227766;
		q_err[6]=0.251188643;
		q_err[7]=0.199526232;
		q_err[8]=0.158489319;
		q_err[9]=0.128647300;
		q_err[10]=0.10545720;
		q_err[11]=0.09508237;
		q_err[12]=0.07943282;
		q_err[13]=0.07539894;
		q_err[14]=0.03981071;
		q_err[15]=0.03162277;
		q_err[16]=0.02511886;
		q_err[17]=0.01995262;
		q_err[18]=0.01629471;
		q_err[19]=0.01584893;
		q_err[20]=0.01470216;
		q_err[21]=0.00910940;
		q_err[22]=0.00791782;
		q_err[23]=0.00601007;
		q_err[24]=0.00590479;
		q_err[25]=0.00559080;
		q_err[26]=0.00518726;
		q_err[27]=0.00437646;
		q_err[28]=0.00405721;
		q_err[29]=0.00393885;
		q_err[30]=0.00387256;
		q_err[31]=0.00341158;
		q_err[32]=0.00337891;
		q_err[33]=0.00296104;
		q_err[34]=0.00282249;
		q_err[35]=0.00261006;
		q_err[36]=0.00247926;
		q_err[37]=0.00227605;
		q_err[38]=0.00217025;
		q_err[39]=0.00215136;
		q_err[40]=0.00188753;
	}

	//input file
	char read_id[MAX_L], read_se[MAX_L], qual_id[MAX_L], qual_se[MAX_L];
	int read_len;

	gzFile read_fs;
	read_fs = gzopen (InFq.c_str(), "rb");
	if (read_fs==NULL)
	{
		cerr<<"Cann't open files\t"<<InFq<<endl;
		return 1;
	}

	gzgets(read_fs, read_id, MAX_L);
	gzgets(read_fs, read_se, MAX_L);
	string read = read_se;
	read_len = read.length()-1;
	gzrewind(read_fs);

	//output file
	ofstream  OUTFile;
	OUTFile.open (OutStat.c_str());
	if (!OUTFile.good())
	{
		cout<<"Cann't open files\t"<<OutStat<<endl;
		return 1 ;
	}

	while (gzgets (read_fs, read_id, MAX_L)) 
	{
		gzgets (read_fs, read_se, MAX_L);
		gzgets (read_fs, qual_id, MAX_L);
		gzgets (read_fs, qual_se, MAX_L);
		++read_num;
		for (int c=0;c<read_len;c++) 
		{
			++cycle_base[c][UC[read_se[c]]];

			int q = qual_se[c]-(para->LowQint);
			if(q<0) {q=0;}
			++cycle_qual[c][q];
		}
	}
	gzclose(read_fs);


	long base_total = (long) read_len * read_num;
	for(int c=0;c<read_len;c++)
	{
		base[0]+=cycle_base[c][0];
		base[1]+=cycle_base[c][1];
		base[2]+=cycle_base[c][2];
		base[3]+=cycle_base[c][3];
		base[4]+=cycle_base[c][4];

	}
	for(int c=0;c<read_len;c++)
	{
		for(int q=0;q<MAX_Q;q++)
		{
			qual[q]+=cycle_qual[c][q];
		}
	}
	int qMAX=0;
	for(int q=0;q<MAX_Q;q++)
	{
		if (qual[q])
		{
			qMAX=q;
			if (q>=30) q30+=qual[q];
			if (q>=20) q20+=qual[q];
			err+=q_err[q]*qual[q];
		}
	}

	OUTFile<<" the default quality shift value is: -"<<(para->LowQint)<<", "
		<<read_num<<" sequences, "
		<<base_total<<" total length, Max length:"
		<<read_len<<", average length:"
		<<setiosflags(ios::fixed)<<setprecision(2)<<(double)read_len<<endl;
	OUTFile<<"Standard deviations at 0.25:  total "<<100*(sqrt(0.25*(double)base_total)/base_total)
		<<"%, per base "<<100*(sqrt(0.25*(double)read_num)/read_num)
		<<"%"<<endl;
	OUTFile<<"             A     C     G     T     N ";

	for(int q=0;q<=qMAX;q++)
	{
		OUTFile<<setw(4)<<q<<' ';
	}

	OUTFile<<endl;
	OUTFile<<"Total    ";

	for(int b=0;b<=4;b++)  
	{
		OUTFile<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<100*(double)base[b]/base_total<<' ';
	}

	for(int q=0;q<=qMAX;q++) 
	{
		OUTFile<<setw(4)<<(int)(1000*((double)qual[q]/base_total))<<" ";
	}
	OUTFile<<endl;

	for(int c=0;c<read_len;c++)	
	{
		OUTFile<<"base "<<setw(3)<<c+1<<' ';
		for(int b=0;b<=4;b++)  {OUTFile<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<100*(double)cycle_base[c][b]/read_num<<' ';}
		for(int q=0;q<=qMAX;q++)
		{
			OUTFile<<setw(4)<<(int)(1000*((double)cycle_qual[c][q]/read_num))<<' ';
		}
		OUTFile<<endl;
	}

	OUTFile<<endl;
	OUTFile<<"Error Rate\t%GC\tQ20\tQ30"<<endl;
	OUTFile.precision(2);
	OUTFile<<(100.0*err)/(base_total)<<"\t"<<(100.0*(base[1]+base[2]))/(base_total-base[4])<<"\t"<<(100.0*q20)/base_total<<"\t"<<(100.0*q30)/base_total<<endl;
	OUTFile.close();


	return 0;
}



int RunSingleFqAdpter( string InFq, string OutStat,  Para_A24 * para  )
{

	static	long int base[5];//记录每种碱基的总数
	static	long int qual[MAX_Q];//每个质量值的总频率
	static	long int cycle_qual[MAX_L][MAX_Q];//每个cycle的质量分布
	static	long int cycle_base[MAX_L][5];//每个cycle的碱基分布
	static	double   q_err[MAX_Q];//记录质量为q的碱基的错误率
	static	double   err = 0.;
	static	long int q20 = 0;
	static	long int q30 = 0;
	int read_num = 0;
	int empt_num= 0;
	int adpt_num =0;


	int UC[256]={4};
	UC['A']=0; UC['C']=1;UC['G']=2;UC['T']=3 ; UC['N']=4;
	UC['a']=0; UC['c']=1;UC['g']=2;UC['t']=3 ; UC['n']=4;

	read_num = 0;
	q20 = 0;
	q30 = 0;
	err = 0.0;
	base[0]=0 ; base[1]=0 ; base[2]=0 ; base[3]=0 ; base[4]=0 ;

	for (int i=0 ; i<MAX_Q ; i++)
	{
		qual[i]=0;
		q_err[i]=0;
		for (int j=0 ; j<MAX_L ; j++)
		{
			cycle_qual[j][i]=0;
		}
	}

	for (int j=0 ; j<MAX_L ; j++)
	{
		cycle_base[j][0]=0;
		cycle_base[j][1]=0;
		cycle_base[j][2]=0;
		cycle_base[j][3]=0;
		cycle_base[j][4]=0;
	}

	//计算所有可能的质量值的err
	for(int q=0;q<MAX_Q;q++) {q_err[q]=1.0/(pow((double)10.0,(double)(q*0.1)));}
	{
		q_err[2]=0.121821900;
		q_err[3]=0.501187233;
		q_err[4]=0.398107170;
		q_err[5]=0.316227766;
		q_err[6]=0.251188643;
		q_err[7]=0.199526232;
		q_err[8]=0.158489319;
		q_err[9]=0.128647300;
		q_err[10]=0.10545720;
		q_err[11]=0.09508237;
		q_err[12]=0.07943282;
		q_err[13]=0.07539894;
		q_err[14]=0.03981071;
		q_err[15]=0.03162277;
		q_err[16]=0.02511886;
		q_err[17]=0.01995262;
		q_err[18]=0.01629471;
		q_err[19]=0.01584893;
		q_err[20]=0.01470216;
		q_err[21]=0.00910940;
		q_err[22]=0.00791782;
		q_err[23]=0.00601007;
		q_err[24]=0.00590479;
		q_err[25]=0.00559080;
		q_err[26]=0.00518726;
		q_err[27]=0.00437646;
		q_err[28]=0.00405721;
		q_err[29]=0.00393885;
		q_err[30]=0.00387256;
		q_err[31]=0.00341158;
		q_err[32]=0.00337891;
		q_err[33]=0.00296104;
		q_err[34]=0.00282249;
		q_err[35]=0.00261006;
		q_err[36]=0.00247926;
		q_err[37]=0.00227605;
		q_err[38]=0.00217025;
		q_err[39]=0.00215136;
		q_err[40]=0.00188753;
	}


	igzstream ADF ((para->InFq2).c_str(),ifstream::in);
	if(!ADF.good())
	{
		cerr << "open Adpter File error: "<<para->InFq2<<endl;
		return 1;
	}
	string adpt_id, adpt_se;
	getline(ADF,adpt_id);	adpt_id.erase(0,1);
	getline(ADF,adpt_se);
	int adpt_len = adpt_se.length();
	ADF.close();

	string OUT_ListAdp=(para->OutFq1)+".adapter.list.gz";
	ogzstream OUT_List (OUT_ListAdp.c_str());
	string OUT_statAdp=(para->OutFq1)+".adapter.stat";
	ofstream  OUT_Stat (OUT_statAdp.c_str());


	string title = "#reads_id   reads_len   reads_start   reads_end   adapter_id   adapter_len   adapter_start   adapter_end   align_len   mismatch";
	OUT_List<<title<<endl;
	OUT_Stat<<"adapter_id\tpolluted_reads\tempty_reads\tadapter_sequence"<<endl;


	//input file
	char read_id[MAX_L], read_se[MAX_L], qual_id[MAX_L], qual_se[MAX_L];
	int read_len;

	gzFile read_fs;
	read_fs = gzopen (InFq.c_str(), "rb");
	if (read_fs==NULL)
	{
		cerr<<"Cann't open files\t"<<InFq<<endl;
		return 1;
	}

	gzgets(read_fs, read_id, MAX_L);
	gzgets(read_fs, read_se, MAX_L);
	string read = read_se;
	read_len = read.length()-1;
	gzrewind(read_fs);

	//output file
	ofstream  OUTFile;
	OUTFile.open (OutStat.c_str());
	if (!OUTFile.good())
	{
		cout<<"Cann't open files\t"<<OutStat<<endl;
		return 1 ;
	}

	int Find_Len=read_len-9;

	while (gzgets (read_fs, read_id, MAX_L)) 
	{
		gzgets (read_fs, read_se, MAX_L);
		gzgets (read_fs, qual_id, MAX_L);
		gzgets (read_fs, qual_se, MAX_L);
		read_num++;
		for (int c=0;c<read_len;c++)
		{
			++cycle_base[c][UC[read_se[c]]];
			int q = qual_se[c]-(para->LowQint);
			if(q<0) {q=0;}
			++cycle_qual[c][q];
		}


		bool find=true;
		int a1 = adpt_len-15;
		int r1 = 0;
		int len;
		int mis;

		for ( ;a1>0 ; a1--)
		{
			int len1 = adpt_len - a1;
			len=(len1<read_len)? len1:read_len;
			mis = 0;
			int map[MAX_L];
			map[0]=0;
			for (int c=0; c<len; c++)
			{
				if (adpt_se[a1+c]==read_se[c]) {map[mis]++;}
				else {mis++;map[mis]=0;}
			}			


			int max_map=map[0];
			for (int c=1; c<=mis; c++)
			{
				if (map[c]>max_map) {max_map=map[c];}
			}

			if ( (max_map>14)    || (len>(mis*5)) )
			{
				find = false ;
				break;
			}

		}

		if (find)
		{
			for (r1=0; r1<Find_Len; r1++)
			{
				int len2 = read_len - r1;
				len=(adpt_len<len2)? adpt_len:len2;
				mis = 0;
				int map[MAX_L];
				map[0]=0;
				for (int c=0; c<len; c++)
				{
					if (adpt_se[c]==read_se[r1+c]) {map[mis]++;}
					else {mis++;map[mis]=0;}
				}

				int max_map=map[0];
				for (int c=1; c<=mis; c++)
				{
					if (map[c]>max_map) {max_map=map[c];}
				}
				if ( (max_map>14)   ||  (len>(mis*5))  )
				{
					find = false;
					break;
				}
			}

		}

		if  (find)
		{
		}
		else
		{
			string RID=read_id;
			RID=RID.substr(1,RID.length()-2);
			OUT_List<<RID<<'\t'<<read_len<<'\t'<<r1<<'\t'<<r1+len-1<<'\t'<<adpt_id<<'\t'<<adpt_len<<'\t'<<a1<<'\t'<<a1+len-1<<'\t'<<len<<'\t'<<mis<<"\n";
			adpt_num++;
			if (r1<=3)
			{
				empt_num++;
			}
		}




	}
	gzclose(read_fs);
	OUT_List.close();


	float adpt_percent = ((float) adpt_num)/read_num*100;
	float empt_percent = ((float) empt_num)/read_num*100;
	OUT_Stat<<adpt_id<<'\t'<<adpt_num<<" ("<<setprecision(3)<<fixed<<adpt_percent<<"%)"<<'\t'<<empt_num<<" ("<<empt_percent<<"%)"<<'\t'<<adpt_se<<endl;
	OUT_Stat<<"\ntotal_reads: "<<read_num<<endl;
	OUT_Stat.close();



	long base_total = (long) read_len * read_num;
	for(int c=0;c<read_len;c++)
	{
		base[0]+=cycle_base[c][0];
		base[1]+=cycle_base[c][1];
		base[2]+=cycle_base[c][2];
		base[3]+=cycle_base[c][3];
		base[4]+=cycle_base[c][4];
	}

	for(int c=0;c<read_len;c++)
	{
		for(int q=0;q<MAX_Q;q++)
		{
			qual[q]+=(cycle_qual[c][q]);
		}
	}
	int qMAX=0;
	for(int q=0;q<MAX_Q;q++)
	{
		if (qual[q])
		{
			qMAX=q;
			if (q>=30) q30+=qual[q];
			if (q>=20) q20+=qual[q];
			err+=q_err[q]*qual[q];
		}
	}

	OUTFile<<" the default quality shift value is: -"<<(para->LowQint)<<", "
		<<read_num<<" sequences, "
		<<base_total<<" total length, Max length:"
		<<read_len<<", average length:"
		<<setiosflags(ios::fixed)<<setprecision(2)<<(double)read_len<<endl;
	OUTFile<<"Standard deviations at 0.25:  total "<<100*(sqrt(0.25*(double)base_total)/base_total)
		<<"%, per base "<<100*(sqrt(0.25*(double)read_num)/read_num)
		<<"%"<<endl;
	OUTFile<<"             A     C     G     T     N ";

	for(int q=0;q<=qMAX;q++)
	{
		OUTFile<<setw(4)<<q<<' ';
	}

	OUTFile<<endl;
	OUTFile<<"Total    ";

	for(int b=0;b<=4;b++)
	{
		OUTFile<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<100*(double)base[b]/base_total<<' ';
	}

	for(int q=0;q<=qMAX;q++)
	{
		OUTFile<<setw(4)<<(int)(1000*((double)qual[q]/base_total))<<" ";
	}
	OUTFile<<endl;

	for(int c=0;c<read_len;c++)	
	{
		OUTFile<<"base "<<setw(3)<<c+1<<' ';
		for(int b=0;b<=4;b++)  {OUTFile<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<100*(double)cycle_base[c][b]/read_num<<' ';}
		for(int q=0;q<=qMAX;q++)
		{
			OUTFile<<setw(4)<<(int)(1000*((double)cycle_qual[c][q]/read_num))<<' ';
		}
		OUTFile<<endl;
	}

	OUTFile<<endl;
	OUTFile<<"Error Rate\t%GC\tQ20\tQ30"<<endl;
	OUTFile.precision(2);
	OUTFile<<(100.0*err)/(base_total)<<"\t"<<(100.0*(base[1]+base[2]))/(base_total-base[4])<<"\t"<<(100.0*q20)/base_total<<"\t"<<(100.0*q30)/base_total<<endl;
	OUTFile.close();


	return 0;
}






int PlotfqCheck ( vector  <string> infile , string outfile)
{
	string CCC=outfile+"/base.tmp";
	string DDD=outfile+"/base.gnuplot";
	string PS=outfile+"/base.ps";
	string PNG=outfile+"/base.png";
	string PDF=outfile+"/base.pdf";
	ofstream OUTInfo (CCC.c_str());
	ofstream OUTPlot (DDD.c_str());


	string CCC_S=outfile+"/qual.tmp";
	string DDD_S=outfile+"/qual.gnuplot";
	string PS_S=outfile+"/qual.ps";
	string PNG_S=outfile+"/qual.png";
	string PDF_S=outfile+"/qual.pdf";
	ofstream OUTInfo_S (CCC_S.c_str());
	ofstream OUTPlot_S (DDD_S.c_str());


	string Flag;
	int coutFlag=0;
	double NumA,NumC,NumT,NumG,NumN,NumX;
	int FileNum=infile.size();
	string EndFlag="";

	int QQ ;
	for (int kk=0 ; kk< FileNum ; kk++)
	{
		igzstream IN (infile[kk].c_str(),ifstream::in);
		while (!IN.eof())
		{
			string line;
			getline(IN,line);
			if (line.length()<=0)  { continue ;}
			istringstream isone (line,istringstream::in);
			isone>>Flag;
			if  (Flag !="base")
			{
				continue  ;
			}
			coutFlag++;
			isone>>NumX>>NumA>>NumC>>NumG>>NumT>>NumN;
			OUTInfo<<coutFlag<<"\t"<<NumA<<"\t"<<NumC<<"\t"<<NumG<<"\t"<<NumT<<"\t"<<NumN<<endl;
			for (int jjjk=0 ; jjjk<42 ; jjjk++)
			{
				isone>>QQ;
				OUTInfo_S<<coutFlag<<"\t"<<jjjk<<"\t"<<QQ<<endl;
			}

		}
		IN.close();
		if (kk!=(FileNum-1))
		{
			EndFlag=EndFlag+Int2Str(coutFlag)+" 50\n";
		}
	}
	OUTInfo.close();
	OUTInfo_S.close();


	OUTPlot<<"reset;"<<endl;
	OUTPlot<<"set out \'"<<PS<<"\';"<<endl;
	OUTPlot<<"set xlabel 'Position along reads';\nset grid front lc rgb 'gray';\nset terminal postscript portrait color size 8, 5;\nset ylabel 'Percent';\nset title 'Base percentage composition along reads';\nset yrange [-2:52];\nset xrange [0:"<<coutFlag<<"+1];\n"<<"plot '"<<CCC<<"\' u 1:2 w l lw 3 t 'A', '' u 1:3 w l lw 3 t 'C', '' u 1:4 w l lw 3 t 'G', '' u 1:5 w l lw 3 t 'T', '' u 1:6 w l lw 3 t 'N', '-' w i lc rgb 'blue' lw 3 notitle;\n"<<EndFlag<<"e\nreset;"<<endl;
	OUTPlot.close();
	string path="gnuplot  "+DDD+" ;   convert " +PS+"  "+PNG+"; convert "  +PS+"  "+PDF  +" ; rm -rf "+CCC+" "+DDD+" "+PS+"; ";
	//string path="gnuplot  "+DDD+" ;   convert " +PS+"  "+PNG+" ; rm -rf "+CCC+" "+DDD+" "+PS+"; ";
	//	system(path.c_str()) ;

	popen(path.c_str(),"r") ;


	OUTPlot_S<<"reset;"<<endl;
	OUTPlot_S<<"set out \'"<<PS_S<<"\';"<<endl;
	OUTPlot_S<<"set xlabel 'Position along reads';\nset cbrange [0:100];\nunset colorbox;\nset palette defined (0 '#ffffff', 10 '#00ff00', 30 '#ffff00', 50 '#ff0000', 70 '#800000', 100 '#000000');\nset grid front lc rgb 'gray';\nset terminal postscript portrait color size 8, 5;\nset ylabel 'Quality';\nset title 'Base percentage composition along reads';\nset yrange [-2:42];\nset key off;\nset xrange [0:"<<coutFlag<<"+1];\n"<<"plot '"<<CCC_S<<"\' u 1:2:($3/10) w image;\nreset;"<<endl;
	OUTPlot_S.close();
	//path="gnuplot  "+DDD_S+" ;   convert " +PS_S+"  "+PNG_S+ " ; rm -rf "+CCC_S+" "+DDD_S+" "+PS_S+"; ";
	path="gnuplot  "+DDD_S+" ;   convert " +PS_S+"  "+PNG_S+"; convert "  +PS_S+"  "+PDF_S  +" ; rm -rf "+CCC_S+" "+DDD_S+" "+PS_S+"; ";
	//std::system(path.c_str()) ;
	popen(path.c_str(),"r") ;
}

int FQ_Check_main(int argc, char **argv)
//int main (int argc, char *argv[ ])
{

	Para_A24 * para_A24 = new Para_A24;
	Para_A24 * P2In = new Para_A24;
	para_A24->LowQint=64;	P2In->LowQint=64;

	if( parse_cmd_FqCheck(argc, argv, para_A24, P2In )==0)
	{
		delete  para_A24 ;
		delete P2In ;
		return 0 ;
	}
	string path=(para_A24->OutFq1);
	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);
	if (ext != "fqcheck")
	{
		path=path+".fqcheck" ;
	}
	else
	{
		(para_A24->OutFq1)=path.substr(0,path.rfind('.') ==string::npos ? path.length() : path.rfind('.'));
	}
	para_A24->LowQint=GetShiftQ((para_A24->InFq1));

	vector<string> FqCheck;
	FqCheck.push_back(path);
	if ( (para_A24->InFq2).empty() )
	{
		RunSingleFq(  (para_A24->InFq1) , path , para_A24 );
	}
	else
	{
		RunSingleFqAdpter(  (para_A24->InFq1) , path  , para_A24   );
	}

	if  ((P2In->InFq1).empty())
	{

	}
	else
	{
		path=(P2In->OutFq1);
		string ext2 =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);
		if (ext2 != "fqcheck")
		{
			path=path+".fqcheck" ;
		}
		else
		{
			(P2In->OutFq1)=path.substr(0,path.rfind('.') ==string::npos ? path.length() : path.rfind('.'));
		}

		P2In->LowQint=GetShiftQ((P2In->InFq1));

		FqCheck.push_back(path);
		if ( (P2In->InFq2).empty() )
		{
			RunSingleFq(  (P2In->InFq1) , path , P2In );
		}
		else
		{
			RunSingleFqAdpter(  (P2In->InFq1) , path  , P2In   );
		}

	}

	if (FqCheck.size()<2)
	{
		cout <<"\twarning: only 1.fqchek , no Png Figure Out"<<endl;
		delete para_A24 ;
		delete P2In ;
		return (0);

	}

	FILE   *stream;
	char   buf[1024]={'\0'};
	string cc="which  gnuplot   2> /dev/null ";
	memset( buf, '\0', sizeof(buf) );
	stream=popen(cc.c_str(),"r") ;
	fread( buf, sizeof(char), sizeof(buf), stream);
	string binPath=buf;
	if (binPath == "")
	{
		cout <<"\twarning::can't find the [gnuplot] in your $PATH ; no Png Figure Out"<<endl;
	}
	else
	{
		string cc="which  convert  2> /dev/null ";
		memset( buf, '\0', sizeof(buf) );
		stream=popen(cc.c_str(),"r") ;
		fread( buf, sizeof(char), sizeof(buf), stream);
		binPath=buf;
		if (binPath == "" )
		{
			cout <<"\twarning:can't find the [convert] in your $PATH ; no png Figure Out"<<endl;
		}
		else
		{
			string path=(para_A24->OutFq1);
			if  (path.rfind('/') ==string::npos)
			{
				path="./";
			}
			else
			{
				path=path.substr(0, path.rfind('/') + 1);
			}
			PlotfqCheck (  FqCheck ,path );
		}
	}

	delete para_A24 ;
	delete P2In ;
	return (0);
}

#endif  // FQ_check_H_




