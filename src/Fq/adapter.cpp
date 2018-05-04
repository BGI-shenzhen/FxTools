#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <math.h>
#include <zlib.h>
#define _A 0
#define _C 1
#define _G 2
#define _T 3
#define _N 4
#define MAX_L 256
#define MAX_Q 128
#define Q_SHIFT -64
using namespace std;

void usage();

int main (int argc, char *argv[ ]) {

	int opt;
	int afrag,rfrag,lfrag,sfrag,cfrag;
	afrag=rfrag=lfrag=sfrag=cfrag=0;
	string adpt_file, read_file; //input file
	string list_file, stat_file; //output file
	string fqcheck_file; //output file
	while ((opt=getopt(argc, argv, "a:r:l:s:c:h"))!=-1) {
		switch (opt) {
			case 'a': afrag=1; adpt_file=optarg; break;    //input file
			case 'r': rfrag=1; read_file=optarg; break;    //input file
			case 'l': lfrag=1; list_file=optarg; break;    //output file
			case 's': sfrag=1; stat_file=optarg; break;    //output file
			case 'c': cfrag=1; fqcheck_file=optarg; break; //output file
			case 'h': usage(); break;
			default:  usage();
		}
	}
	if (!(afrag&&rfrag&&lfrag&&sfrag&&cfrag)) usage();

	static long int base[5];//记录每种碱基的总数
	static long int qual[MAX_Q];//每个质量值的总频率
	static long int cycle_qual[MAX_L][MAX_Q];//每个cycle的质量分布
	static long int cycle_base[MAX_L][5];//每个cycle的碱基分布
	static double   q_err[MAX_Q];//记录质量为q的碱基的错误率
	static double   err = 0.;
	static long int q20 = 0;
	static long int q30 = 0;
	int read_num = 0;
	int adpt_num = 0;
	int empt_num = 0;


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
	string adpt_id, adpt_se;
	int read_len;
	int adpt_len;
	gzFile read_fs;
	read_fs = gzopen (read_file.c_str(), "rb");
	fstream adpt_fs;
	adpt_fs.open (adpt_file.c_str(), ios::in);
	if (read_fs==NULL || adpt_fs==NULL) {
		cout<<"Cann't open files!"<<endl;
		usage();
	}
	gzgets(read_fs, read_id, MAX_L);
	gzgets(read_fs, read_se, MAX_L);
	string read = read_se;
	read_len = read.length()-1;
	gzrewind(read_fs);
	getline (adpt_fs, adpt_id, '\n');
	getline (adpt_fs, adpt_se, '\n');
	adpt_id.erase(0,1);
	adpt_len = adpt_se.length();
	adpt_fs.close();

	//output file
	gzFile list_fs;
	list_fs = gzopen (list_file.c_str(), "wb");
	fstream stat_fs;
	stat_fs.open (stat_file.c_str(), ios::out);
	fstream fqcheck_fs;
	fqcheck_fs.open (fqcheck_file.c_str(), ios::out);
	if (list_fs==NULL || stat_fs==NULL || fqcheck_fs==NULL) {
		cout<<"Cann't open files!"<<endl;
		usage();
	}
	string title = "#reads_id   reads_len   reads_start   reads_end   adapter_id   adapter_len   adapter_start   adapter_end   align_len   mismatch\n";
	gzwrite (list_fs, (char *)title.c_str(), title.size());
	stat_fs<<"adapter_id\tpolluted_reads\tempty_reads\tadapter_sequence"<<endl;

	int f = 0;
	std::ostringstream ss;
	while (gzgets (read_fs, read_id, MAX_L)) {
		gzgets (read_fs, read_se, MAX_L);
		gzgets (read_fs, qual_id, MAX_L);
		gzgets (read_fs, qual_se, MAX_L);
		++read_num;
		for (int c=0;c<read_len;c++) {
			switch (read_se[c]) {
				case 'A': ++cycle_base[c][_A];break;
				case 'C': ++cycle_base[c][_C];break;
				case 'G': ++cycle_base[c][_G];break;
				case 'T': ++cycle_base[c][_T];break;
				case 'N': ++cycle_base[c][_N];break;
			}
		}
		for (int c=0;c<read_len;c++) {
			int q = qual_se[c]+Q_SHIFT;
			if(q<0) {q=0;}
			++cycle_qual[c][q];
		}

		int find=0;
		int a1 = adpt_len-15;
		int r1 = 0;
		int len;
		int mis;
		for (r1=0; r1<=read_len-10;) {
			int len1 = adpt_len - a1;
			int len2 = read_len - r1;
			len = (len1<len2)? len1:len2;
			mis = 0;
			int map[MAX_L];
			map[0]=0;
			for (int c=0; c<len; c++) {
				if (adpt_se[a1+c]==read_se[r1+c]) {map[mis]++;}
				else {mis++;map[mis]=0;}
			}
			int max_map=0;
			for (int c=0; c<=mis; c++) {
				if (map[c]>max_map) {max_map=map[c];}
			}
			if ((len*0.2>mis) or (max_map>=15)) {
				find = 1;
				break;
			}
			if (a1>0) {a1--;}
			else      {r1++;}
		}
		if (find) {
			string id = read_id;
			id.erase(0,1);
			id.erase(id.length()-1,1);
			ss<<id<<'\t'<<read_len<<'\t'<<r1<<'\t'<<r1+len-1<<'\t'
				<<adpt_id<<'\t'<<adpt_len<<'\t'<<a1<<'\t'<<a1+len-1<<'\t'
				<<len<<'\t'<<mis<<endl;
			++adpt_num;
			if (r1<=3) {++empt_num;}
			f++;
			if (f==500) {
				f=0;
				std::string sss = ss.str();
				gzwrite(list_fs, (char *)sss.c_str(), sss.size());
				ss.str("");
			}
		}
	}
	gzclose(read_fs);

	std::string sss = ss.str();
	gzwrite(list_fs, (char *)sss.c_str(), sss.size());
	ss.str("");
	gzclose(list_fs);

	float adpt_percent = ((float) adpt_num)/read_num*100;
	float empt_percent = ((float) empt_num)/read_num*100;
	stat_fs<<adpt_id<<'\t'<<adpt_num<<" ("<<setprecision(3)<<fixed<<adpt_percent<<"%)"<<'\t'<<empt_num<<" ("<<empt_percent<<"%)"<<'\t'<<adpt_se<<endl;
	stat_fs<<"\ntotal_reads: "<<read_num<<endl;
	stat_fs.close();

	long base_total = (long) read_len * read_num;
	for(int c=0;c<read_len;c++){
		base[_A]+=cycle_base[c][_A];
		base[_C]+=cycle_base[c][_C];
		base[_G]+=cycle_base[c][_G];
		base[_T]+=cycle_base[c][_T];
		base[_N]+=cycle_base[c][_N];
	}
	for(int c=0;c<read_len;c++){
		for(int q=0;q<=MAX_Q;q++) {
			qual[q]+=cycle_qual[c][q];
		}
	}
	int qMAX=0;
	for(int q=0;q<MAX_Q;q++) {
		if (qual[q]) {
			qMAX=q;
			if (q>=30) q30+=qual[q];
			if (q>=20) q20+=qual[q];
			err+=q_err[q]*qual[q];
		}
	}

	fqcheck_fs<<" the default quality shift value is: "<<Q_SHIFT<<", "
		<<read_num<<" sequences, "
		<<base_total<<" total length, Max length:"
		<<read_len<<", average length:"
		<<setiosflags(ios::fixed)<<setprecision(2)<<(double)read_len<<endl;
	fqcheck_fs<<"Standard deviations at 0.25:  total "<<100*(sqrt(0.25*(double)base_total)/base_total)
		<<"%, per base "<<100*(sqrt(0.25*(double)read_num)/read_num)
		<<"%"<<endl;
	fqcheck_fs<<"             A     C     G     T     N ";
	for(int q=0;q<=qMAX;q++) {fqcheck_fs<<setw(4)<<q<<' ';}
	fqcheck_fs<<endl;
	fqcheck_fs<<"Total    ";
	for(int b=_A;b<=_N;b++)  {fqcheck_fs<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<100*(double)base[b]/base_total<<' ';}
	for(int q=0;q<=qMAX;q++) {fqcheck_fs<<setw(4)<<(int)(1000*((double)qual[q]/base_total))<<" ";}
	fqcheck_fs<<endl;
	for(int c=0;c<read_len;c++)	{
		fqcheck_fs<<"base "<<setw(3)<<c+1<<' ';
		for(int b=_A;b<=_N;b++)  {fqcheck_fs<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<100*(double)cycle_base[c][b]/read_num<<' ';}
		for(int q=0;q<=qMAX;q++) {fqcheck_fs<<setw(4)<<(int)(1000*((double)cycle_qual[c][q]/read_num))<<' ';}
		fqcheck_fs<<endl;
	}
	fqcheck_fs<<endl;
	fqcheck_fs<<"Error Rate\t%GC\tQ20\tQ30"<<endl;
	fqcheck_fs.precision(2);
	fqcheck_fs<<(100.0*err)/(base_total)<<"\t"<<(100.0*(base[_C]+base[_G]))/(base_total-base[_N])<<"\t"<<(100.0*q20)/base_total<<"\t"<<(100.0*q30)/base_total<<endl;
	fqcheck_fs.close();
	return(0);

}

void usage() {
	cout << "\nUsage: filter_adapter -a <*.fa> -r <*.fq.gz> -l <*.adapter.list> -s <*.adapter.stat> -c <*.fqcheck>\n"
		<< "  -a <str>   input fasta file of adapters\n"
		<< "  -r <str>   input fastq file of reads\n"
		<< "  -l <str>   output adapter list file\n"
		<< "  -s <str>   output adapter statistics file\n"
		<< "  -c <str>   output fqcheck file\n"
		<< endl ;
	exit(1);
}
