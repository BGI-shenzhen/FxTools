#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <map>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <gzstream.h>
#include "../ALL/comm.h"

#define MAX_LEN_FILE 128
#define MAXFILE 256
#define VERSION "1.02"
using namespace std;

string Basename(string filepath){
	char basename[MAX_LEN_FILE]; 
	if(filepath.size() == 0){
		return NULL;
	}else{
		size_t index = filepath.find_last_of("/");
		if(index == string::npos){
			return NULL;
		}else{
			int i = 0;
			for(i;filepath[index+i+1]!='\0'; i++)	basename[i]=filepath[index+i+1];
			basename[i]='\0';
			return basename;
		}
	}
}

string Basename(string filepath, int k){
	if(filepath.size() == 0){
		return NULL;
	}else{
		size_t index = filepath.find_last_of("/");
		if(index == string::npos){
			return NULL;
		}else{
			string name = filepath.substr(index+1, filepath.size());
			index = name.find(".tr.");
			if(index == string::npos){
				return name;
			}else{
				string last_name = name.substr(0,index);
				return last_name+".gz";
			}
		}
	}
}
inline string GetId(string head){
	int i = 0;
	size_t len = head.size();
	for(i = 0; i<len;i++){
		if(head[i] == '_'){
			break;
		}
	}
	return head.substr(i+1, len-i-3);
}

/*
   int Usage(){
   cerr << "\n";
   cerr << "DataFilter: a program for filtering pooling and index fastq data" << endl;
   cerr << "Version:" << VERSION << endl;
   cerr << "Usage:\t" << "./DataFilter <command> [options]" << endl;
   cerr << "Command: \tpooling     pooling library data filter" << endl;
   cerr << "         \tindex       index library data filter"   << endl;
   cerr << "         \tSE          SE index library data filter"<< endl;
   cerr << "Please choose the right type of filter for you project\n" << endl;
   return 1;
   }
   */

int pooling_usage(){
	cerr << "\nUsage:" << endl;
	cerr << "pooling <1.fq.tr.gz> <2.fq.tr.gz> <1.adapter.list.gz> <2.adapter.list.gz> <outdir> <sample_list_of_pooling> <needrawdata 1:0>\n" << endl;
	return 1;
}

int index_usage(){
	cerr << "\nUsage:" << endl;
	cerr << "rmAdapter <1.fq.gz> <2.fq.gz> <1.adapter.list.gz> <2.adapter.list.gz> <outdir>\n" << endl;
	return 1;
}

int se_usage(){
	cerr << "\nUsage:" << endl;
	cerr << "rmAdapterSE <1.fq.gz> <1.adapter.list.gz> <outdir>\n" << endl;
}

bool FQ_Pooling_main(int argc, char *argv[]){
	//bool Pooling(int argc, char *argv[]){
	if(argc < 7){
		pooling_usage();
		return false;
	}
	string fq1 = argv[1];
	string fq2 = argv[2];
	string adp1 = argv[3];
	string adp2 = argv[4];
	string od = argv[5];
	string list = argv[6];
	int tempshiftQ=GetShiftQ( fq1 );
	cout<<"ShiftQ is "<<tempshiftQ<<endl;
	//	string index_id = argv[7];
	int yes = 0;
	//	if(argc == 9){
	//		yes = atoi(argv[8]);
	//	}
	if(argc == 8){
		yes = atoi(argv[7]);
	}
	igzstream gz_fq1(fq1.c_str());
	igzstream gz_fq2(fq2.c_str());
	igzstream gz_ad1(adp1.c_str());
	igzstream gz_ad2(adp2.c_str());
	ifstream lt(list.c_str());
	string line("");
	vector<string> line1(4);
	vector<string> line2(4);
	map<string,int> dap1;
	map<string,int> dap2;
	map<string,string> LIST;
	vector<string> pool_index;
	map<string,long> raw_gc;
	map<string,long> raw_q20;
	map<string,long> raw_read;
	map<string,long> raw_bases;
	map<string,long> clean_gc;
	map<string,long> clean_q20;
	map<string,long> clean_read;
	map<string,long> clean_bases;
	getline(gz_ad1,line);
	while(getline(gz_ad1,line)){
		istringstream isline(line);
		string head;
		int i =0;
		isline >> head;
		for(i=0;i<head.size();i++){
			if(head[i] == '#'  ||  head[i] == ' ') break;
		}
		dap1[head.substr(0,i)]=1;
	}
	getline(gz_ad2,line);
	while(getline(gz_ad2,line)){
		istringstream isline(line);
		string head;
		int i=0;
		isline >> head;
		for(i=0;i<head.size();i++){
			if(head[i]=='#'  ||  head[i] == ' ' ) break;
		}
		dap2[head.substr(0,i)]=1;
	}
	ogzstream handle[MAXFILE];
	ogzstream hhandle[MAXFILE];
	ogzstream raw_handle[MAXFILE];
	ogzstream raw_hhandle[MAXFILE];
	int jk = 0;
	while(getline(lt,line)){
		istringstream isline(line);
		string index,sample;
		isline >> sample >> index;
		//	index=index_id+"_"+index;
		pool_index.push_back(index);
		LIST[index] = sample;
		raw_gc[sample]=0;
		raw_q20[sample]=0;
		raw_read[sample] = 0;
		raw_bases[sample] = 0;
		clean_gc[sample]=0;
		clean_q20[sample]=0;
		clean_bases[sample] = 0;
		clean_read[sample]=0;
		string path = od + "/" + sample;
		if(access(path.c_str(),F_OK) != 0){
			mkdir(path.c_str(),0755);
		}

		string temp = path+"/"+Basename(fq1, 0);
		handle[jk].open(temp.c_str());

		temp = path+"/"+Basename(fq2, 0);
		hhandle[jk].open(temp.c_str());

		string raw_path = od+"/raw/"+sample;
		string rp = od+"/raw";
		if(yes == 1){
			if(access(rp.c_str(),F_OK)!=0){
				mkdir(rp.c_str(),0755);
			}
			if(access(raw_path.c_str(),F_OK) !=0 ){
				mkdir(raw_path.c_str(),0755);
			}
			string ttp = raw_path+"/"+Basename(fq1, 0);
			raw_handle[jk].open(ttp.c_str());
			ttp=raw_path+"/"+Basename(fq2, 0);
			raw_hhandle[jk].open(ttp.c_str());
		}

		jk++;
	}
	while(getline(gz_fq1,line1[0])){
		getline(gz_fq2,line2[0]);

		getline(gz_fq1,line1[1]);
		getline(gz_fq2,line2[1]);

		getline(gz_fq1,line1[2]);
		getline(gz_fq2,line2[2]);

		getline(gz_fq1,line1[3]);
		getline(gz_fq2,line2[3]);

		string id1 = GetId(line1[0]);
		string id2 = GetId(line2[0]);
		int length = line1[3].size();
		int ind=-1,kkj=0;
		vector<string>::iterator itt = pool_index.begin();
		for(itt;itt != pool_index.end();itt++){
			if(*itt==id1 && id1==id2){
				ind=kkj;
				break;
			}
			kkj++;
		}
		if(ind == -1){
			continue;
		}
		raw_read[pool_index[ind]]+=2;
		raw_bases[pool_index[ind]]+=length*2;


		if(yes == 1){
			raw_handle[ind] << line1[0] << "\n" << line1[1] << "\n" << line1[2] << "\n" << line1[3] << endl;
			raw_hhandle[ind] << line2[0] << "\n" << line2[1] << "\n" << line2[2] << "\n" << line2[3] << endl;
		}

		int gc = 0;
		int qual = 0;
		int q20 = 0;
		int n_Num=0;
		for(int j = 0; j<length;j++){
			if(line1[1][j] == 'G' || line1[1][j] == 'C' || line1[1][j] == 'g' || line1[1][j] == 'c'){
				gc++;
			}
			else if (line1[1][j] == 'N' || line1[1][j] == 'n' )
			{
				n_Num++;
			}
			int q = line1[3][j]-tempshiftQ;
			if(q <= 10){
				qual++;
			}
			if(q >= 20){
				q20++;
			}
		}

		int ggc=0;
		int qqual=0;
		int qq20 = 0;
		int qn_Num=0;
		for(int j = 0; j<length;j++){
			if(line2[1][j] == 'G' || line2[1][j] == 'C' || line2[1][j] == 'g' || line2[1][j] == 'c'){
				ggc++;
			}
			else if ( line2[1][j] == 'N' || line2[1][j] == 'n')
			{
				qn_Num++;
			}
			int q = line2[3][j]-tempshiftQ;
			if(q <= 10){
				qqual++;
			}
			if(q >= 20){
				qq20++;
			}
		}
		raw_gc[pool_index[ind]]=raw_gc[pool_index[ind]]+gc+ggc;
		raw_q20[pool_index[ind]]=raw_q20[pool_index[ind]]+q20+qq20;
		if( (qqual >= (length*0.3)) || (qual >= (length*0.3)) ||  (n_Num*20> length) ||  (qn_Num*20> length) ){
			continue;
		}

		char fqtag1[256];
		for(int j=1; (line1[0][j]!='#' &&  line1[0][j]!=' ');j++){
			fqtag1[j-1] = line1[0][j];
			fqtag1[j]='\0';
		}
		char fqtag2[256];
		for(int j=1;(line2[0][j]!='#' && line2[0][j]!=' ' );j++){
			fqtag2[j-1] = line2[0][j];
			fqtag2[j]='\0';
		}
		map<string,int>::iterator adit1 = dap1.find(fqtag1);
		map<string,int>::iterator adit2 = dap2.find(fqtag2);
		if(adit1 != dap1.end()||adit2 != dap2.end()){
			continue;
		}

		clean_gc[pool_index[ind]]=clean_gc[pool_index[ind]]+gc+ggc;
		clean_q20[pool_index[ind]]=clean_q20[pool_index[ind]]+q20+qq20;
		clean_read[pool_index[ind]]+=2;
		clean_bases[pool_index[ind]]+=length*2;
		handle[ind] << line1[0] << "\n" << line1[1] << "\n" << line1[2] << "\n" << line1[3] << endl;
		hhandle[ind] << line2[0] << "\n" << line2[1] << "\n" << line2[2] << "\n" << line2[3] << endl;
	}

	for(int ko=0;ko<pool_index.size();ko++){
		ofstream raw_info_file;
		ofstream clean_info_file;
		string nameR = od+"/"+LIST[pool_index[ko]]+"/raw_"+Basename(fq1, 0)+".total.info";
		raw_info_file.open(nameR.c_str());
		string nameC = od+"/"+LIST[pool_index[ko]]+"/clean_"+Basename(fq1, 0)+".total.info";
		clean_info_file.open(nameC.c_str());

		float raw_q20_rate = (float)raw_q20[pool_index[ko]]/raw_bases[pool_index[ko]]*100;
		float raw_gc_rate = (float)raw_gc[pool_index[ko]]/raw_bases[pool_index[ko]]*100;
		float clean_q20_rate = (float)clean_q20[pool_index[ko]]/clean_bases[pool_index[ko]]*100;
		float clean_gc_rate = (float)clean_gc[pool_index[ko]]/clean_bases[pool_index[ko]]*100;
		raw_info_file << "fqfile\t" << "Reads(M)\t" << "Bases(G)\t" << "GC_Rate(%)\t" << "Q20_Rate(%)" << endl;
		raw_info_file << Basename(fq1, 0) << "\t" << raw_read[pool_index[ko]]/1000000.0 << "\t" << raw_bases[pool_index[ko]]/1000000000.0 << "\t" << raw_gc_rate << "\t" << raw_q20_rate << endl;

		clean_info_file << "fqfile\t" << "Reads(M)\t" << "Bases(G)\t" << "GC_Rate(%)\t" << "Q20_Rate(%)" << endl;
		clean_info_file << Basename(fq1, 0) << "\t" << clean_read[pool_index[ko]]/1000000.0 << "\t" << clean_bases[pool_index[ko]]/1000000000.0 << "\t" << clean_gc_rate << "\t" << clean_q20_rate << endl;
		raw_info_file.close();
		clean_info_file.close();
		handle[ko].close();
		hhandle[ko].close();
		if(yes==1){
			raw_handle[ko].close();
			raw_hhandle[ko].close();
		}
	}

	return true;
}

bool FQ_Se_main(int argc, char *argv[]){
	if(argc != 4){
		se_usage();
		return false;
	}
	string fq1 = argv[1];
	int tempshiftQ=GetShiftQ( fq1 );
	cout<<"ShiftQ is "<<tempshiftQ<<endl;
	string adp1 = argv[2];
	string od = argv[3];
	igzstream gz_fq1(fq1.c_str());
	igzstream gz_ad1(adp1.c_str());
	string outname1 = od+"/"+Basename(fq1);
	ogzstream handle(outname1.c_str());
	string line("");
	vector<string> line1(4);
	map<string,int> dap1;
	long raw_gc=0;
	long raw_q20=0;
	long raw_read=0;
	long raw_bases=0;
	long clean_gc=0;
	long clean_q20=0;
	long clean_read=0;
	long clean_bases=0;
	getline(gz_ad1,line);
	while(getline(gz_ad1,line)){
		istringstream isline(line);
		string head;
		int i;
		isline >> head;
		for(i=0;i<head.size();i++){
			if(head[i]=='#' ||  head[i]==' ' ) break;
		}
		dap1[head.substr(0,i)]=1;
	}
	while(getline(gz_fq1,line1[0])){

		getline(gz_fq1,line1[1]);

		getline(gz_fq1,line1[2]);

		getline(gz_fq1,line1[3]);
		int read_len = line1[1].size();


		int gc = 0;
		int qual = 0;
		int q20 = 0;
		int N_num=0;
		for(int j = 0; j<read_len;j++){
			if(line1[1][j] == 'G' || line1[1][j] == 'C' || line1[1][j] == 'g' || line1[1][j] == 'c'){
				gc++;
			}
			else if ( (line1[1][j] == 'N')  ||  ( (line1[1][j] == 'n')) )
			{
				N_num++;
			}
			int q = line1[3][j]-tempshiftQ;
			if(q <=10){
				qual++;
			}
			if(q >= 20){
				q20++;
			}
		}

		raw_gc=raw_gc+gc;
		raw_q20=raw_q20+q20;
		raw_read+=1;
		raw_bases+=read_len;

		if( (qual >= read_len/3)  || ( N_num*20>  read_len )){
			continue;
		}
		char fqtag1[256];
		for(int j=1;(line1[0][j]!='#'  &&  line1[0][j]!=' ' ) ;j++){
			fqtag1[j-1] = line1[0][j];
			fqtag1[j]='\0';
		}
		map<string,int>::iterator adit1 = dap1.find(fqtag1);
		if(adit1 != dap1.end()){
			continue;
		}
		clean_gc=clean_gc+gc;
		clean_q20=clean_q20+q20;
		clean_read+=1;
		clean_bases+=read_len;
		handle << line1[0] << "\n" << line1[1] << "\n" << line1[2] << "\n" << line1[3] << endl;
	}
	ofstream raw_info_file;
	ofstream clean_info_file;
	string nameR = od+"/raw_"+Basename(fq1)+".total.info";
	raw_info_file.open(nameR.c_str());
	string nameC = od+"/clean_"+Basename(fq1)+".total.info";
	clean_info_file.open(nameC.c_str());

	float raw_q20_rate = (float)raw_q20/raw_bases*100;
	float raw_gc_rate = (float)raw_gc/raw_bases*100;
	float clean_q20_rate = (float)clean_q20/clean_bases*100;
	float clean_gc_rate = (float)clean_gc/clean_bases*100;

	raw_info_file << "fqfile\t" << "Reads\t" << "Bases\t" << "GC_Rate(%)\t" << "Q20_Rate(%)" << endl;
	raw_info_file << Basename(fq1) << "\t" << raw_read/1000000.0 << "\t" << raw_bases/1000000000.0 << "\t" << raw_gc_rate << "\t" << raw_q20_rate << endl;

	clean_info_file << "fqfile\t" << "Reads\t" << "Bases\t" << "GC_Rate(%)\t" << "Q20_Rate(%)" << endl;
	clean_info_file << Basename(fq1) << "\t" << clean_read/1000000.0 << "\t" << clean_bases/1000000000.0 << "\t" << clean_gc_rate << "\t" << clean_q20_rate << endl;
	raw_info_file.close();
	clean_info_file.close();
	handle.close();
	return true;
}

bool FQ_Index_main(int argc, char *argv[]) {
//int  main(int argc, char *argv[]) {
	//bool Index(int argc, char *argv[]){
	if(argc != 6){
		//	cout << "Write index usage" << endl;
		index_usage();
		return false;
	}
	string fq1 = argv[1];
	int tempshiftQ=GetShiftQ( fq1 );
	cout<<"ShiftQ is "<<tempshiftQ<<endl;
	string fq2 = argv[2];
	string adp1 = argv[3];
	string adp2 = argv[4];
	string od = argv[5];
	igzstream gz_fq1(fq1.c_str());
	igzstream gz_fq2(fq2.c_str());
	igzstream gz_ad1(adp1.c_str());
	igzstream gz_ad2(adp2.c_str());
	string outname1 = od+"/"+Basename(fq1);
	string outname2 = od+"/"+Basename(fq2);
	ogzstream handle(outname1.c_str());
	ogzstream hhandle(outname2.c_str());
	string line("");
	vector<string> line1(4);
	vector<string> line2(4);
	map<string,int> dap1;
	map<string,int> dap2;
	long raw_gc=0;
	long raw_q20=0;
	long raw_read=0;
	long raw_bases=0;
	long clean_gc=0;
	long clean_q20=0;
	long clean_read=0;
	long clean_bases=0;
	getline(gz_ad1,line);
	while(getline(gz_ad1,line)){
		istringstream isline(line);
		string head;
		int i;
		isline >> head;
		for(i=0;i<head.size();i++){
			if(head[i]=='#' ||  head[i]==' ') break;
		}
		dap1[head.substr(0,i)]=1;
	}
	getline(gz_ad2,line);
	while(getline(gz_ad2,line)){
		istringstream isline(line);
		string head;
		int i;
		isline >> head;
		for(i=0;i<head.size();i++){
			if(head[i]=='#' ||  head[i]==' ' ) break;
		}
		dap2[head.substr(0,i)]=1;
	}
	while(getline(gz_fq1,line1[0])){
		getline(gz_fq2,line2[0]);

		getline(gz_fq1,line1[1]);
		getline(gz_fq2,line2[1]);

		getline(gz_fq1,line1[2]);
		getline(gz_fq2,line2[2]);

		getline(gz_fq1,line1[3]);
		getline(gz_fq2,line2[3]);
		int read_len = line1[1].size();


		int gc = 0;
		int qual = 0;
		int q20 = 0;
		int n_Num =0;
		for(int j = 0; j<read_len;j++){
			if(line1[1][j] == 'G' || line1[1][j] == 'C' || line1[1][j] == 'g' || line1[1][j] == 'c'){
				gc++;
			}
			else if ( line1[1][j] == 'N' || line1[1][j] == 'n')
			{
				n_Num++;
			}
			int q = line1[3][j]-tempshiftQ;
			if(q <= 10){
				qual++;
			}
			if(q >= 20){
				q20++;
			}
		}

		int ggc=0;
		int qqual=0;
		int qq20 = 0;
		int qn_Num = 0;
		for(int j = 0; j<read_len;j++){
			if(line2[1][j] == 'G' || line2[1][j] == 'C' || line2[1][j] == 'g' || line2[1][j] == 'c'){
				ggc++;
			}
			else if ( line2[1][j] == 'N' || line2[1][j] == 'n')
			{
				qn_Num++;
			}

			int q = line2[3][j]-tempshiftQ;
			if(q <= 10){
				qqual++;
			}
			if(q >= 20){
				qq20++;
			}
		}

		raw_gc=raw_gc+gc+ggc;
		raw_q20=raw_q20+q20+qq20;
		raw_read+=2;
		raw_bases+=read_len*2;

		if( (qqual >= (read_len*0.3)) || (qual >= (read_len*0.3))  || (qn_Num*20> read_len) ||  (n_Num*20> read_len ) ){
			continue;
		}
		char fqtag1[256];
		for(int j=1; (line1[0][j]!='#'  &&  line1[0][j]!=' ') ;j++){
			fqtag1[j-1] = line1[0][j];
			fqtag1[j]='\0';
		}
		char fqtag2[256];
		for(int j=1;(line2[0][j]!='#' &&  line2[0][j]!=' ' );j++){
			fqtag2[j-1] = line2[0][j];
			fqtag2[j]='\0';
		}
		map<string,int>::iterator adit1 = dap1.find(fqtag1);
		map<string,int>::iterator adit2 = dap2.find(fqtag2);
		if(adit1 != dap1.end()||adit2 != dap2.end()){
			continue;
		}
		clean_gc=clean_gc+gc+ggc;
		clean_q20=clean_q20+q20+qq20;
		clean_read+=2;
		clean_bases+=read_len*2;
		handle << line1[0] << "\n" << line1[1] << "\n" << line1[2] << "\n" << line1[3] << endl;
		hhandle << line2[0] << "\n" << line2[1] << "\n" << line2[2] << "\n" << line2[3] << endl;
	}
	ofstream raw_info_file;
	ofstream clean_info_file;
	string nameR = od+"/raw_"+Basename(fq1)+".total.info";
	raw_info_file.open(nameR.c_str());
	string nameC = od+"/clean_"+Basename(fq1)+".total.info";
	clean_info_file.open(nameC.c_str());

	float raw_q20_rate = (float)raw_q20/raw_bases*100;
	float raw_gc_rate = (float)raw_gc/raw_bases*100;
	float clean_q20_rate = (float)clean_q20/clean_bases*100;
	float clean_gc_rate = (float)clean_gc/clean_bases*100;

	raw_info_file << "fqfile\t" << "Reads\t" << "Bases\t" << "GC_Rate(%)\t" << "Q20_Rate(%)" << endl;
	raw_info_file << Basename(fq1) << "\t" << raw_read/1000000.0 << "\t" << raw_bases/1000000000.0 << "\t" << raw_gc_rate << "\t" << raw_q20_rate << endl;

	clean_info_file << "fqfile\t" << "Reads\t" << "Bases\t" << "GC_Rate(%)\t" << "Q20_Rate(%)" << endl;
	clean_info_file << Basename(fq1) << "\t" << clean_read/1000000.0 << "\t" << clean_bases/1000000000.0 << "\t" << clean_gc_rate << "\t" << clean_q20_rate << endl;
	raw_info_file.close();
	clean_info_file.close();
	handle.close();
	hhandle.close();
	return true;
}

/*////
  int main(int argc, char *argv[]){
  if(argc < 2){
  return Usage();
  }
  if(strcmp(argv[1],"pooling") == 0){
  return Pooling(argc-1, argv+1);
  }else if(strcmp(argv[1],"index") == 0){
  return Index(argc-1,argv+1);
  }else if(strcmp(argv[1],"SE") == 0){
  return Se(argc-1,argv+1);
  }
  else{
  return Usage();
  }
  return 0;
  }
///*////

