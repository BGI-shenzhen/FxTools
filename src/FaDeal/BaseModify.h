#ifndef FQ_BaseModify_H_
#define FQ_BaseModify_H_
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
typedef long long  llong ;

int  usage_sort_FC18()
{
	cout <<""
		"\n"
		"\tUsage: BaseModify -i <in.fa> -o <out.fa>  -s <A.snp>\n"
		"\n"
		"\t\t-i   <str>   Input fa File to Modify\n"
		"\t\t-s   <str>   Input File with chang Site[chr 100 N G]\n"
		"\t\t-o   <str>   Output Modify fa File\n"
		"\n"
		"\t\t-h           Show this help\n"
		"\n";
	return 1;
}

int parse_Acmd_FC18(int argc, char **argv ,In3str1v * para_FC18 )
{
	if (argc <=2 ) {usage_sort_FC18();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InPut"  || flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_FC18->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut" || flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FC18->InStr2=argv[i];
		}       
		else if (flag  ==  "Site" || flag  == "s")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FC18->InStr3=argv[i];
		} 
		else if (flag  == "help" || flag  == "h")
		{
			usage_sort_FC18();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para_FC18->InStr2).empty() || (para_FC18->InStr1).empty() || (para_FC18->InStr3).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	return 1 ;
}



//int main(int argc, char *argv[])
int FA_BaseModify_main(int argc, char *argv[])
{
	In3str1v * para_FC18 = new In3str1v;

	if(parse_Acmd_FC18(argc, argv, para_FC18 )==0)
	{
		delete  para_FC18 ;
		return 1;    
	}

	int linecut= FaCutLine ((para_FC18->InStr1));

	(para_FC18->InStr2)=add_Asuffix((para_FC18->InStr2));
	ogzstream  OUT ((para_FC18->InStr2).c_str());
	if(!OUT.good())
	{
		cerr << "open OUT File error: "<<(para_FC18->InStr2)<<endl;
		delete para_FC18 ;return 1;
	}

	igzstream SNP ((para_FC18->InStr3).c_str(),ifstream::in);
	if (SNP.fail())
	{
		cerr << "open SNP File error: "<<(para_FC18->InStr3)<<endl;
		delete para_FC18  ; return  1;
	}

	map <string,map <llong ,char> >  Site;

	while(!SNP.eof())
	{
		string  line ;
		getline(SNP,line);
		if (line.length()<=0)  { continue ; }
		istringstream isone (line,istringstream::in);
		string chr ,depth,cp ;

		llong position ;
		char ref  ,allele1 ;
		isone>>chr>>position>>ref>>allele1;

		map  <string,map <llong ,char> >  :: iterator it=Site.find(chr);
		if (it != Site.end())
		{
			(it->second).insert(map <llong,char>  :: value_type(position,allele1)) ;
		}
		else
		{
			map <llong,char > gene_cds_str;
			gene_cds_str[position]=allele1;
			Site[chr]=gene_cds_str;
		}
	}
	SNP.close();

	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen((para_FC18->InStr1).c_str(),"r");
	seq = kseq_init(fp);

	while ((l = kseq_read(seq)) >= 0 )
	{
		string RefSeq=(seq->seq.s);
		string chr=(seq->name.s);

		map < string,map <llong,char> > :: iterator it=Site.find(chr);

		if (it != Site.end())
		{

			map <llong, char > :: iterator Y ;
			for(  Y=(it->second).begin() ; Y!= (it->second).end(); Y++ )
			{
				RefSeq[(Y->first)-1]=Y->second ;
			}
		}
		string ID=chr;
		if (seq->comment.l)
		{
			string comment=(seq->comment.s);
			ID+="\t"+comment;
		}
		Display( RefSeq ,ID , OUT , linecut);
	}

	kseq_destroy(seq);
	gzclose(fp);
	OUT.close();

	delete para_FC18 ;
	return 0;
}

#endif // BaseModify

///////// swimming in the sky and flying in the sea ////////////
