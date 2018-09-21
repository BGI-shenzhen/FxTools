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
#include "../ALL/gzstream.h"
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
		"\tUsage: BaseModify -i <in.fa> -o <out.fa>  -s <modify.list>\n"
		"\n"
		"\t\t-i   <str>   input FASTA\n"
		"\t\t-o   <str>   output file\n"
		"\t\t-s   <str>   file containing sites to be modified\n"
		"\t\t             formatted as [seq_id  site  original_base  modification_base]\n"
		"\n"
		"\t\t-h           show more details for help\n"
		"\n";
	return 1;
}

void More_HelpFA_13()
{
	cout<<"\n\
\n\
\t\t1.  BaseModify -i <in.fa> -o AAA -s modify_list\n\
\t\tthis will modify the bases in modify_list in input FASTA and output the result to a compressed file named AAA in current directory. For example, if we would like to modify two bases, one is on site 10 in seq1 from A to G, and the other is on site 20 in seq2 from C to T, the modify_list would show as:\n\
\t\tseq1 10 A G\n\
\t\tseq2 20 C T\n\
\n\
\n";
}


int parse_Acmd_FC18(int argc, char **argv ,In3str1v * para_FC18 )
{
	if (argc <=1 ) {usage_sort_FC18();return 0;}

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
			More_HelpFA_13();return 0;
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
