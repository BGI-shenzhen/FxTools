#ifndef CDS2Pep_H_
#define CDS2Pep_H_
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

int  print_AusageA09()
{
	cout <<""
		"\n"
		"\tUsage: CDS2Pep  -i <inCDS.fa>  -o <outPep.fa>\n"
		"\n"
		"\t\t-i    <str>   Input CDS fa \n"
		"\t\t-o    <str>   Output pep file\n"
		"\n"
		"\t\t-w            Warning Out the wrong CDS SeqID\n"
		"\t\t-h            show this help\n" 
		"\n";
	return 1;
}


int parse_AcmdA09(int argc, char **argv , In3str1v * paraA09 )
{
	if (argc <=2 ) {print_AusageA09();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InCDS" ||  flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraA09->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut"||  flag  == "o" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			paraA09->InStr2=argv[i];
		}
		else if (flag  == "Warning" ||  flag  == "w")
		{
			paraA09->TF=false ;
		}
		else if (flag  == "help"||  flag  == "h")
		{
			print_AusageA09();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((paraA09->InStr1).empty()  ||  (paraA09->InStr2).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	(paraA09->InStr2)=add_Asuffix(paraA09->InStr2) ;

	return 1 ;
}

int  CDS2pep  ( map <string,string> CODE , string inStr , string & OutStr )
{
	int length=inStr.length();
	string temp_vv ;
	for( int ii=0 ; ii<length; ii+=3)
	{
		string base=inStr.substr(ii,3);	
		temp_vv=CODE[base];
		if (temp_vv == "")
		{
			temp_vv="X";
		}
		OutStr+=temp_vv ;
	}

	if( (temp_vv != "U")  &&  (temp_vv != "*") )
	{
		if ((length)%3!=0) 
		{
			return 0 ;
		}
		else if (temp_vv ==  "")
		{
			return -1 ;
		}
		else
		{
			return -2 ;
		}
	}
	return 1 ;
}


int FA_CDS2Pep_main(int argc, char *argv[])
{
	In3str1v * paraA09 = new In3str1v;
	if( parse_AcmdA09(argc, argv, paraA09 )==0)
	{
		delete  paraA09 ;
		return 1;
	}

	map <string,string> CODE ;

	CODE["GCA"]="A"; CODE["GCC"]="A"; CODE["GCG"]="A"; CODE["GCT"]="A";  // #Alanine
	CODE["TGC"]="C"; CODE["TGT"]="C";                                    // #Cysteine
	CODE["GAC"]="D"; CODE["GAT"]="D";                                    // #AsparticAcid
	CODE["GAA"]="E"; CODE["GAG"]="E";                                    // #GlutamicAcid
	CODE["TTC"]="F"; CODE["TTT"]="F";                                    // #Phenylalanine
	CODE["GGA"]="G"; CODE["GGC"]="G"; CODE["GGG"]="G"; CODE["GGT"]="G";  // #Glycine
	CODE["CAC"]="H"; CODE["CAT"]="H";                                    // #Histidine
	CODE["ATA"]="I"; CODE["ATC"]="I"; CODE["ATT"]="I";                   // #Isoleucine
	CODE["AAA"]="K"; CODE["AAG"]="K";                                    // #Lysine
	CODE["CTA"]="L"; CODE["CTC"]="L"; CODE["CTG"]="L"; CODE["CTT"]="L";
	CODE["TTA"]="L"; CODE["TTG"]="L";                                    // #Leucine
	CODE["ATG"]="M";                                                     // #Methionine
	CODE["AAC"]="N"; CODE["AAT"]="N";                                    // #Asparagine
	CODE["CCA"]="P"; CODE["CCC"]="P"; CODE["CCG"]="P"; CODE["CCT"]="P";  // #Proline
	CODE["CAA"]="Q"; CODE["CAG"]="Q";                                    // #Glutamine
	CODE["CGA"]="R"; CODE["CGC"]="R"; CODE["CGG"]="R"; CODE["CGT"]="R";
	CODE["AGA"]="R"; CODE["AGG"]="R";                                    // #Arginine
	CODE["TCA"]="S"; CODE["TCC"]="S"; CODE["TCG"]="S"; CODE["TCT"]="S";
	CODE["AGC"]="S"; CODE["AGT"]="S";                                    // #Serine
	CODE["ACA"]="T"; CODE["ACC"]="T"; CODE["ACG"]="T"; CODE["ACT"]="T";  // #Threonine
	CODE["GTA"]="V"; CODE["GTC"]="V"; CODE["GTG"]="V"; CODE["GTT"]="V";  // #Valine
	CODE["TGG"]="W";                                                     // #Tryptophan
	CODE["TAC"]="Y"; CODE["TAT"]="Y";                                    // #Tyrosine
	CODE["TAA"]="*"; CODE["TAG"]="*"; CODE["TGA"]="*";                   // #Stop
	CODE["NNN"]="X";                                                     // #NA_NA

	int linecut= FaCutLine ((paraA09->InStr1));

	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen((paraA09->InStr1).c_str(), "r");
	seq = kseq_init(fp);
	ogzstream  OUT ((paraA09->InStr2).c_str());
	if(!OUT.good())
	{
		cerr << "open OUT File error: "<<(paraA09->InStr2)<<endl;
	}
	while ((l = kseq_read(seq)) >= 0)
	{
		string Ref=(seq->seq.s) ;
		string OutStr ;
		string ID= seq->name.s ;
		transform(Ref.begin(), Ref.end(), Ref.begin(),::toupper);
		if (seq->comment.l) 
		{
			ID=ID+"\t"+seq->comment.s ;
		}
		int AA=CDS2pep (CODE ,Ref , OutStr ) ;

		if (!(paraA09->TF))
		{
			if (AA==1)
			{
			}
			else if (AA==0)
			{
				cerr<<ID<<"\t"<<"Length Not 3X"<<endl;
			}
			else if (AA==-1)
			{
				cerr<<ID<<"\t"<<"End Not 'U' empty"<<endl;
			}			
			else
			{
				cerr<<ID<<"\t"<<"End Not 'U' other"<<endl;
			}
		}
		Display(OutStr , ID , OUT , linecut );
	}
	OUT.close();
	kseq_destroy(seq);
	gzclose(fp);
	delete paraA09 ;
	return 0;
}
#endif /// CDS2Pep_H_ 
///////// swimming in the sky and flying in the sea ////////////

