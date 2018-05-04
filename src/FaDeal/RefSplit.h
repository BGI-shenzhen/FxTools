#ifndef RefSplit_H_
#define RefSplit_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include "../ALL/comm.h"
#include <cstdlib>
#include <stdio.h>
#include "../ALL/kseq.h"
#include "../ALL/DataClass.h"

using namespace std;

int  print_RefSplit()
{
	cout <<""
		"\n"
		"\tUsage: split  -i <in.fa>  -o   <outDir>\n"
		"\n"
		"\t\t-i    <str>   Input fa File to split\n"
		"\n"
		"\t\t-o    <str>   Output Dir path [PWD]\n"
		"\t\t-s    <int>   Fixed Num Seq for each subFile[1]\n"
		"\t\t-f    <int>   Fixed Num of subFile\n"
		"\n"
		"\t\t-g            No Not gzip the result\n"
		"\n"
		"\t\t-h            show this help\n" 
		"\n";
	return 1;
}

int parse_RScmd(int argc, char **argv , In3str1v * paraA09)
{
	if (argc <=2 ) {print_RefSplit();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check."<< endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InFile"  ||   flag  == "i" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraA09->InStr1=argv[i];
		}
		else if (flag  ==  "OutDir"  ||   flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			paraA09->InStr2=argv[i];
		}
		else if (flag  ==  "SeqNum"  ||   flag  == "s")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			paraA09->InInt=atoi(argv[i]);
			paraA09->TF2=true ;
		}
		else if (flag  ==  "FileNum"   ||   flag  == "f")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			paraA09->InInt=atoi(argv[i]);
			paraA09->TF2=false ;
		}
		else if (flag  == "help"  ||   flag  == "h")
		{
			print_RefSplit();return 0;
		}
		else if (flag  == "NoGzip" ||  flag  == "g")
		{
			paraA09->TF=false ;
		}

		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if  ((paraA09->InStr1).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	return 1 ;

}



int FA_CutSeq_main(   In3str1v * para_A26 )
{

	igzstream  IN ((para_A26->InStr1).c_str(),ifstream::in);
	if(!IN.good())
	{
		cerr << "open InputFile error: "<<(para_A26->InStr1)<<endl;
		 return 0;
	}

	string ext =(para_A26->InStr1).substr((para_A26->InStr1).rfind('/') ==string::npos ? (para_A26->InStr1).length() : (para_A26->InStr1).rfind('/') + 1);
	if (ext=="")
	{
		ext=(para_A26->InStr1);
	}
	ext=replace_all(ext,".fasta.gz","");
	ext=replace_all(ext,".fasta","");
	ext=replace_all(ext,".fa.gz","");
	ext=replace_all(ext,".fa","");
	if (ext=="")
	{           
		ext="Out" ;
	}
	(para_A26->InStr2)=(para_A26->InStr2)+ext+"_cut/";
	mkdir((para_A26->InStr2).c_str() , 0755 ) ;

	string seqAA ;
	getline(IN, seqAA, '>');
	string chr_name ;
	if ((para_A26->InInt)<1)
	{
		cout<<"cut off must be positive"<<endl;
		return 0;
	}
	else if ( (para_A26->InInt) == 1 )
	{
		while(!IN.eof())
		{
			string chr_line;
			getline(IN, chr_line , '\n') ;
			getline(IN, seqAA , '>') ;
			if  (chr_line.length()==0) {continue ;}
			istringstream isone (chr_line,istringstream::in);
			isone>>chr_name ;
			string outfile=(para_A26->InStr2)+chr_name+".fa.gz";
			ogzstream  OUT (outfile.c_str());
			if (OUT.fail())
			{
				cerr << "open OUT File error: "<<outfile<<endl;
				return  0;
			}
			OUT<<">"<<chr_line<<"\n"<<seqAA;
			OUT.close();
			OUT.clear();
		}
	}
	else
	{
		int A=1;
		while(!IN.eof())
		{
			string chr_line ;
			getline(IN, chr_line , '\n') ;
			getline(IN, seqAA , '>') ;
			string outfile=(para_A26->InStr2)+ext+"_"+Int2Str(A)+".fa.gz";
			ogzstream  OUT (outfile.c_str());
			if (OUT.fail())
			{
				cerr << "open OUT File error: "<<outfile<<endl;
				 return  0;
			}            
			OUT<<">"<<chr_line<<"\n"<<seqAA ;

			for ( int i=1 ; i< (para_A26->InInt) ; i++ )
			{
				if (!IN.eof())
				{
					getline(IN, chr_line , '\n');
					getline(IN, seqAA , '>') ;
					OUT<<">"<<chr_line<<"\n"<<seqAA;
				}
				else
				{
					break ;
				}
			}

			OUT.close();
			OUT.clear();
			A++;
		}

	}

	return 0;
}


int FA_CutSeqNOGzip_main(   In3str1v * para_A26 )
{

	igzstream  IN ((para_A26->InStr1).c_str(),ifstream::in);
	if(!IN.good())
	{
		cerr << "open InputFile error: "<<(para_A26->InStr1)<<endl;
		return 0;
	}

	string ext =(para_A26->InStr1).substr((para_A26->InStr1).rfind('/') ==string::npos ? (para_A26->InStr1).length() : (para_A26->InStr1).rfind('/') + 1);
	if (ext=="")
	{
		ext=(para_A26->InStr1);
	}
	ext=replace_all(ext,".fasta.gz","");
	ext=replace_all(ext,".fasta","");
	ext=replace_all(ext,".fa.gz","");
	ext=replace_all(ext,".fa","");
	if (ext=="")
	{           
		ext="Out" ;
	}
	(para_A26->InStr2)=(para_A26->InStr2)+ext+"_cut/";
	mkdir((para_A26->InStr2).c_str() , 0755 ) ;

	string seqAA ;
	getline(IN, seqAA, '>');
	string chr_name ;
	if ((para_A26->InInt)<1)
	{
		cout<<"cut off must be positive"<<endl;
		 return 0;
	}
	else if ( (para_A26->InInt) == 1 )
	{
		while(!IN.eof())
		{
			string chr_line;
			getline(IN, chr_line , '\n') ;
			getline(IN, seqAA , '>') ;
			if  (chr_line.length()==0) {continue ;}
			istringstream isone (chr_line,istringstream::in);
			isone>>chr_name ;
			string outfile=(para_A26->InStr2)+chr_name+".fa";
			ofstream  OUT (outfile.c_str());
			if (OUT.fail())
			{
				cerr << "open OUT File error: "<<outfile<<endl;
				 return  0;
			}
			OUT<<">"<<chr_line<<"\n"<<seqAA;
			OUT.close();
			OUT.clear();
		}
	}
	else
	{
		int A=1;
		while(!IN.eof())
		{
			string chr_line ;
			getline(IN, chr_line , '\n') ;
			getline(IN, seqAA , '>') ;
			string outfile=(para_A26->InStr2)+ext+"_"+Int2Str(A)+".fa";
			ofstream  OUT (outfile.c_str());
			if (OUT.fail())
			{
				cerr << "open OUT File error: "<<outfile<<endl;
				 return  0;
			}            
			OUT<<">"<<chr_line<<"\n"<<seqAA ;

			for ( int i=1 ; i< (para_A26->InInt) ; i++ )
			{
				if (!IN.eof())
				{
					getline(IN, chr_line , '\n');
					getline(IN, seqAA , '>') ;
					OUT<<">"<<chr_line<<"\n"<<seqAA;
				}
				else
				{
					break ;
				}
			}

			OUT.close();
			OUT.clear();
			A++;
		}

	}

	return 0;
}


int FA_CutFileNOGzip_main( In3str1v * para_A26 )
{

	string ext =(para_A26->InStr1).substr((para_A26->InStr1).rfind('/') ==string::npos ? (para_A26->InStr1).length() : (para_A26->InStr1).rfind('/') + 1);
	if (ext=="")
	{
		ext=(para_A26->InStr1);
	}
	ext=replace_all(ext,".fa.gz","");
	ext=replace_all(ext,".fa","");
	if (ext=="")
	{           
		ext="Out" ;
	}
	(para_A26->InStr2)=(para_A26->InStr2)+ext+"_cut/";
	mkdir((para_A26->InStr2).c_str() , 0755 ) ;

	if ((para_A26->InInt)<2)
	{
		cout<<"cut off must be biger 1"<<endl;
		return 0;
	}
	else
	{
		gzFile fp;
		kseq_t *seq;
		int l;
		fp = gzopen((para_A26->InStr1).c_str(), "r");
		seq = kseq_init(fp);
		llong SumLeng=0;
		map <string,llong>  ChrLeng ;

		while ((l = kseq_read(seq)) >= 0)
		{
			string  chr=seq->name.s ;
			llong  chr_length=(seq->seq.l);
			SumLeng+=chr_length;
			ChrLeng[chr]=chr_length ;
		}

		kseq_destroy(seq);
		gzclose(fp);

		llong Sub_len=llong(SumLeng/(para_A26->InInt));
		cout<<"SumLeng: "<<SumLeng<<"\tSubFileLeng: "<<Sub_len<<endl; 
		igzstream  IN ((para_A26->InStr1).c_str(),ifstream::in);
		if(!IN.good())
		{
			cerr << "open InputFile error: "<<(para_A26->InStr1)<<endl;
			 return 0;
		}

		string seqAA ;
		string chr_name;
		getline(IN, seqAA, '>');
		int A=1;

		while(!IN.eof())
		{
			string chr_line ;
			string outfile=(para_A26->InStr2)+ext+"_"+Int2Str(A)+".fa";
			ofstream  OUT (outfile.c_str());

			if (OUT.fail())
			{
				cerr << "open OUT File error: "<<outfile<<endl;
				return  0;
			}
			llong Cur_len=0; 

			while(Cur_len<Sub_len)
			{
				if (!IN.eof())
				{
					getline(IN, chr_line , '\n');
					getline(IN, seqAA , '>');
					OUT<<">"<<chr_line<<"\n"<<seqAA;
					istringstream isone (chr_line,istringstream::in);
					isone>>chr_name ;
					Cur_len+=ChrLeng[chr_name];
				}
				else
				{
					break ;
				}
			}
			OUT.close();
			OUT.clear();
			A++;
		}

	}

	return 0;
}


int FA_CutFile_main ( In3str1v * para_A26 )
{

	string ext =(para_A26->InStr1).substr((para_A26->InStr1).rfind('/') ==string::npos ? (para_A26->InStr1).length() : (para_A26->InStr1).rfind('/') + 1);
	if (ext=="")
	{
		ext=(para_A26->InStr1);
	}
	ext=replace_all(ext,".fa.gz","");
	ext=replace_all(ext,".fa","");
	if (ext=="")
	{           
		ext="Out" ;
	}
	(para_A26->InStr2)=(para_A26->InStr2)+ext+"_cut/";
	mkdir((para_A26->InStr2).c_str() , 0755 ) ;

	if ((para_A26->InInt)<2)
	{
		cout<<"cut off must be biger 1"<<endl;
		return 0;
	}
	else
	{
		gzFile fp;
		kseq_t *seq;
		int l;
		fp = gzopen((para_A26->InStr1).c_str(), "r");
		seq = kseq_init(fp);
		llong SumLeng=0;
		map <string,llong>  ChrLeng ;

		while ((l = kseq_read(seq)) >= 0)
		{
			string  chr=seq->name.s ;
			llong  chr_length=(seq->seq.l);
			SumLeng+=chr_length;
			ChrLeng[chr]=chr_length ;
		}

		kseq_destroy(seq);
		gzclose(fp);

		llong Sub_len=llong(SumLeng/(para_A26->InInt));
		cout<<"SumLeng: "<<SumLeng<<"\tSubFileLeng: "<<Sub_len<<endl; 
		igzstream  IN ((para_A26->InStr1).c_str(),ifstream::in);
		if(!IN.good())
		{
			cerr << "open InputFile error: "<<(para_A26->InStr1)<<endl;
			 return 0;
		}

		string seqAA ;
		getline(IN, seqAA, '>');
		int A=1;

		while(!IN.eof())
		{
			string chr_line ;
			string chr_name;
			string outfile=(para_A26->InStr2)+ext+"_"+Int2Str(A)+".fa.gz";
			ogzstream  OUT (outfile.c_str());

			if (OUT.fail())
			{
				cerr << "open OUT File error: "<<outfile<<endl;
				return  0;
			}
			llong Cur_len=0; 

			while(Cur_len<Sub_len)
			{
				if (!IN.eof())
				{
					getline(IN, chr_line , '\n');
					getline(IN, seqAA , '>');
					OUT<<">"<<chr_line<<"\n"<<seqAA;
					istringstream isone (chr_line,istringstream::in);
					isone>>chr_name ;
					Cur_len+=ChrLeng[chr_name];
				}
				else
				{
					break ;
				}
			}
			OUT.close();
			OUT.clear();
			A++;
		}

	}

	return 0;
}


int FA_Split_main(int argc, char *argv[])
{
	In3str1v * paraA09 = new In3str1v;
	paraA09->InInt=1;
	if( parse_RScmd(argc, argv, paraA09 )==0)
	{
		delete  paraA09 ;
		return 1;
	}
	if ( (paraA09->InStr2).empty())
	{
		(paraA09->InStr2)="./" ;
	}

	if  ((paraA09->TF2)  &&  ( paraA09->TF))
	{
		 FA_CutSeq_main( paraA09  );
	}
	else if ( (paraA09->TF2)  &&  (!paraA09->TF))
	{
		FA_CutSeqNOGzip_main ( paraA09  );
	}
	else if  ((!paraA09->TF2)  &&  ( paraA09->TF))
	{
		FA_CutFile_main  ( paraA09  );
	}
	else
	{
		FA_CutFileNOGzip_main ( paraA09  );
	}

	delete paraA09 ;
	return 0;
}
#endif /// RefSplit_H_ 
///////// swimming in the sky and flying in the sea ////////////



///////// swimming in the sky and flying in the sea ////////////

///////// swimming in the sky and flying in the sea ////////////

