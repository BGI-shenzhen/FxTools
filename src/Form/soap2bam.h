#ifndef soap2bam_H_
#define soap2bam_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <cstdlib>
#include <gzstream.h>

#include <zlib.h>
#include <stdio.h>


#include "../ALL/comm.h"
#include "../ALL/DataClass.h"

#include "./msort/msort.h"
#include "./msort/sort_funs.c"
#include "./msort/stdhashc.cc"
#include "./msort/msort_h.h"

#include <htslib/sam.h>
#include <htslib/kstring.h>

#define UINT32_MAX  (4294967295U)


//KSEQ_INIT(gzFile, gzread)

using namespace std;
typedef long long  llong ;
int msort_main(int argc, char **argv) ;


int  print_Ausage_Formt01()
{
	cout <<""
		"\n"
		"\tUsage: soap2bam  -i <in.soap> -s  <out.sam> \n"
		"\tUsage: soap2bam  -i <in.soap> -b  <out.bam>  -d Ref.fa\n"
		"\n"
		"\t\t-i  <str>   Input SortSoap file \n"
		"\t\t-o  <str>   Output Bam file\n"
		"\n"
		"\t\t-d  <str>   Input Ref.fa or Ref.dict for bam Head\n"
		"\t\t-s  <str>   Output Sam file\n"
		"\t\t-p          if soap is PairOut,for flag\n"
		"\t\t-Q  <int>   shift the seqQ [+31/-31] [0]\n"
		"\t\t-g          if soap is no original\n"
		"\t\t            same ReadID No Neighbors\n"
		"\n"
		"\t\t-h          show this help [hewm2008 v1.02]\n"
		"\n";
	return 1;
}


int parse_Acmd_Formt01(int argc, char **argv , Para_Formt01 * para_Formt01)
{
	if (argc <=2 ) {print_Ausage_Formt01();return 0;}

	for(int i = 1; i < argc; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InSoap"  ||  flag  == "i"	)
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Formt01->input_file=argv[i];
		}           
		else if (flag  ==  "OutSam" ||  flag  == "s" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Formt01->OutSamFile=argv[i];
		}
		else if (flag  ==  "OutBam" ||  flag  == "b")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Formt01->OutBamFile=argv[i];
		}
		else if (flag  ==  "Dict" ||  flag  == "d")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Formt01->Dict=argv[i];
		}
		else if (flag  ==  "ShiftQ " ||  flag  == "Q")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Formt01->shiftQ=atoi(argv[i]);
		}
		else if (flag  ==  "Pair"  ||  flag  == "p")
		{
			(para_Formt01->PE)=1;
		}
		else if (flag  ==  "NoOri"  ||  flag  == "g")
		{
			(para_Formt01->sort)=1;
		}       
		else if (flag  == "help" ||  flag  == "h")
		{
			print_Ausage_Formt01();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}


	if  ((para_Formt01->input_file).empty()  )
	{
		cerr<< "lack argument for the must [-InSoap] "<<endl ;
		return 0;
	}
	if ( (para_Formt01->OutSamFile).empty()    &&    (para_Formt01->OutBamFile).empty() )
	{
		cerr<< "lack argument for the must [-OutSam]  or [-OutBam] "<<endl ;
		return 0;
	}
	else if (   (! (para_Formt01->OutSamFile).empty() ))
	{
		(para_Formt01->OutSamFile)=add_Asuffix( para_Formt01->OutSamFile );
	}
	else
	{
		if ((para_Formt01->Dict).empty())
		{
			cerr<< "[-OutBam] must together with [-Dict] "<<endl ;
			return 0;
		}
		string path=(para_Formt01->OutBamFile);
		string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);		
		if (ext != "bam")
		{
			(para_Formt01->OutBamFile)=path+".bam" ;
		}
	}
	return 1 ;
}




//*///
bam_hdr_t * sam_hdr_sanitise(bam_hdr_t *h)
	//static bam_hdr_t * sam_hdr_sanitise(bam_hdr_t *h)
{
	if (!h)
		return NULL;

	if (h->l_text == 0)
		return h;

	uint32_t i, lnum = 0;
	char *cp = h->text, last = '\n';
	for (i = 0; i < h->l_text; i++) {
		if (cp[i] == 0)
			break;

		if (last == '\n') {
			lnum++;
			if (cp[i] != '@') {
				hts_log_error("Malformed SAM header at line %u", lnum);
				bam_hdr_destroy(h);
				return NULL;
			}
		}

		last = cp[i];
	}

	if (i < h->l_text) { // Early nul found.  Complain if not just padding.
		uint32_t j = i;
		while (j < h->l_text && cp[j] == '\0') j++;
		if (j < h->l_text)
			hts_log_warning("Unexpected NUL character in header. Possibly truncated");
	}
	if (last != '\n') {
		hts_log_warning("Missing trailing newline on SAM header. Possibly truncated");

		if (h->l_text == UINT32_MAX) {
			hts_log_error("No room for extra newline");
			bam_hdr_destroy(h);
			return NULL;
		}

		if (i >= h->l_text - 1) {
			cp = (char*) (realloc(h->text, (size_t) h->l_text+2));
			if (!cp) {
				bam_hdr_destroy(h);
				return NULL;
			}
			h->text = cp;
		}
		cp[i++] = '\n';

		if (h->l_text < i)
			h->l_text = i;
		cp[h->l_text] = '\0';
	}

	return h;
}
////*/


//*/////
bam_hdr_t *Fa_hdr_read( string RefPath)
{
	gzFile fp;
	kseq_t *seq;
	bam_hdr_t *h = NULL;

	kstring_t str = { 0, 0, NULL };

	int l;
	fp = gzopen(RefPath.c_str(), "r");
	seq = kseq_init(fp);

	while ((l = kseq_read(seq)) >= 0)
	{
		stringstream   sstrm;
		sstrm  <<  (seq->seq.l);

		kputs("@SQ\tSN:", &str);
		kputs(seq->name.s, &str);
		kputs("\tLN:", &str);
		kputl(seq->seq.l, &str);
		kputc('\n', &str);
	}

	if (str.l == 0) kputsn("", 0, &str);
	h = sam_hdr_parse(str.l, str.s);
	h->l_text = str.l; h->text = str.s;

	kseq_destroy(seq);
	gzclose(fp);

	return sam_hdr_sanitise(h);
}

///*///*/

int sam2bam2( kstring_t str , bam_hdr_t *h, bam1_t *b)
{
	int ret;
	ret = sam_parse1(&str, h, b);
	if (ret < 0) {
		cerr<<"same the wrong at this "<<str.s ;
	}
	return 1;
}





int matmis( string line ,string & Base ,int & data )
{    
	vector<string> inf;
	split(line,inf,"->");
	if (inf.size()<2)
	{
		return 0 ;
	}
	else
	{
		Base=inf[0];
		int L=inf[1].length();
		int i=0 ;
		for ( i=0; i<L ; i++ )
		{
			char N= (inf[1][i]) ;
			if ( N  == 'A' ||  N == 'C' ||  N == 'T'  ||  N == 'G' )
			{
				break ;
			}
		}
		string B=inf[1].substr(0,i);
		data=atoi(B.c_str());
		return 1 ;
	}
	//    if ($ARGV[0] =~ /^([ACGT])->(\d+)/i);
}

int mating ( SamLine * A ,  SamLine * B )
{
	int Insert=0;
	if (  (A->chr) != "*" && ( (A->chr) ==(B->chr))  )  //# then calculate $isize
	{
		llong  x1=( (A->Flag) & 0x10 ) ? (A->position)+((A->seq).length()) : (A->position);
		llong  x2=( (B->Flag) & 0x10 ) ? (B->position)+((B->seq).length()) : (B->position);
		Insert=x2-x1;
	}
	//# update mate coordinate
	if  ( (B->chr) != "*" )
	{
		A->coor=B->position;
		A->isize=Insert;
		A->XorD= ((B->chr) == (A->chr)) ? "=" : (B->chr) ; 
		if ((B->Flag) & 0x10 )
		{
			A->Flag |= 0x20 ;
		}
	}
	else
	{
		A->Flag |= 0x8;
	}

	if  ( (A->chr) != "*" )
	{
		B->coor=A->position;
		B->isize=0-Insert;
		B->XorD= ((A->chr) == (B->chr)) ? "=" : (A->chr) ; 
		if ((A->Flag) & 0x10 )
		{
			B->Flag |= 0x20 ;
		}
	}
	else
	{
		B->Flag |= 0x8;
	}    
	return 1;
}


int soap2samQ(string line , Para_Formt01 * para_Formt01 , SamLine *sam )
{
	vector<string> inf;
	split(line,inf," \t");
	int Length=inf.size();
	if ( Length < 9 || (inf[0].empty()) )
	{
		return 0 ;
	}
	sam->rm();
	int A=inf[3][0];
	if ( A>57 || A<48 ) // fix SOAP-2.1.x bugs
	{
		vector<string> tmp (Length-1);
		tmp[0]=inf[0]; tmp[1]=inf[1];   tmp[2]=inf[2]; 
		for(int ii=4 ; ii<Length ; ii++ )
		{
			tmp[ii-1]=tmp[ii];
		}
		inf.clear();
		inf=tmp ;
	}
	// readID //
	int RID_length=inf[0].length() ;
	//    sam->RID=getID(inf[0]);
	//*
	sam->RID=(inf[0]); 
	if ( RID_length >2 )
	{
		if ( inf[0][RID_length-2]== '/'  &&  ( (inf[0][RID_length-1]== '1')  || (inf[0][RID_length-1]== '2') ))
		{
			sam->RID=inf[0].substr(0, RID_length-2);                   
		}
	}
	///*///

	// initial flag (will be updated later)
	sam->Flag =0 ; 

	(sam->Flag) |= 1 | 1<<(inf[4] == "a" ? 6 : 7);
	if  (para_Formt01->PE)
	{
		(sam->Flag)  |= 2;
	}
	if (inf[6] == "-" )
	{           
		(sam->Flag) |= 0x10 ;
	}

	// # read & quality
	sam->seq=inf[1];
	int RLength=inf[1].length();
	int QLength=inf[2].length();
	sam->Qseq = (QLength> RLength ) ? inf[2].substr(0,RLength): inf[2] ;
	//*/////
	QLength=(sam->Qseq).length();
	for(string::size_type ix=0; ix<QLength; ix++)
	{
		(sam->Qseq)[ix]-=(para_Formt01->shiftQ) ; ////// chang the "#"  Quli to "B" Quli
	}
	///*////
	// cigar 
	sam->cigar= Int2Str(RLength)+"M";

	// # coor
	sam->chr = inf[7]; sam->position = atoi(inf[8].c_str());

	// mapQ
	sam->mapQ  = (inf[3] == "1" ) ? 30 : 0;

	// # mate coordinate
	//$s->[6] = '*'; $s->[7] = $s->[8] = 0;

	// # aux 
	sam->NM_i="NM:i:"+inf[9];
	sam->MD="";
	A=atoi(inf[9].c_str());
	if (A)
	{
		string Base ;
		int data ;
		map <int ,string  > MAP ;
		for (int ii=10 ; ii<Length ; ii++ )
		{
			if (matmis( inf[ii], Base , data))
			{
				MAP.insert(map < int,string > :: value_type (data,Base));
			}
		}
		int a=0;
		map <int , string> :: iterator it=MAP.begin();
		for(it=MAP.begin() ; it!=MAP.end(); it++)
		{
			int c =(it->first) - a ;
			(sam->MD) += Int2Str(c) + (it->second);
			a += (c  + 1);
		}
		(sam->MD) += Int2Str(RLength-a);
	}
	else
	{
		(sam->MD)=Int2Str(RLength);
	}
	(sam->MD)="MD:Z:"+(sam->MD);
	sam->IF=true ;
	return  1;
}
///////////////////





int soap2sam(string line , Para_Formt01 * para_Formt01 , SamLine *sam )
{
	vector<string> inf;
	split(line,inf," \t");
	int Length=inf.size();
	if ( Length < 9 || (inf[0].empty()) )
	{
		return 0 ;
	}
	sam->rm();
	int A=inf[3][0];
	if ( A>57 || A<48 ) // fix SOAP-2.1.x bugs
	{
		vector<string> tmp (Length-1);
		tmp[0]=inf[0]; tmp[1]=inf[1];   tmp[2]=inf[2]; 
		for(int ii=4 ; ii<Length ; ii++ )
		{
			tmp[ii-1]=tmp[ii];
		}
		inf.clear();
		inf=tmp ;
	}
	// readID //
	int RID_length=inf[0].length() ;
	//    sam->RID=getID(inf[0]);
	//*
	sam->RID=(inf[0]); 
	if ( RID_length >2 )
	{
		if ( inf[0][RID_length-2]== '/'  &&  ( (inf[0][RID_length-1]== '1')  || (inf[0][RID_length-1]== '2') ))
		{
			sam->RID=inf[0].substr(0, RID_length-2);                   
		}
	}
	///*///

	// initial flag (will be updated later)
	sam->Flag =0 ; 

	(sam->Flag) |= 1 | 1<<(inf[4] == "a" ? 6 : 7);
	if  (para_Formt01->PE)
	{
		(sam->Flag)  |= 2;
	}
	if (inf[6] == "-" )
	{           
		(sam->Flag) |= 0x10 ;
	}

	// # read & quality
	sam->seq=inf[1];
	int RLength=inf[1].length();
	int QLength=inf[2].length();
	sam->Qseq = (QLength> RLength ) ? inf[2].substr(0,RLength): inf[2] ;


	// cigar 
	sam->cigar= Int2Str(RLength)+"M";

	// # coor
	sam->chr = inf[7]; sam->position = atoi(inf[8].c_str());

	// mapQ
	sam->mapQ  = (inf[3] == "1" ) ? 30 : 0;

	// # mate coordinate
	//$s->[6] = '*'; $s->[7] = $s->[8] = 0;

	// # aux 
	sam->NM_i="NM:i:"+inf[9];
	sam->MD="";
	A=atoi(inf[9].c_str());
	if (A)
	{
		string Base ;
		int data ;
		map <int ,string  > MAP ;
		for (int ii=10 ; ii<Length ; ii++ )
		{
			if (matmis( inf[ii], Base , data))
			{
				MAP.insert(map < int,string > :: value_type (data,Base));
			}
		}
		int a=0;
		map <int , string> :: iterator it=MAP.begin();
		for(it=MAP.begin() ; it!=MAP.end(); it++)
		{
			int c =(it->first) - a ;
			(sam->MD) += Int2Str(c) + (it->second);
			a += (c  + 1);
		}
		(sam->MD) += Int2Str(RLength-a);
	}
	else
	{
		(sam->MD)=Int2Str(RLength);
	}
	(sam->MD)="MD:Z:"+(sam->MD);
	sam->IF=true ;
	return  1;
}
///////////////////








int soap2sam_mainfun ( Para_Formt01 * para_Formt01 )
{
	ogzstream OUT ((para_Formt01->OutSamFile).c_str());
	if(!OUT.good())
	{
		cerr << "open InputFile error: "<<(para_Formt01->OutSamFile)<<endl;
		return 0;
	}

	if (!(para_Formt01->Dict).empty())
	{
		string path=(para_Formt01->Dict);
		string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);

		if (ext=="dict")
		{
			Write_Sam_head (para_Formt01->Dict,OUT) ;
		}
		else
		{
			int AA=path.rfind('/') ==string::npos ? path.length() : path.rfind('/')+1 ;
			string File_Pre =path.substr(0,AA );
			ext =path.substr(AA);
			ext=replace_all(ext,".fasta.gz","");
			ext=replace_all(ext,".fasta","");
			ext=replace_all(ext,".fa.gz","");
			ext=replace_all(ext,".fa","");
			string newDictPath=File_Pre+ext+".dict";

			if ( access(newDictPath.c_str(), 0) == 0 )
			{
				Write_Sam_head (newDictPath , OUT) ;
			}
			else
			{
				string newDictPathV2=(para_Formt01->Dict)+".dict";
				if ( access(newDictPathV2.c_str(), 0) == 0 )
				{
					Write_Sam_head (newDictPathV2 , OUT) ;
				}
				else
				{
					string newDictPathV3=(para_Formt01->Dict)+".chrlist";
					if ( access(newDictPathV3.c_str(), 0) == 0 )
					{
						ifstream IN (newDictPathV3.c_str(),ifstream::in);
						string  line ;
						getline(IN,line);
						while(!IN.eof())
						{
							getline(IN,line);
							if (line.length()<=0)  { continue ; }
							vector<string> inf;
							split(line,inf," \t");
							OUT<<"@SQ\tSN:"<<inf[0]<<"\tLN:"<<inf[1]<<endl;
						}
						IN.close();
					}
					else
					{
						gzFile fp;
						kseq_t *seq;
						int l;
						fp = gzopen((para_Formt01->Dict).c_str(), "r");
						seq = kseq_init(fp);
						while ((l = kseq_read(seq)) >= 0)
						{
							OUT<<"@SQ\tSN:"<<(seq->name.s)<<"\tLN:"<<(seq->seq.l)<<endl;
						}
						kseq_destroy(seq);
						gzclose(fp);
					}
				}
			}
		}
	}




	if ( (para_Formt01->shiftQ)==0)
	{



		if (para_Formt01->sort)
		{
			char * A = const_cast<char*>((para_Formt01->input_file).c_str());
			//char * B = const_cast<char*>(tmpo.c_str());
			char * ssTmp[3]={(char *)"msort", (char *)"-k1", A };
			file_t *db;
			string  outPut ;
			db=ReadParaAnd ( 3 , ssTmp , outPut) ;
			sort_file_lines(db);   

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;

			for(int i=0;i<db->line_size;i++)
			{
				string line=db->lines[i].str ;
				if ( soap2sam (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID))) 
					{
						mating( sam1 ,  sam2  );
						sam2->Print(OUT); sam2->rm();              
						sam1->Print(OUT); sam1->rm();
					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;
							sam2->Print(OUT);
							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				sam2->Print(OUT); sam2->rm();
			}
			free_file_lines(db); 
			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
		}
		else
		{
			igzstream IN ((para_Formt01->input_file).c_str(),ifstream::in);

			if(!IN.good())
			{
				cerr << "open InputFile error: "<<(para_Formt01->input_file)<<endl;
				return 1;
			}

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;

			while(!IN.eof())
			{
				string  line ;
				getline(IN,line);
				if ( soap2sam (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID) )) 
					{
						mating( sam1 ,  sam2  );
						sam2->Print(OUT); sam2->rm();              
						sam1->Print(OUT); sam1->rm();
					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;
							sam2->Print(OUT);
							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				sam2->Print(OUT); sam2->rm();
			}

			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
			IN.close();    
		}




	}

	else
	{



		if (para_Formt01->sort)
		{
			char * A = const_cast<char*>((para_Formt01->input_file).c_str());
			char * ssTmp[3]={(char *)"msort", (char *)"-k1", A };
			file_t *db;
			string  outPut ;
			db=ReadParaAnd ( 3 , ssTmp , outPut) ;
			sort_file_lines(db);   

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;

			for(int i=0;i<db->line_size;i++)
			{
				string line=db->lines[i].str ;
				if ( soap2samQ (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID))) 
					{
						mating( sam1 ,  sam2  );
						sam2->Print(OUT); sam2->rm();              
						sam1->Print(OUT); sam1->rm();
					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;
							sam2->Print(OUT);
							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				sam2->Print(OUT); sam2->rm();
			}
			free_file_lines(db); 
			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
		}
		else
		{
			igzstream IN ((para_Formt01->input_file).c_str(),ifstream::in);

			if(!IN.good())
			{
				cerr << "open InputFile error: "<<(para_Formt01->input_file)<<endl;
				return 1;
			}

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;

			while(!IN.eof())
			{
				string  line ;
				getline(IN,line);
				if ( soap2samQ (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID) )) 
					{
						mating( sam1 ,  sam2  );
						sam2->Print(OUT); sam2->rm();              
						sam1->Print(OUT); sam1->rm();
					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;
							sam2->Print(OUT);
							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				sam2->Print(OUT); sam2->rm();
			}

			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
			IN.close();    
		}



	}

	OUT.close();
	return 1;
}




int soap2bam_mainfun ( Para_Formt01 * para_Formt01 )
{

	htsFile *OUTBam = hts_open((para_Formt01->OutBamFile).c_str(), "wb");
	bam_hdr_t *header=NULL;

	string path=(para_Formt01->Dict);

	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);

	if (ext=="dict")
	{
		ifstream IN (path.c_str(),ifstream::in);
		kstring_t str = { 0, 0, NULL };
		while(!IN.eof())
		{
			string  line ;
			getline(IN,line);
			if (line.length()<=0)  { continue ; }
			kputs(line.c_str(), &str);
			kputc('\n', &str);
		}
		IN.close();
		//		if (str.l == 0) kputsn("", 0, &str);
		header= sam_hdr_parse(str.l, str.s);
		header->l_text = str.l; header->text = str.s;
		header=sam_hdr_sanitise(header);
	}
	else
	{
		int AA=path.rfind('/') ==string::npos ? path.length() : path.rfind('/')+1 ;
		string File_Pre =path.substr(0,AA);
		ext =path.substr(AA);
		ext=replace_all(ext,".fasta.gz","");
		ext=replace_all(ext,".fasta","");
		ext=replace_all(ext,".fa.gz","");
		ext=replace_all(ext,".fa","");
		ext=replace_all(ext,".dict","");
		string newDictPath=File_Pre+ext+".dict";
		if ( access(newDictPath.c_str(), 0) == 0 )
		{

			ifstream IN (newDictPath.c_str(),ifstream::in);
			kstring_t str = { 0, 0, NULL };
			while(!IN.eof())
			{
				string  line ;
				getline(IN,line);
				if (line.length()<=0)  { continue ; }
				kputs(line.c_str(), &str);
				kputc('\n', &str);

			}
			IN.close();
			header= sam_hdr_parse(str.l, str.s);
			header->l_text = str.l; header->text = str.s;
			header=sam_hdr_sanitise(header);
		}

		string newDictPathV2=(para_Formt01->Dict)+".dict";
		if ( access(newDictPathV2.c_str(), 0) == 0 )
		{
			ifstream IN (newDictPathV2.c_str(),ifstream::in);
			kstring_t str = { 0, 0, NULL };
			while(!IN.eof())
			{
				string  line ;
				getline(IN,line);
				if (line.length()<=0)  { continue ; }
				kputs(line.c_str(), &str);
				kputc('\n', &str);
			}
			IN.close();
			header= sam_hdr_parse(str.l, str.s);
			header->l_text = str.l; header->text = str.s;
			header=sam_hdr_sanitise(header);
		}
		else
		{
			string newDictPathV3=(para_Formt01->Dict)+".chrlist";
			if ( access(newDictPathV3.c_str(), 0) == 0 )
			{

				ifstream IN (newDictPathV3.c_str(),ifstream::in);
				kstring_t str = { 0, 0, NULL };
				string  line ;
				getline(IN,line);
				while(!IN.eof())
				{
					getline(IN,line);
					if (line.length()<=0)  { continue ; }
					vector<string> inf;
					split(line,inf," \t");
					string Nline="@SQ\tSN:"+inf[0]+"\tLN:"+inf[1]+"\n";
					kputs(Nline.c_str(), &str);
				}
				IN.close();

				header= sam_hdr_parse(str.l, str.s);
				header->l_text = str.l; header->text = str.s;
				header=sam_hdr_sanitise(header);

			}
			else
			{
				header=Fa_hdr_read(para_Formt01->Dict);
			}
		}
	}


	if (sam_hdr_write(OUTBam, header) < 0)
	{
		fprintf(stderr, "Error writing output.\n");
		exit(-1);
	}

	bam1_t *aln = bam_init1();




	if ( (para_Formt01->shiftQ)==0)
	{


		if (para_Formt01->sort)
		{
			char * A = const_cast<char*>((para_Formt01->input_file).c_str());
			//char * B = const_cast<char*>(tmpo.c_str());
			char * ssTmp[3]={(char *)"msort", (char *)"-k1", A };
			file_t *db;
			string  outPut ;
			db=ReadParaAnd ( 3 , ssTmp , outPut) ;
			sort_file_lines(db);   

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;
			string OutStrA ,OutStrB;

			for(int i=0;i<db->line_size;i++)
			{
				string line=db->lines[i].str ;
				if ( soap2sam (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID))) 
					{
						mating( sam1 ,  sam2  );


						kstring_t str2 = { 0, 0, NULL };
						sam2->OUT2str(OutStrB);
						kputs(OutStrB.c_str(), &str2);
						sam2bam2( str2 , header , aln);
						int aa=sam_write1(OUTBam, header, aln);
						sam2->rm();     

						kstring_t str1 = { 0, 0, NULL };
						sam1->OUT2str(OutStrA);
						kputs(OutStrA.c_str(), &str1);
						sam2bam2( str1 , header , aln);
						aa=sam_write1(OUTBam, header, aln);
						sam1->rm();     

					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;
							kstring_t str2 = { 0, 0, NULL };
							sam2->OUT2str(OutStrB);
							kputs(OutStrB.c_str(), &str2);
							sam2bam2( str2 , header , aln);
							int aa=sam_write1(OUTBam, header, aln);
							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				kstring_t str2 = { 0, 0, NULL };
				sam2->OUT2str(OutStrB);
				kputs(OutStrB.c_str(), &str2);
				sam2bam2( str2 , header , aln);
				int aa=sam_write1(OUTBam, header, aln);
				sam2->rm();  
			}
			free_file_lines(db); 
			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
		}

		else
		{
			igzstream IN ((para_Formt01->input_file).c_str(),ifstream::in);

			if(!IN.good())
			{
				cerr << "open InputFile error: "<<(para_Formt01->input_file)<<endl;
				return 1;
			}

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;
			string OutStrA,OutStrB;

			while(!IN.eof())
			{
				string  line ;
				getline(IN,line);
				if ( soap2sam (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID) )) 
					{
						mating( sam1 ,  sam2  );


						kstring_t str2 = { 0, 0, NULL };
						sam2->OUT2str(OutStrB);
						kputs(OutStrB.c_str(), &str2);
						sam2bam2( str2 , header , aln);
						int aa=sam_write1(OUTBam, header, aln);
						sam2->rm();     

						kstring_t str1 = { 0, 0, NULL };
						sam1->OUT2str(OutStrA);
						kputs(OutStrA.c_str(), &str1);
						sam2bam2( str1 , header , aln);
						aa=sam_write1(OUTBam, header, aln);
						sam1->rm();     

					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;

							kstring_t str2 = { 0, 0, NULL };
							sam2->OUT2str(OutStrB);
							kputs(OutStrB.c_str(), &str2);
							sam2bam2( str2 , header , aln);
							int aa=sam_write1(OUTBam, header, aln);
							sam2->rm();  

							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				kstring_t str2 = { 0, 0, NULL };
				sam2->OUT2str(OutStrB);
				kputs(OutStrB.c_str(), &str2);
				sam2bam2( str2 , header , aln);
				int aa=sam_write1(OUTBam, header, aln);
				sam2->rm();     
			}

			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
			IN.close();    
		}


	}
	else
	{




		if (para_Formt01->sort)
		{
			char * A = const_cast<char*>((para_Formt01->input_file).c_str());
			//char * B = const_cast<char*>(tmpo.c_str());
			char * ssTmp[3]={(char *)"msort", (char *)"-k1", A };
			file_t *db;
			string  outPut ;
			db=ReadParaAnd ( 3 , ssTmp , outPut) ;
			sort_file_lines(db);   

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;
			string OutStrA ,OutStrB;

			for(int i=0;i<db->line_size;i++)
			{
				string line=db->lines[i].str ;
				if ( soap2samQ (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID))) 
					{
						mating( sam1 ,  sam2  );


						kstring_t str2 = { 0, 0, NULL };
						sam2->OUT2str(OutStrB);
						kputs(OutStrB.c_str(), &str2);
						sam2bam2( str2 , header , aln);
						int aa=sam_write1(OUTBam, header, aln);
						sam2->rm();     

						kstring_t str1 = { 0, 0, NULL };
						sam1->OUT2str(OutStrA);
						kputs(OutStrA.c_str(), &str1);
						sam2bam2( str1 , header , aln);
						aa=sam_write1(OUTBam, header, aln);
						sam1->rm();     

					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;
							kstring_t str2 = { 0, 0, NULL };
							sam2->OUT2str(OutStrB);
							kputs(OutStrB.c_str(), &str2);
							sam2bam2( str2 , header , aln);
							int aa=sam_write1(OUTBam, header, aln);
							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				kstring_t str2 = { 0, 0, NULL };
				sam2->OUT2str(OutStrB);
				kputs(OutStrB.c_str(), &str2);
				sam2bam2( str2 , header , aln);
				int aa=sam_write1(OUTBam, header, aln);
				sam2->rm();  
			}
			free_file_lines(db); 
			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
		}

		else
		{
			igzstream IN ((para_Formt01->input_file).c_str(),ifstream::in);

			if(!IN.good())
			{
				cerr << "open InputFile error: "<<(para_Formt01->input_file)<<endl;
				return 1;
			}

			SamLine  *sam1 = new SamLine ;
			SamLine  *sam2 = new SamLine ;
			SamLine  *samtmp = new SamLine ;
			string OutStrA,OutStrB;

			while(!IN.eof())
			{
				string  line ;
				getline(IN,line);
				if ( soap2samQ (line ,  para_Formt01 , sam1 ))
				{
					if ( (sam2->IF)  && ( (sam1->RID) == (sam2->RID) )) 
					{
						mating( sam1 ,  sam2  );



						kstring_t str2 = { 0, 0, NULL };
						sam2->OUT2str(OutStrB);
						kputs(OutStrB.c_str(), &str2);
						sam2bam2( str2 , header , aln);
						int aa=sam_write1(OUTBam, header, aln);
						sam2->rm();

						kstring_t str1 = { 0, 0, NULL };
						sam1->OUT2str(OutStrA);
						kputs(OutStrA.c_str(), &str1);
						sam2bam2( str1 , header , aln);
						aa=sam_write1(OUTBam, header, aln);
						sam1->rm();

					}
					else
					{
						if ((sam2->IF))
						{
							sam2->Flag |= 0x8 ;

							kstring_t str2 = { 0, 0, NULL };
							sam2->OUT2str(OutStrB);
							kputs(OutStrB.c_str(), &str2);
							sam2bam2( str2 , header , aln);
							int aa=sam_write1(OUTBam, header, aln);
							sam2->rm();  

							sam2->Flag |= 0x8 ;
						}
						samtmp->copy(sam2);
						sam2->copy(sam1); 
						sam1->copy(samtmp);
					}
				}
			}

			if ((sam2->IF))
			{
				kstring_t str2 = { 0, 0, NULL };
				sam2->OUT2str(OutStrB);
				kputs(OutStrB.c_str(), &str2);
				sam2bam2( str2 , header , aln);
				int aa=sam_write1(OUTBam, header, aln);
				sam2->rm();
			}

			delete sam1 ;
			delete sam2 ;
			delete samtmp ;
			IN.close();
		}







	}

	bam_hdr_destroy(header);
	bam_destroy1(aln);
	sam_close(OUTBam);
	return 1;
}





int Soap2Xam_main(int argc, char *argv[])
{
	Para_Formt01 * para_Formt01 = new Para_Formt01;
	if (parse_Acmd_Formt01(argc, argv , para_Formt01 )==0)
	{
		delete  para_Formt01 ;
		return 1;
	}


	if (!((para_Formt01->OutBamFile).empty()))
	{
		soap2bam_mainfun (  para_Formt01 );
	}
	else
	{
		soap2sam_mainfun (  para_Formt01 );
	}

	delete para_Formt01 ;
	return 0;

}


///////// swimming in the sky and flying in the sea ////////////
//


#endif // soap2bam_H_


////////////////////////swimming in the sea & flying in the sky //////////////////




#ifndef sameFun_H_
#define SameFun_H_


#endif 
///////// swimming in the sky and flying in the sea ////////////


