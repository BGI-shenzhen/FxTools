#ifndef RefGetCDS_H_
#define RefGetCDS_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "../ALL/comm.h"
#include "../ALL/DataClass.h"
#include <cstdlib>
#include <gzstream.h>
#include <zlib.h>
#include <stdio.h>
#include "../ALL/kseq.h"
#include "CDS2Pep.h"

//KSEQ_AINIT(gzFile, gzread)

using namespace std;
typedef long long  llong ;

int  print_AusageFA04()
{
	cout <<""
		"\n"
		"\tUsage: getCdsPep -r <in.fa> -g <in.gff>  -o <out> \n"
		"\n"
		"\t\t-r       <str>   Input Ref fa\n"
		"\t\t-g       <str>   Input Ref gff\n"
		"\t\t-o       <str>   Output cds/pep seq file(prefix)\n"        
		"\n"
		"\t\t-s               OutPut the 4D Site(fourfold degenerate site)\n"
		"\n"
		"\t\t-h               show this help\n" 
		"\n";
	return 1;
}


int parse_AcmdFA04(int argc, char **argv , In3str1v * paraFA04 )
{
	if (argc <=2 ) {print_AusageFA04();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "Ref" ||  flag  == "r")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut"  ||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
		}
		else if (flag  ==  "Gff"||  flag  == "g")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr3=argv[i];
		}
		else if (flag  == "help" ||  flag  == "h")
		{
			print_AusageFA04();return 0;
		}
		else if (flag  == "4DSite"||  flag  == "s")
		{
			paraFA04->TF=false ;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((paraFA04->InStr1).empty() || (paraFA04->InStr2).empty() ||  (paraFA04->InStr3).empty()   )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}


	return 1 ;
}

//int main(int argc, char *argv[])
int FA_GetCDS_main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	if( parse_AcmdFA04(argc, argv, paraFA04 )==0)
	{
		delete  paraFA04 ;
		return 1 ;
	}

	map <string,map <string ,bool> > GeneList;
	map <string,map <llong,llong> >  CDSList ;

	igzstream  MRG ((paraFA04->InStr3).c_str(),ifstream::in);
	if(!MRG.good())
	{
		cerr << "open InputFile error: "<<(paraFA04->InStr3)<<endl;
		delete  paraFA04 ; return 1;
	}

	while(!MRG.eof())
	{
		string  line ;
		getline(MRG,line);
		if (line.length()<=0  || line[0] == '#' )  { continue  ; }
		istringstream isone (line,istringstream::in);
		string chr , flag , CDS , ZhengFu ,geneID ;
		llong Start,End ;
		isone>>chr>>flag>>CDS ;
		if (CDS  != "CDS" )  { continue  ; }
		isone>>Start>>End>>flag>>ZhengFu>>flag>>geneID ;
		vector<string> inf;
		vector<string> Temp;
		split(geneID,inf,",;");
		split(inf[0],Temp,"=");
		string GeneID= Temp[Temp.size()-1] ;
		for ( int jk=1;  jk<inf.size() ; jk++)
		{
			vector<string> Temptmp2;
			split(inf[jk],Temptmp2,"=");
			if (Temptmp2[0] == "Parent")
			{
				GeneID= Temptmp2[Temptmp2.size()-1] ;
			}
		}

		map  <string,map <llong,llong> >  :: iterator it=CDSList.find(GeneID);
		if (it == CDSList.end() )
		{
			bool A=false ;
			if ( ZhengFu == "-" )
			{
				A=true ;
			}
			map <string,map <string ,bool> > :: iterator X =GeneList.find(chr);
			if (X==GeneList.end())
			{
				map <string,bool> First;
				First[GeneID]=A;
				GeneList.insert(map <string,map <string ,bool> > :: value_type(chr,First));
			}
			else
			{
				(X->second).insert(map <string ,bool>  :: value_type(GeneID ,A ));
			}
			map <llong,llong> DD;
			DD[Start]=End ;
			CDSList.insert(map <string,map <llong,llong> > ::value_type(GeneID,DD));
			//CDSList
		}
		else
		{
			(it->second).insert(map <llong,llong>  :: value_type(Start,End)) ;
		}
	}

	MRG.close();


	int linecut= FaCutLine ((paraFA04->InStr1));


	gzFile fp;
	kseq_t *seq;
	int l;

	//    (paraFA04->InStr2)=add_Asuffix(paraFA04->InStr2) ;
	string outpath=(paraFA04->InStr2)+".cds.fa.gz";
	fp = gzopen((paraFA04->InStr1).c_str(), "r");
	seq = kseq_init(fp);    
	ogzstream OUT ((outpath).c_str());

	if(!OUT.good())
	{
		cerr << "open OUT File error: "<<(outpath)<<endl;
		delete  paraFA04 ; return  1;
	}

	map <char ,char >  Complement;
	Complement['A']='T'; Complement['C']='G';  Complement['T']='A'; Complement['G']='C'; Complement['N']='N';
	Complement['a']='t'; Complement['c']='g';  Complement['t']='a'; Complement['g']='c'; Complement['n']='n';
	Complement['M']='K'; Complement['R']='Y';  Complement['W']='W'; Complement['S']='S'; Complement['Y']='R';
	Complement['K']='M';
	if(paraFA04->TF)
	{
		while ((l = kseq_read(seq)) >= 0)
		{
			string chr=seq->name.s ;
			map <string,map <string,bool> > :: iterator it=GeneList.find(chr);
			if (it != GeneList.end())
			{
				string Ref= seq->seq.s ;
				map <string,bool> ::iterator  Y ;
				for(  Y=(it->second).begin() ; Y!= (it->second).end(); Y++ )
				{
					bool TF=Y->second;
					map  <string,map <llong,llong> > :: iterator mCDS=CDSList.find(Y->first);
					string BB ;
					map <llong,llong> :: iterator PUPU ;
					for(  PUPU=(mCDS->second).begin() ; PUPU!= (mCDS->second).end(); PUPU++ )
					{
						llong START=(PUPU->first)-1 ;
						int Leng=(PUPU->second)-START;
						BB=BB+Ref.substr(START,Leng) ;
					}
					if (TF)
					{
						reverse(BB.begin(), BB.end());
						int leng=BB.size();
						for (int i=0 ; i<leng ; i++)
						{
							BB[i]=Complement[BB[i]];
						}
					}
					transform(BB.begin(), BB.end(), BB.begin(),::toupper);
					string ID=(Y->first) ;
					Display(BB , ID ,  OUT ,linecut  );
				}
			}
		}
	}
	else
	{



		string out4DTV=(paraFA04->InStr2)+".4Dsite.gz";
		ogzstream OUT4DTV ((out4DTV).c_str());

		map <string,bool> CODE4DTV ;
		map <string,bool> ::iterator map_CODE4DTV_it ;
		CODE4DTV["GC"]=true ;
		CODE4DTV["GG"]=true ;
		CODE4DTV["CT"]=true ;
		CODE4DTV["CC"]=true ;
		CODE4DTV["CG"]=true ;
		CODE4DTV["TC"]=true ;
		CODE4DTV["AC"]=true ;
		CODE4DTV["GT"]=true ;

		while ((l = kseq_read(seq)) >= 0)
		{
			string chr=seq->name.s ;
			map <string,map <string,bool> > :: iterator it=GeneList.find(chr);
			if (it != GeneList.end())
			{
				string Ref= seq->seq.s ;
				map <string,bool> ::iterator  Y ;
				for(  Y=(it->second).begin() ; Y!= (it->second).end(); Y++ )
				{
					bool TF=Y->second ;
					map  <string,map <llong,llong> > :: iterator mCDS=CDSList.find(Y->first);
					string BB ;
					map <llong,llong> :: iterator PUPU ;
					for(  PUPU=(mCDS->second).begin() ; PUPU!= (mCDS->second).end(); PUPU++ )
					{
						llong START=(PUPU->first)-1 ;
						int Leng=(PUPU->second)-START;
						BB=BB+Ref.substr(START,Leng) ;
					}
					
					if (TF)
					{
						reverse(BB.begin(), BB.end());
						int leng=BB.size();
						for (int i=0 ; i<leng ; i++)
						{
							BB[i]=Complement[BB[i]];
						}
					}

					transform(BB.begin(), BB.end(), BB.begin(),::toupper);
					string ID=(Y->first);
					Display(BB ,ID, OUT ,linecut);
		
					


					map <int,string> Site4DTV ;
					int length=BB.length();
					string temp_vv ="";
					for(int ii=0 ; ii<length; ii+=3)
					{
						string base=BB.substr(ii,3);
						string base3CDS=base.substr(0,2); ;
						map_CODE4DTV_it=CODE4DTV.find(base3CDS);
						if (map_CODE4DTV_it==CODE4DTV.end())
						{
							continue ;
						}

						int jjk=ii+2;
						Site4DTV[jjk]=base;
					}

					if  (Site4DTV.empty())
					{
						continue ;
					}

					if (TF)
					{
						map<llong,llong>::reverse_iterator rPUPU ;
					for(  rPUPU=(mCDS->second).rbegin() ; rPUPU!= (mCDS->second).rend(); rPUPU++ )
					{
						llong START=(rPUPU->first);
						llong END=(rPUPU->second);
						int iik=0;
						map <int,string>  :: iterator  It_4D ;
						for (   ; END>=START ;   END--)
						{
							It_4D=Site4DTV.find(iik);
							if (It_4D!=Site4DTV.end())
							{
								OUT4DTV<<chr<<"\t"<<END<<"\t"<<ID<<":-\t"<<iik<<"\t"<<It_4D->second<<"\n";
							}
							iik++;
						}
					}



					}
					else
					{
					for(  PUPU=(mCDS->second).begin() ; PUPU!= (mCDS->second).end(); PUPU++ )
					{
						llong START=(PUPU->first);
						llong END=(PUPU->second);
						int iik=0;
						map <int,string>  :: iterator  It_4D ;
						for (   ; START<=END ;   START++)
						{
							It_4D=Site4DTV.find(iik);
							if (It_4D!=Site4DTV.end())
							{
								OUT4DTV<<chr<<"\t"<<START<<"\t"<<ID<<":+\t"<<iik<<"\t"<<It_4D->second<<"\n";
							}
							iik++;
						}
					}

					}
				}
			}
		}




		OUT4DTV.close();  

	}


	OUT.close();
	kseq_destroy(seq);
	gzclose(fp);
	string PP= (paraFA04->InStr2)+".pep.fa.gz";
	char * A = const_cast<char*>((PP).c_str());
	char * B = const_cast<char*>((outpath).c_str());
	char * ssTmp[5]={(char *)"CDS2Pep", (char *)"-InCDS", B ,(char *)"-OutPut",A};
	FA_CDS2Pep_main( 5 , ssTmp ) ;
	delete paraFA04 ;
	return 0;
}
#endif // RefGetCDS_H_ //
///////// swimming in the sky and flying in the sea ////////////

