#ifndef ChangPosi_h_
#define ChangPosi_h_

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

//KSEQ_AINIT(gzFile, gzread)
using namespace std;

typedef long long  llong ;

int  print_AusageFABP()
{
	cout <<""
		"\n"
		"\tUsage: ChangPosi -i <in.snp> -l <Ref.fa.merlist> -o <Out.snp>\n"
		"\n"
		"\t\t-i     <str>   Input the Position File to back Scaffolds\n"
		"\t\t-l     <str>   Input the merger list produced by regenerate\n"
		"\t\t-o     <str>   OutPut file\n"
		"\n"
		"\t\t-h             show this help\n"
		"\n";
	return 1;
}


int parse_AcmdFABP (int argc, char **argv , In3str1v  * paraFABP)
{
	if (argc <=2 ) {print_AusageFABP();return 0;}

	for(int i = 1; i < argc ;  i++)
	{
		if(argv[i][0] != '-' )
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InPut" ||  flag  == "i" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFABP->InStr1=argv[i];
		}
		else if (flag  ==  "MerList"   ||  flag  == "l")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFABP->InStr2=argv[i];
		}
		else if (flag  ==  "OutPut"  ||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFABP->InStr3=argv[i];
		}
		else if (flag  == "help"  ||  flag  == "h")
		{
			print_AusageFABP();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if  ((paraFABP->InStr1).empty() ||  (paraFABP->InStr2).empty()  ||  (paraFABP->InStr3).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	(paraFABP->InStr3)=add_Asuffix(paraFABP->InStr3) ;
	return 1 ;
}

int FA_ChangPosi_main(int argc, char *argv[])
	//	int main(int argc, char *argv[])
{
	In3str1v  * paraFABP = new In3str1v ;
	if( parse_AcmdFABP(argc, argv, paraFABP )==0)
	{
		delete  paraFABP ;
		return 1 ;
	}

	map <string , map <llong , pair <string,llong > > >  Name_sca ;

	igzstream LIST ((paraFABP->InStr2).c_str(),ifstream::in);
	igzstream SNP ((paraFABP->InStr1).c_str(),ifstream::in);
	ogzstream OUT ((paraFABP->InStr3).c_str());

	if(!LIST.good())
	{
		cerr << "open IN File error: "<<(paraFABP->InStr2)<<endl;
		return 1;
	}
	if(!SNP.good())
	{
		cerr << "open IN File error: "<<(paraFABP->InStr1)<<endl;
		return 1;
	}
	if(!OUT.good())
	{
		cerr << "open OUT File error: "<<(paraFABP->InStr3)<<endl;
		return 1;
	}


	map <string , map <llong , pair <string,llong > > > :: iterator  map_it ;

	while(!LIST.eof())
	{
		string  line ;
		getline(LIST,line);
		if (line.length()<=0)  { continue  ; }
		if (line[0] == '#' )  { continue  ; }
		string sca_name , chr_name;
		llong sca_start , sca_end , chr_start, chr_end;
		istringstream isone (line,istringstream::in);
		isone>>sca_name>>sca_start>>sca_end>>chr_name>>chr_start>>chr_end ;
		map_it=Name_sca.find(chr_name);
		pair <string,llong > Pare_data = make_pair(sca_name,chr_end);

		if (map_it==Name_sca.end())
		{
			map <llong ,pair <string,llong > > Site2Sca;
			Site2Sca[chr_start]=Pare_data ;
			Name_sca.insert(map <string , map <llong , pair <string,llong > > >  :: value_type (chr_name, Site2Sca ));
		}
		else
		{
			(map_it->second).insert(map <llong , pair <string,llong > >  :: value_type (chr_start, Pare_data));
		}
	}
	LIST.close();

	map <llong , pair <string,llong > > :: iterator Sed_Map_it ;
	while(!SNP.eof())
	{
		string  line ;
		getline(SNP,line);
		if (line.length()<=0)  { continue  ;}
		if (line[0] == '#' )  {OUT<<line<<"\n"; continue  ; }
		vector<string> inf;
		split(line,inf,"\t ");
		string chr_name=inf[0];
		map_it=Name_sca.find(chr_name);
		if (map_it==Name_sca.end())
		{
			cerr<<"something wrong at this chr \t"<<line<<endl;
		}
		else
		{
			llong Site  ; 
			llong newSite=-1 ;
			string  newName;
			istringstream isone2 (inf[1],istringstream::in);
		 	isone2>>Site ;
			Sed_Map_it=(map_it->second).begin();
			for ( Sed_Map_it=(map_it->second).begin() ; Sed_Map_it!=(map_it->second).end(); Sed_Map_it++ )
			{
				if (Site >= (Sed_Map_it->first)  &&   Site<=((Sed_Map_it->second).second)  )
				{
					newSite=Site-(Sed_Map_it->first)+1;
					newName=((Sed_Map_it->second).first);
					break ;
				}
			}

			if  (newSite == -1 )
			{
				cerr<<"something wrong at this site\t"<<line<<endl;
			}
			else
			{
				OUT<<newName<<"\t"<<newSite ;
				int vet_size=inf.size();
				for (int j=2; j< vet_size ; j++)
				{
					OUT<<"\t"<<inf[j];
				}
				OUT<<"\n";
			}
		}
	}
	SNP.close();
	OUT.close();


	delete paraFABP ;
	return 0;
}

///////// swimming in the sky and flying in the sea ////////////
#endif //  ChangPosi_h_
