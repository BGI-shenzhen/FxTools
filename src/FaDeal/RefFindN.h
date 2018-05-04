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
#include <algorithm>
#include <zlib.h>
#include <stdio.h>
#include "../ALL/kseq.h"

using namespace std;

int  print_AusageA29()
{
	cout <<""
		"\n"
		"\tUsage: findN  -i <in.fa> \n"
		"\n"
		"\t\t-i    <str>   InPut fa for find N Region\n"
		"\n"
		"\t\t-o    <str>   OutPut file,otherwise[STDOUT]\n"
		"\n"
		"\t\t-h            show this help\n" 
		"\n";
	return 1;
}


int parse_AcmdA29(int argc, char **argv , In3str1v * paraA29 )
{
	if (argc <=2 ) {print_AusageA29();return 0;}

	for(int i = 1; i < argc  ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InPut" ||  flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg( flag ) ; return 0;}
			i++;
			paraA29->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut"  ||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg( flag ) ;return 0;}
			i++;
			paraA29->InStr2=argv[i];
		}
		else if (flag  == "help" ||  flag  == "h")
		{
			print_AusageA29();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((paraA29->InStr1).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	return 1 ;
}




int FA_NRegion_main(int argc, char *argv[])
{
	In3str1v * paraA29 = new In3str1v;
	if( parse_AcmdA29(argc, argv, paraA29 )==0)
	{
		delete  paraA29 ;
		return 0;
	}
	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen((paraA29->InStr1).c_str(), "r");
	seq = kseq_init(fp);

	std::ostream * OUT = &cout;
	std::ofstream fout;
	if (!(paraA29->InStr2).empty())
	{
		fout.open((paraA29->InStr2).c_str());
		OUT = &fout;
	}

	while ((l = kseq_read(seq)) >= 0) 
	{
		string Ref=(seq->seq.s) ;
		size_t N_tail_position = 0 ;
		size_t N_head_position=Ref.find_first_of( "Nn", 0 );
		while( N_head_position != string::npos)
		{
			N_tail_position=Ref.find_first_of("AaTtCcGg",N_head_position+1);

			if (N_tail_position!= string::npos )
			{
				(*OUT)<<(seq->name.s)<<"\t"<<N_head_position+1<<"\t"<<N_tail_position<<"\n";
			}
			else
			{
				(*OUT)<<(seq->name.s)<<"\t"<<N_head_position+1<<"\t"<<(seq->seq.l)<<"\n";
				break ;
			}
			N_head_position=Ref.find_first_of("Nn",N_tail_position);
		}
	} 

	kseq_destroy(seq);
	gzclose(fp);
	delete paraA29 ;
	return 0;

}

///////// swimming in the sky and flying in the sea ////////////
