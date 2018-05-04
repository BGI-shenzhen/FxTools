#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <zlib.h>
#include "./src/ALL/gzstream.C"
#include "./src/ALL/kseq.h"
#include "./src/ALL/comm.h"
#include "./src/Form/FormTools.h"
#include "./src/FaDeal/FaTools.h"
#include "./src/Fq/FqTools.h"
#include "./src/ALL/HelpTemp.h"


using namespace std;

int Form_Tools_main(int argc, char *argv[]) ;
int Fa_Tools_main(int argc, char *argv[]) ;
int FQ_Tools_main(int argc, char *argv[]) ;
int help_main ( int argc , char *argv[]  ) ;

static int  AllTools_usage ()
{
	cout <<"Program: FxTools\nVersion: 0.16\thewm2008@gmail.com\t"<<__DATE__<<endl;

	cout<<""
		"\n"
		"\tUsage:\n\n"
		"\t\tFatools        Tools For Fasta\n"
		"\t\tFqtools        Tools For Fastq\n"
		"\t\tFormtools      Tools For Form convert\n"
		"\n"
		"\t\tHelp           Show help in detail\n"
		"\n";
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) { return AllTools_usage(); }
	else if (strcmp(argv[1], "Fatools") == 0) { return Fa_Tools_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Fqtools") == 0) { return FQ_Tools_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Formtools") == 0) { return Form_Tools_main(argc-1, argv+1) ; }
	else if (strcmp(argv[1], "Help") == 0 || strcmp(argv[1], "help") == 0 || strcmp(argv[1], "?")== 0 || ( argv[1][0] == '-' &&( argv[1][1] =='h' || argv[1][1] =='H' || argv[1][1] =='?' ) )  || strcmp(argv[1], "less") == 0 )
	{
		return help_main (argc , argv );
	}
	else
	{
		cerr<<"FxTools [main] unrecognized command "<<argv[1]<<endl;
		return 1;
	}

return 0;
}


