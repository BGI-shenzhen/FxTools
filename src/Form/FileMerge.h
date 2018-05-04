#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <cstdlib>
#include "../ALL/comm.h"
#include "../ALL/DataClass.h"
#include <gzstream.h>
using namespace std;

///////////////////

int  print_Other03()
{
	cout <<""
		"\n"
		"\tUsage: merge -i <soap.list> -o <out.soap> -c 9\n"
		"\n"
		"\t\t-i    <str>   Input sort SubFile list\n"
		"\t\t-o    <str>   Output Sort file for merge\n"
		"\t\t\n"
		"\t\t-c    <int>   specification of the column used for merging[9]\n"
//		"\t\t-Row       <int>   the row Number of the sort[9]\n"
		"\t\t-h            show this help\n" 
		"\n";
	return 1;
}

int parse_Other03(int argc, char **argv , In3str1v * para_Other03 )
{
	if (argc <=2 ) {print_Other03();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InList" ||  flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_Other03->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut" ||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_Other03->InStr2=argv[i];
		}
		else if (flag  == "Row"  || flag  == "By" ||  flag  == "c")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			(para_Other03->InInt)=atoi(argv[i]);
		}
		else if (flag  == "help" ||  flag  == "h")
		{
			print_Other03();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}

	}
	if  ((para_Other03->InStr1).empty() || ((para_Other03->InStr2).empty()))
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	(para_Other03->InStr2)=add_Asuffix(para_Other03->InStr2) ;
	return 1 ;
}
////////////////////////////

void  run_Aread (igzstream & IN , long & position , ogzstream & OUT ,string & head_Aone, int Raw=9 )
	//void  run_Aread ( igzstream & IN , long & position , ofstream & OUT ,string & head_Aone )
{
	string line ,temp ;
	long posi_Anow ;
	OUT<<head_Aone<<endl ;

	for (int i =1 ; i>0 ; i++)
	{
		if (getline(IN,line))
		{
			istringstream isone (line,istringstream::in);
			for (int kk=1 ; kk<Raw ; kk++ )
			{
				isone>>temp ;
			}
			isone>>posi_Anow ;

			if (posi_Anow<=position)
			{
				OUT<<line<<endl ;
			}
			else
			{
				head_Aone=line; 
				position=posi_Anow ;
				break ;
			}
		}
		else
		{
			IN.close();
			position=-2;
			break ;
		}
	}    
}


int chose_Afile ( const long Posi[] , const int count )
{
	long start=Posi[0] ;
	int min_Aflag=0;
	for ( int i=1 ; i <count  ; i++ )
	{
		if (  start <0 )
		{
			min_Aflag=i ;
			start=Posi[i];
		}
		else if ( Posi[i]>-1 &&  Posi[i]<start)
		{
			min_Aflag=i ;
			start=Posi[i];
		}
	}
	if  (Posi[min_Aflag]<0 )
	{
		min_Aflag=-2;
	}
	return  min_Aflag ;    
}



int Other_Merge_main(int argc,char *argv[])
{
	In3str1v  * para_Other03 = new In3str1v ;
	(para_Other03->InInt)=9;
	if (parse_Other03( argc, argv, para_Other03)==0)
	{
		delete para_Other03 ;
		return 1;
	}

	ogzstream OUT ((para_Other03->InStr2).c_str());
	if(!OUT.good())
	{
		cerr << "open OutFile error: "<<(para_Other03->InStr2)<<endl;
		return 1;
	}

	ReadList ( (para_Other03->InStr1), (para_Other03->List));
	int File_count=(para_Other03->List).size();

	igzstream *Sam = new  igzstream[File_count] ;
	//      long Posi[File_count];
	long *Posi  = new  long[File_count];
	string temp ,line ;
	long min_now  ;
	//      string Pring[File_count];
	string *Pring = new string[File_count];
	//      string head1 ,head2 ;
	for (int i=0; i<File_count ; i++)
	{
		Sam[i].open((para_Other03->List)[i].c_str(),ifstream::in) ;
		if  (Sam[i].good())
		{
			cout<<"sortfile\t"<<(para_Other03->List)[i]<<"\tsort Raw\t"<<(para_Other03->InInt)<<endl;
			getline(Sam[i],line);
			Pring[i]=line ;
			istringstream isone (line,istringstream::in);
			for (int kk=1 ; kk<(para_Other03->InInt); kk++ )
			{
				isone>>temp ;
			}
			isone>>min_now;
			Posi[i]=min_now;
		}
		else
		{
			cerr<<(para_Other03->List)[i]<<"\tcan't open"<<endl ;
			Pring[i]=line ;
			Posi[i]=-2 ;
		}
	}
	//        OUT<<head2<<endl;
	int file_run_cout=chose_Afile(Posi , File_count);
	while( file_run_cout > -1 )
	{        
		run_Aread(Sam[file_run_cout],Posi[file_run_cout],OUT,Pring[file_run_cout] , (para_Other03->InInt) );
		file_run_cout=chose_Afile ( Posi, File_count ) ;
	}
	delete [] Sam ; delete [] Posi  ; delete [] Pring ;

	OUT.close();
	delete para_Other03 ;
	return 0 ;
}



////////////////////////swimming in the sea & flying in the sky //////////////////



//////////////////
////////////////////////swimming in the sea & flying in the sky //////////////////

///////// swimming in the sky and flying in the sea ////////////
