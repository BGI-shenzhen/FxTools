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
#include "../ALL/gzstream.h"
using namespace std;

///////////////////

int  print_Other03()
{
	cout <<""
		"\n"
		"\tUsage: Merge -i <file.list> -o <out.sort > -c <int>\n"
		"\n"
		"\t\t-i    <str>   input list of sorting file\n"
		"\t\t-o    <str>   output sorted merged result\n"
		"\t\t\n"
		"\t\t-c    <int>   pecification of one or two columns used for sorting the merged file.[9]\n"
		"\t\t-h            show more details for help\n" 
		"\n";
	return 1;
}

void More_HelpFM334()
{
	cout<<"\n\
\n\
\t\t1.	merge -i file.list -c X -o AAA\n\
\t\tThis will merge the files in the file.list, sort by column X and output the result to a compressed file AAA in current directory.\n\
\t\t(1.1)	The input files should also be sorted by column X, otherwise the output file simply captures the input files together without sorting.\n\
\t\t(1.2)	The values in column X of input files should be numeric.\n\
\n\
\t\t2.	merge -i file.list -c X,Y -o AAA\n\
\t\tThis will merge the files in the file.list, sort first by column X and then by column Y and output the result to a compressed file AAA in current directory.\n\
\t\t(2.1) The input files should also be sorted by column X and Y, otherwise the output file simply captures the input files together without sorting.\n\
\t\t(2.2) The values in column X and Y of input files should be string and numeric respectively.\n\
\n\
\n";
}

int parse_Other03(int argc, char **argv , In3str1v * para_Other03 )
{
	if (argc <=1 ) {print_Other03();return 0;}

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
		else if (flag  == "Row"  || flag  == "f" ||  flag  == "c")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			(para_Other03->InStr3)=argv[i];
		}
		else if (flag  == "help" ||  flag  == "h")
		{
			More_HelpFM334();return 0;
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
{
	string line ,temp ;
	long posi_Anow ;
	OUT<<head_Aone<<'\n' ;
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
				OUT<<line<<'\n';
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



void  run_Aread (igzstream & IN , long & position , string & chrID , ogzstream & OUT ,string & head_Aone, int RawFist ,  int  RawSecond )
{
	string chrID_This;
	long posi_This;
	OUT<<head_Aone<<'\n';

	vector<string> TempV;
	for (int i =1 ; i>0 ; i++)
	{
		if (getline(IN,head_Aone))
		{
			TempV.clear();
			split(head_Aone,TempV," \t");
			chrID_This=TempV[RawFist];  posi_This=atoi(TempV[RawSecond].c_str());

			if (posi_This<=position  &&  chrID_This  == chrID)
			{
				OUT<<head_Aone<<'\n';
			}
			else
			{
				position=posi_This ;
				chrID=chrID_This;
				break ;
			}
		}
		else
		{
			chrID="*";
			position=-8;
			break ;
		}
	}

}

int chose_Bfile ( const long Posi[] , string  ChrID[] ,  const int count)
{
	string IDD ="*";
	int i=0;

	int start = -2;

	for (   ; i <count  ; i++)
	{
		if (ChrID[i]!="*")
		{
			IDD=ChrID[i];
			break;
		}
	}

	for (   ; i <count  ; i++)
	{
		if (ChrID[i]!="*")
		{
			if (IDD>ChrID[i])
			{
				IDD=ChrID[i];
			}
		}
	}

	//cerr<<IDD<<"\t"<<i<<"\t";
	int min_Aflag=-1;
	for ( i=0 ; i <count  ; i++)
	{
		if  ( ChrID[i]!=IDD )
		{
			continue ;
		}
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
//	cout<<ChrID[0]<<"\t"<<ChrID[1]<<"\t"<<Posi[0]<<"\t"<<Posi[1]<<"\t#\t"<<min_Aflag<<endl;
	if  (Posi[min_Aflag]<-7)
	{
		min_Aflag=-2;
	}

	return  min_Aflag ;
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
	(para_Other03->InStr3)="9";
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
	long *Posi  = new  long[File_count];
	string temp ,line ;
	long min_now  ;
	string *Pring = new string[File_count];

	if  ((para_Other03->InStr3).find(",")==string::npos)
	{
		(para_Other03->InInt)=atoi((para_Other03->InStr3).c_str());
		for (int i=0; i<File_count ; i++)
		{
			Sam[i].open((para_Other03->List)[i].c_str(),ifstream::in) ;
			if (Sam[i].good())
			{
				cout<<"sortfile\t"<<(para_Other03->List)[i]<<"\tsort Raw\t"<<(para_Other03->InInt)<<endl;
				getline(Sam[i],line);
				while(line[0]=='#')
				{
					OUT<<line<<endl;
					getline(Sam[i],line);
				}
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
		int file_run_cout=chose_Afile(Posi , File_count);
		while( file_run_cout > -1 )
		{        
			run_Aread(Sam[file_run_cout],Posi[file_run_cout],OUT,Pring[file_run_cout] , (para_Other03->InInt) );
			file_run_cout=chose_Afile ( Posi, File_count ) ;
		}
	}
	else
	{

		vector <string> TempV;
		split((para_Other03->InStr3),TempV,",");

		int InInt2=atoi(TempV[1].c_str())-1;
		(para_Other03->InInt)=atoi(TempV[0].c_str())-1;		
		string *ChrID = new string[File_count];
		for (int i=0; i<File_count ; i++)
		{
			Sam[i].open((para_Other03->List)[i].c_str(),ifstream::in) ;
			if  (Sam[i].good())
			{
				cout<<"sortfile\t"<<(para_Other03->List)[i]<<"\tsort Raw\t"<<(para_Other03->InInt)<<","<<InInt2<<endl;
				getline(Sam[i],line);
				while(line[0]=='#')
				{
					OUT<<line<<endl;
					getline(Sam[i],line);
				}
				Pring[i]=line;

				TempV.clear();
				split(line,TempV," \t");
				ChrID[i]=TempV[(para_Other03->InInt)];  Posi[i]=atoi(TempV[InInt2].c_str());
			}
			else
			{
				cerr<<(para_Other03->List)[i]<<"\tcan't open"<<endl ;
				Pring[i]=line ;
				Posi[i]=-8 ;
				ChrID[i]="NANA";
			}
		}

		int file_run_cout=chose_Bfile(Posi,ChrID, File_count);

		while( file_run_cout > -1 )
		{        
			run_Aread(Sam[file_run_cout],Posi[file_run_cout], ChrID[file_run_cout], OUT,Pring[file_run_cout] , (para_Other03->InInt),InInt2);
			file_run_cout=chose_Bfile ( Posi,ChrID, File_count) ;
		}

		delete [] ChrID;
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
