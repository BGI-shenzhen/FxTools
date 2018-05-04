#ifndef RQReform_H_
#define RQReform_H_

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
#include <algorithm>
#include "../ALL/kseq.h"

//KSEQ_AINIT(gzFile, gzread)

using namespace std;

int  print_AusageFQ_reverse()
{
	cout <<""
		"\n"
		"\tUsage: reform -i <in.fq> \n"
		"\n"
		"\t\t-i    <str>   Input fa for upper/lower case\n"
		"\n"
		"\t\t-o    <str>   OutPut file [STDOUT]\n"
		"\t\t-s    <str>   upper/lower the seq base (upper/lower)\n"
		"\t\t-r            reverse the sequence\n"
		"\t\t-c            complement the sequence\n"
		"\n"
		"\t\t-h            show this help [hewm2008]\n"
		"\n";
	return 1;
}


int parse_AcmdFQ_reverse(int argc, char **argv , In3str1v  * paraFQ_reverse  )
{
	if (argc <=2 ) {print_AusageFQ_reverse();return 0;}

	for(int i = 1; i < argc ;  i++)
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
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFQ_reverse->InStr1=argv[i];
		}
		else if (flag  == "Case" ||  flag  == "s" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFQ_reverse->InStr3=argv[i];
			if ( ( (paraFQ_reverse->InStr3) != "upper"  )   &&   ( (paraFQ_reverse->InStr3) != "lower"  ) )
			{
				cerr<<"-Case shoule be  [upper]  or  [lower]" <<endl;
				return 0;
			}
			//			paraFQ_reverse->TF=false;
		}
		else if (flag  ==  "Reverse" ||  flag  == "r")
		{
			paraFQ_reverse->TF2=false;
			//InInt=(paraFQ_reverse->InInt) | 0x1;
		}
		else if (flag  ==  "Complement" ||  flag  == "c")
		{
			paraFQ_reverse->TF=false;
			//(paraFQ_reverse->InInt) | 0x2;
		}
		else if (flag  ==  "OutPut" ||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFQ_reverse->InStr2=argv[i];
		}
		else if (flag  == "help" ||  flag  == "h")
		{
			print_AusageFQ_reverse();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if  ((paraFQ_reverse->InStr1).empty())
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	if (!(paraFQ_reverse->InStr2).empty() )
	{
		(paraFQ_reverse->InStr2)=add_Asuffix(paraFQ_reverse->InStr2) ;
	}
	return 1 ;
}

int FQ_Reform_main(int argc, char *argv[])
//int main(int argc, char *argv[])
{
	In3str1v  * paraFQ_reverse = new In3str1v ;
	if( parse_AcmdFQ_reverse(argc, argv, paraFQ_reverse)==0)
	{
		delete  paraFQ_reverse;
		return 1 ;
	}

	//map <char ,char >  Complement;
	char Complement[256]={'N'};
	Complement['A']='T'; Complement['C']='G';  Complement['T']='A'; Complement['G']='C'; Complement['N']='N';
	Complement['a']='t'; Complement['c']='g';  Complement['t']='a'; Complement['g']='c'; Complement['n']='n';
	Complement['M']='K'; Complement['R']='Y';  Complement['W']='W'; Complement['S']='S'; Complement['Y']='R';
	Complement['K']='M';
	Complement['m']='k'; Complement['r']='y';  Complement['w']='w'; Complement['s']='s'; Complement['y']='r';
	Complement['k']='m';

	
	std::ios::sync_with_stdio(false);
	std::cin.tie(0);


	igzstream IN ((paraFQ_reverse->InStr1).c_str(),ifstream::in); // ifstream  + gz

	if(!IN.good())
	{
		cerr << "open IN File error: "<<paraFQ_reverse->InStr1<<endl;
		return 1;
	}

	string ID_1 ,Ref,temp_1,Quly_1 ;

	ogzstream  OUT ;
	if (!(paraFQ_reverse->InStr2).empty())
	{
		OUT.open((paraFQ_reverse->InStr2).c_str());
		while(!IN.eof())
		{
			getline(IN,ID_1);
			getline(IN,Ref);
			getline(IN,temp_1);
			getline(IN,Quly_1);
			int Refleng=Ref.length();
			if (Refleng<=0)  { continue  ; }

			if ((paraFQ_reverse->InStr3)== "upper" )
			{
				transform(Ref.begin(),Ref.end(),Ref.begin(),::toupper) ;
			}
			else if ((paraFQ_reverse->InStr3)== "lower" )
			{
				transform(Ref.begin(),Ref.end(),Ref.begin(), ::tolower) ;
			}
			if ( !(paraFQ_reverse->TF2) )
			{
				reverse(Ref.begin(), Ref.end());
				reverse(Quly_1.begin(), Quly_1.end());
			}
			if  (!(paraFQ_reverse->TF))
			{
				for (int i=0 ; i<Refleng ; i++)
				{
					Ref[i]=Complement[Ref[i]];
				}
			}
			OUT<<ID_1<<"\n"<<Ref<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n";
		}
		OUT.close();
	}
	else
	{
		while(!IN.eof())
		{
			getline(IN,ID_1);
			getline(IN,Ref);
			getline(IN,temp_1);
			getline(IN,Quly_1);
			int Refleng=Ref.length();
			if (Refleng<=0)  { continue ;}
			
			if ((paraFQ_reverse->InStr3).empty())
			{
			}
			else if ((paraFQ_reverse->InStr3)== "upper"  ||  (paraFQ_reverse->InStr3)== "Upper" )
			{
				transform(Ref.begin(),Ref.end(),Ref.begin(),::toupper) ;
			}
			else if ((paraFQ_reverse->InStr3)== "lower"  ||  (paraFQ_reverse->InStr3)== "Lower" )
			{
				transform(Ref.begin(),Ref.end(),Ref.begin(), ::tolower) ;
			}

			if ( !(paraFQ_reverse->TF2) )
			{
				reverse(Ref.begin(), Ref.end());
				reverse(Quly_1.begin(), Quly_1.end());
			}

			if  (!(paraFQ_reverse->TF))
			{
				for (int i=0 ; i<Refleng ; i++)
				{
					Ref[i]=Complement[Ref[i]];
				}
			}
			fprintf(stdout,"%s\n%s\n%s\n%s\n",ID_1.c_str(),Ref.c_str(),temp_1.c_str(),Quly_1.c_str());
		}
	}

	IN.close();
	delete paraFQ_reverse ;
	return 0;
}

///////// swimming in the sky and flying in the sea ////////////
#endif // RQReform_H_
