#ifndef RefReform_H_
#define RefReform_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include "../ALL/comm.h"
#include "../ALL/DataClass.h"
#include <cstdlib>
#include <gzstream.h>
#include <zlib.h>
#include <stdio.h>
#include "../ALL/kseq.h"

//KSEQ_AINIT(gzFile, gzread)

using namespace std;

int  print_AusageFA06()
{
	cout <<""
		"\n"
		"\tUsage: reform -i <in.fa> \n"
		"\n"
		"\t\t-i           <str>   Input fa for upper/lower case\n"
		"\n"
		"\t\t-o     <str>   OutPut file [STDOUT]\n"
		"\t\t-d             dele remove the comment of each seq\n"
		"\t\t-s     <str>   upper/lower the seq base (upper/lower)\n"
		"\t\t-r             reverse the sequence\n"
		"\t\t-c             complement the sequence\n"
		"\t\t-a             one lines one seq to stored\n"
		"\t\t-e     <int>   How X bp to cut end one line seq\n"
		"\n"
		"\t\t-h             show this help\n"
		"\n";
	return 1;
}


int parse_AcmdFA06(int argc, char **argv , In3str1v  * paraFA06 , int & cutline )
{
	if (argc <=2 ) {print_AusageFA06();return 0;}

	for(int i = 1; i < argc ;  i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InPut" ||    flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA06->InStr1=argv[i];
		}
		else if (flag  == "cutline"  ||  flag  == "e" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			cutline=atoi(argv[i]);
		}
		else if (flag  == "Case" ||   flag  == "s")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA06->InStr3=argv[i];
			if ( ( (paraFA06->InStr3) != "upper"  )   &&   ( (paraFA06->InStr3) != "lower"  )   )
			{
				cerr<<"-Case shoule be  [upper]  or  [lower] " <<endl;
				return 0;
			}
			//			paraFA06->TF=false;
		}
		else if (flag  ==  "NOcomment"  ||   flag  ==  "d")
		{	
			paraFA06->InInt=(paraFA06->InInt) | 0x4;
		}
		else if (flag  ==  "Reverse" || flag  ==  "r" )
		{
			paraFA06->InInt=(paraFA06->InInt) | 0x1;
		}
		else if (flag  ==  "Complement"  || flag  ==  "c")
		{
			paraFA06->InInt=(paraFA06->InInt) | 0x2;
		}
		else if (flag  ==  "oneline" ||   flag  == "a")
		{
			paraFA06->InInt=(paraFA06->InInt) | 0x8;
		}
		else if (flag  ==  "OutPut" || flag  == "o" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA06->InStr2=argv[i];
		}
		else if (flag  == "help" || flag  == "h")
		{
			print_AusageFA06();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if  ((paraFA06->InStr1).empty())
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	if (!(paraFA06->InStr2).empty() )
	{
		(paraFA06->InStr2)=add_Asuffix(paraFA06->InStr2) ;
	}
	return 1 ;
}

int FA_Reform_main(int argc, char *argv[])
//int main(int argc, char *argv[])
{
	In3str1v  * paraFA06 = new In3str1v ;
	int cutline=-1 ;
	if( parse_AcmdFA06(argc, argv, paraFA06 , cutline)==0)
	{
		delete  paraFA06;
		return 1;
	}

	char buf[65536];
	setbuf(stdout, buf);

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);
	//map <char ,char >  Complement;

	char Complement[256]={'N'};	
	Complement['A']='T'; Complement['C']='G';  Complement['T']='A'; Complement['G']='C'; Complement['N']='N';
	Complement['a']='t'; Complement['c']='g';  Complement['t']='a'; Complement['g']='c'; Complement['n']='n';
	Complement['M']='K'; Complement['R']='Y';  Complement['W']='W'; Complement['S']='S'; Complement['Y']='R';
	Complement['K']='M';
	Complement['m']='k'; Complement['r']='y';  Complement['w']='w'; Complement['s']='s'; Complement['y']='r';
	Complement['k']='m';
	Complement['\n']='\n';



	if  ((paraFA06->InInt)==3)
	{
		igzstream  INID ((paraFA06->InStr1).c_str(),ifstream::in);
		if(!INID.good())
		{
			cerr << "open InputFile error: "<<(paraFA06->InStr1)<<endl;
			delete paraFA06 ; return 0;
		}

		string seq;
		getline(INID, seq ,'>');

		if ((paraFA06->InStr2).empty())
		{
			while(!INID.eof())
			{
				string line ;
				getline(INID, line,'\n');
				if (line.empty())
				{
					continue;
				}
				getline(INID, seq ,'>');
				int END=seq.length()-1;	
///*///
				int START=0;
				while(START<END)
				{
					while(seq[START]=='\n')
					{
						START++;
					}
					while(seq[END]=='\n')
					{
						END--;
					}
					char Atmp=Complement[seq[END]];
					seq[END]=Complement[seq[START]];
					seq[START]=Atmp;
					START++;
					END--;
				}
				if (END==START)
				{
					seq[END]=Complement[seq[START]];
				}
				fprintf(stdout,">%s\n%s",(line).c_str(),seq.c_str());
//				*///
			}
		}
		else
		{
			ogzstream  OUT2 ((paraFA06->InStr2).c_str());
			while(!INID.eof())
			{
				string line;
				getline(INID, line,'\n');
				if (line.empty())
				{
					continue;
				}
				getline(INID, seq ,'>');
				int END=seq.length()-1;
				int START=0;
				while(START<END)
				{
					while(seq[START]=='\n')
					{
						START++;
					}
					while(seq[END]=='\n')
					{
						END--;
					}
					char Atmp=Complement[seq[END]];
					seq[END]=Complement[seq[START]];
					seq[START]=Atmp;
					START++;
					END--;
				}
				if (END==START)
				{
					seq[END]=Complement[seq[START]];
				}
				OUT2<<">"<<line<<"\n"<<seq;
			}


			OUT2.close();

		}
		INID.close();


		delete paraFA06 ;
		return 0;
	}
	else if ((paraFA06->InInt)==2)
	{
		igzstream  INID ((paraFA06->InStr1).c_str(),ifstream::in);
		if(!INID.good())
		{
			cerr << "open InputFile error: "<<(paraFA06->InStr1)<<endl;
			delete paraFA06 ; return 0;
		}
		string seq;
		getline(INID, seq ,'>');

		if ((paraFA06->InStr2).empty())
		{
			while(!INID.eof())
			{
				string line;
				getline(INID, line,'\n');
				if (line.empty())
				{
					continue;
				}
				getline(INID, seq ,'>');
				int Refleng=seq.length();
				for (int i=0 ; i<Refleng ; i++)
				{
					seq[i]=Complement[seq[i]];
				}
				fprintf(stdout,">%s\n%s",(line).c_str(),seq.c_str());
			}
		}
		else
		{
			ogzstream  OUT2 ((paraFA06->InStr2).c_str());
			while(!INID.eof())
			{
				string line;
				getline(INID, line,'\n');
				if (line.empty())
				{
					continue;
				}
				getline(INID, seq ,'>');
				int Refleng=seq.length();
				for (int i=0 ; i<Refleng ; i++)
				{
					seq[i]=Complement[seq[i]];
				}
				OUT2<<">"<<line<<"\n"<<seq;
			}
			OUT2.close();
		}


		INID.close();
		delete paraFA06 ;
		return 0;
	}
	else if ((paraFA06->InInt)==1)
	{

		igzstream  INID ((paraFA06->InStr1).c_str(),ifstream::in);
		if(!INID.good())
		{
			cerr << "open InputFile error: "<<(paraFA06->InStr1)<<endl;
			delete paraFA06 ; return 0;
		}

		string seq;
		getline(INID, seq ,'>');

		if ((paraFA06->InStr2).empty())
		{
			while(!INID.eof())
			{
				string line,ID;
				getline(INID, line,'\n');
				if (line.empty())
				{
					continue;
				}
				getline(INID, seq ,'>');
				int END=seq.length()-1;
				int START=0;
				while(START<END)
				{
					while(seq[START]=='\n')
					{
						START++;
					}
					while(seq[END]=='\n')
					{
						END--;
					}
					char Atmp=seq[END];
					seq[END]=seq[START];
					seq[START]=Atmp;
					START++;
					END--;
				}
				fprintf(stdout,">%s\n%s",(line).c_str(),seq.c_str());
			}


		}
		else
		{

			ogzstream  OUT2 ((paraFA06->InStr2).c_str());
			while(!INID.eof())
			{
				string line;
				getline(INID, line,'\n');
				if (line.empty())
				{
					continue;
				}
				getline(INID, seq ,'>');
				int END=seq.length()-1;
				int START=0;
				while(START<END)
				{
					while(seq[START]=='\n')
					{
						START++;
					}
					while(seq[END]=='\n')
					{
						END--;
					}
					char Atmp=seq[END];
					seq[END]=seq[START];
					seq[START]=Atmp;
					START++;
					END--;
				}
				OUT2<<">"<<line<<"\n"<<seq;
			}

			OUT2.close();

		}
		INID.close();
		delete paraFA06 ;
		return 0;
	}


	int linecut=cutline;
	if (cutline<0)
	{
		linecut	= FaCutLine ((paraFA06->InStr1));
	}

	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen((paraFA06->InStr1).c_str(), "r");
	seq = kseq_init(fp);
	ogzstream  OUT ;
	if (!(paraFA06->InStr2).empty())
	{
		OUT.open((paraFA06->InStr2).c_str());
	}

	while ((l = kseq_read(seq)) >= 0)
	{
		string Ref=(seq->seq.s) ;
		string chr=(seq->name.s);
		string ID=chr ;
		int Refleng=(seq->seq.l);

		if (seq->comment.l)
		{
			string comment=(seq->comment.s);
			ID+="\t"+comment;
		}

		if  ((paraFA06->InInt) & 0x4)
		{
			ID=chr ;
		}

		if ((paraFA06->InStr3).empty())
		{
		}
		else if ((paraFA06->InStr3)== "upper" )
		{
			transform(Ref.begin(),Ref.end(),Ref.begin(),::toupper) ;
		}
		else if ((paraFA06->InStr3)== "lower" )
		{
			transform(Ref.begin(),Ref.end(),Ref.begin(), ::tolower) ;
		}

		if  ( ((paraFA06->InInt) & 0x1)   &&  ((paraFA06->InInt) & 0x2) )
		{
			int Min=Refleng/2;
			if (Refleng%2!=0)
			{
				Ref[Min]=Complement[Ref[Min]];
			}
			int CCDD=Refleng-1;
			for (int i=0 ; i<Min ; i++)
			{
				char N=Complement[Ref[i]];
				Ref[i]=Complement[Ref[CCDD-i]]; Ref[CCDD-i]=N;
			}
		}
		else if ( (paraFA06->InInt) & 0x1 )
		{
			reverse(Ref.begin(), Ref.end());
		}
		else if  ( (paraFA06->InInt) & 0x2  )
		{
			for (int i=0 ; i<Refleng ; i++)
			{
				Ref[i]=Complement[Ref[i]];
			}
		}
		if  ((paraFA06->InInt) & 0x8)
		{
			if ((paraFA06->InStr2).empty())
			{
				//				cout<<">"<<ID<<endl;				cout<<Ref<<endl;
				fprintf(stdout,">%s\n%s\n",ID.c_str(),Ref.c_str());
			}
			else
			{
				OUT<<">"<<ID<<"\n";
				OUT<<Ref<<"\n";
			}
		}
		else
		{
			int Endline=int(Refleng/linecut);
			if ((paraFA06->InStr2).empty())
			{
				//cout<<">"<<ID<<endl;

				fprintf(stdout,">%s\n",ID.c_str());
				for (int i=0 ; i< Endline ; i++)
				{
					//	cout<<Ref.substr(i*linecut,linecut)<<"\n";
					fprintf(stdout,"%s\n",(Ref.substr(i*linecut,linecut)).c_str());
				}
				long  AA=Endline*linecut ;

				if (Refleng > AA)
				{
					//					cout<<Ref.substr(AA)<<endl;
					fprintf(stdout,"%s\n",(Ref.substr(AA)).c_str());
					//
				}
			}
			else
			{
				Display( Ref , ID ,  OUT , linecut  ) ;
			}
		}
	}

	if (!(paraFA06->InStr2).empty())
	{
		OUT.close();
	}

	kseq_destroy(seq);
	gzclose(fp);
	delete paraFA06 ;
	return 0;
}

///////// swimming in the sky and flying in the sea ////////////
#endif // RefReform_H_
