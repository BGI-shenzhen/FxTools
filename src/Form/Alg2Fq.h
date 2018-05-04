#include <iostream>
#include <fstream>
#include <gzstream.h>
#include "../Soap/xam2soap/TSamCtrl.h"
#include "../Soap/xam2soap/Ttools.h"
#include "../ALL/comm.h"

using namespace std;


int  print_Ausage_A32()
{
	cout <<""
		"\n"
		"\tUsage: Alg2Fq  -InBam <in.bam>  -OutPut <Out.fq>\n"
		"\n"
		"\t\t-InSoap    <str>   Input Soap Format file\n"
		"\t\t-InSam     <str>   Input Sam Format file\n"
		"\t\t-InBam     <str>   Input Bam Format file\n"
		"\t\t-OutPut    <str>   OutPut the Fq Format file\n"
		"\n"
		"\t\t-UnMap             only output UnMapp read [NA]\n"
		"\t\t-SamNoHead         If InFile is Sam with NoHead\n"
		"\t\t-help              show this help\n" 
		"\n";
	return 1;
}


int parse_Acmd_A32(int argc, char **argv, ParaClass * para_A32)
{
	if (argc <=2 ) {print_Ausage_A32();return 0 ;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InSam" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_A32->InPut1=argv[i];
			(para_A32->InInt1)=1;
		}
		else if (flag  == "InBam" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_A32->InPut2=argv[i];
			para_A32->InInt1=2;
		}
		else if (flag  == "InSoap" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_A32->InStr3=argv[i];
			para_A32->InInt2=4;
		}
		else if (flag  ==  "OutPut")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_A32->OutPut1=argv[i];
		}
		else if (flag  == "SamNoHead")
		{
			para_A32->InInt1+=32;
		}
		else if (flag  == "UnMap")
		{
			para_A32->TF=false ;
		}
		else if (flag  == "help")
		{
			print_Ausage_A32();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if  ((para_A32->OutPut1).empty() || ( (para_A32->InInt1)==0 &&  (para_A32->InInt2)==0 )  )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	(para_A32->OutPut1)= add_Asuffix ( (para_A32->OutPut1) ) ;

	return 1 ;
}

void sam2fq( std::vector<std::string> vec , ogzstream & OUT , map <char,char> Samp )
{            
	int flag = atoi(vec[1].c_str());
	string outprint="";
	if (flag & (0x1 << 6)) {
		outprint = "@"+vec[0] + "/1" + "\n"; 
	} else if (flag & (0x1 << 7)) {
		outprint = "@"+vec[0] + "/2" + "\n"; 
	} else {
		outprint = "@"+vec[0] + "\n";
	}
	if ( flag & 16 )
	{
		reverse(vec[9].begin(), vec[9].end());        
		reverse(vec[10].begin(), vec[10].end());
		int leng=vec[9].size();
		for (int i=0 ; i<leng ; i++)
		{
			vec[9][i]=Samp[vec[9][i]];
		}
	}
	outprint=outprint+vec[9]+"\n+\n"+vec[10];
	OUT<<outprint<<endl;
}

int Xam2fq_main(int argc, char **argv)
	//int main(int argc, char **argv)
{

	ParaClass * para_A32 = new ParaClass;
	if (parse_Acmd_A32(argc, argv,para_A32  )==0)
	{
		delete  para_A32 ;
		return  0;
	}
	ogzstream OUT(para_A32->OutPut1.c_str()) ;
	if (!OUT.good())
	{
		cerr<<"Can't open OutFile "<<para_A32->OutPut1<<endl;
		return 0;
	}

	map <char,char> Samp;
	Samp['A']='T'; Samp['C']='G';
	Samp['T']='A'; Samp['G']='C';
	Samp['N']='N';

	char in_mode[5] ={ 0 };
	string sampath ;
	//char * sampath=argv[2];

	if ((para_A32->InInt2)==4)
	{
		igzstream INSOAP ((para_A32->InStr3).c_str(),ifstream::in);
		if (!INSOAP.good())
		{
			cerr<<"Can't open OutFile "<<para_A32->InStr3<<endl;
			delete para_A32 ; return  0;
		}

		while(!INSOAP.eof())
		{
			string  line ;
			getline(INSOAP,line);
			if (line.length()<=0)  { continue  ; }
			istringstream isone (line,istringstream::in);
			string ID , Seq ,Qseq , ab  ,zf ,printOUTFQ ;
			int hit , Rlength ;
			isone>>ID>>Seq>>Qseq>>hit>>ab>>Rlength>>zf ;
			int ID_length=ID.length()-2 ;
			if( ID_length>0 && ID[ID_length]=='/' && (ID[ID_length+1]=='1' || ID[ID_length+1]=='2' ) )
			{
				printOUTFQ="@"+ID+"\n";
			}
			else if (ab == "a")
			{
				printOUTFQ="@"+ID+"/1\n";
			}
			else 
			{
				printOUTFQ="@"+ID+"/2\n";
			}

			if (zf == "-")
			{
				reverse(Seq.begin(), Seq.end());        
				reverse(Qseq.begin(), Qseq.end());
				for (int i=0 ; i<Rlength ; i++)
				{
					Seq[i]=Samp[Seq[i]];
				}
			}
			printOUTFQ=printOUTFQ+Seq+"\n+\n"+Qseq;
			OUT<<printOUTFQ<<endl;
		}
		INSOAP.close();
	}

	if ( (para_A32->InInt1) >16)
	{

		if   ( (para_A32->InPut1).empty() )
		{
			cerr<<"since -SamNoHead should together with -InSam "<<endl;
			OUT.close();
			delete para_A32 ;
			return 0;
		}

		igzstream INSAM ((para_A32->InPut1).c_str(),ifstream::in);
		if (!INSAM.good())
		{
			cerr<<"Can't open OutFile "<<para_A32->InPut1<<endl;
			delete para_A32 ; return  0 ;
		}

		if (!(para_A32->TF))
		{
			while(!INSAM.eof())
			{
				string  line ;
				getline(INSAM,line);
				if (line.length()<=0)  { continue  ; }
				std::vector<std::string> vec;
				TStringSplit(line, '\t', vec);
				if ((vec[5] == "*"  && vec[8]=="0" ) || (vec.size()<13 ))
				{
					sam2fq( vec , OUT  , Samp );
				}
			}

		}
		else
		{
			while(!INSAM.eof())
			{
				string  line ;
				getline(INSAM,line);
				if (line.length()<=0)  { continue  ; }
				std::vector<std::string> vec;
				TStringSplit(line, '\t', vec);
				sam2fq( vec , OUT  , Samp );
			}
		}
		INSAM.close();
		OUT.close();
		delete para_A32 ;
		return 0;
	}


	if ((para_A32->InInt1)==1)
	{
		in_mode[0]='r';
		sampath=(para_A32->InPut1);        
	}
	else if  ((para_A32->InInt1)==2 )
	{
		in_mode[0]='r';in_mode[1]='b';
		sampath=(para_A32->InPut2);
	}
	else
	{
		delete  para_A32 ;
		OUT.close();
		return  0;
	}

	TSamCtrl samhandle;
	samhandle.open(sampath.c_str(),in_mode);
	string line;
	if (!(para_A32->TF))
	{
		while(samhandle.readline(line)!=-1)
		{
			std::vector<std::string> vec;
			TStringSplit(line, '\t', vec);
			if ((vec[5] == "*"  && vec[8]=="0" ) || (vec.size()<13 ))
			{
				sam2fq( vec , OUT  , Samp );
			}
		}
	}
	else
	{
		while(samhandle.readline(line)!=-1)
		{
			std::vector<std::string> vec;
			TStringSplit(line, '\t', vec);
			sam2fq( vec , OUT , Samp );
		}
	}

	OUT.close();
	delete para_A32 ;
	return 0;
}
///////// swimming in the sky and flying in the sea ////////////
