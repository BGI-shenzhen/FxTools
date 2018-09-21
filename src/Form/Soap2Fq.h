#include <iostream>
#include <fstream>
#include "../ALL/gzstream.h"
#include "../ALL/comm.h"

using namespace std;


int  print_Ausage_A32()
{
	cout <<""
		"\n"
		"\tUsage: Soap2Fq  -i <in.soap>  -o <Out.fq>\n"
		"\n"
		"\t\t-i    <str>   input Soap\n"
		"\t\t-o    <str>   output FASTQ\n"
		"\n"
		"\t\t-u            only output UnMapp read\n"
		"\t\t-h            show more details for help\n" 
		"\n";
	return 1;
}

void More_HelpFM13()
{
	cout<<"\n\
\n\
\t\t1.  Soap2fq -i <in.bam > -o AAA\n\
\t\tThis will convert all the SOAP to FASTQ and output to a compressed file named AAA in current directory\n\
\n\
\t\t2.  Soap2fq -i <in.bam > -o AAA -u\n\
\t\tThis will convert the unmapped part of SOAP to FASTQ and output to a compressed file named AAA in current directory. If all the reads are mapped, this will give nothing.\n\
\n\
\n";
}


int parse_Acmd_A32(int argc, char **argv, ParaClass * para_A32)
{
	if (argc <=1 ) {print_Ausage_A32();return 0 ;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InSoap"  ||  flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_A32->InStr3=argv[i];
		}
		else if (flag  ==  "OutPut"||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_A32->OutPut1=argv[i];
		}
		else if (flag  == "UnMap" ||  flag  == "u")
		{
			para_A32->TF=false ;
		}
		else if (flag  == "help" ||  flag  == "h")
		{
		   More_HelpFM13();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if  ((para_A32->OutPut1).empty() ||  (para_A32->InStr3).empty()  )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	(para_A32->OutPut1)= add_Asuffix ( (para_A32->OutPut1) ) ;

	return 1 ;
}

int Soap2fq_main(int argc, char **argv)
//int main(int argc, char **argv)
{
	ParaClass * para_A32 = new ParaClass;
	if (parse_Acmd_A32(argc, argv,para_A32  )==0)
	{
		delete  para_A32 ;
		return 1;
	}
	ogzstream OUT(para_A32->OutPut1.c_str()) ;
	if (!OUT.good())
	{
		cerr<<"Can't open OutFile "<<para_A32->OutPut1<<endl;
		return 1;
	}

	int Samp[256]={0};
//	map <char,char> Samp;
	Samp['A']='T'; Samp['C']='G';
	Samp['T']='A'; Samp['G']='C';
	Samp['N']='N';


	igzstream INSOAP ((para_A32->InStr3).c_str(),ifstream::in);
	if (!INSOAP.good())
	{
		cerr<<"Can't open OutFile "<<para_A32->InStr3<<endl;
		delete para_A32 ; return  1;
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



	OUT.close();
	delete para_A32 ;
	return 0;
}
///////// swimming in the sky and flying in the sea ////////////
