#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "../ALL/comm.h"
#include <cstdlib>
#include <gzstream.h>
#include <algorithm>
#include <zlib.h>
#include <stdio.h>
//#include "../ALL/kseq.h"
#include "../ALL/DataClass.h"
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h> 
#include <fcntl.h> 
#include <unistd.h>
//#include <strong>

using namespace std;
typedef long long  llong ;

int  usage_sort()
{
	cout <<""
		"\n"
		"\tUsage: sort -i  <in.fa> -o  <out.fa>\n"
		"\n"
		"\t\t-i     <str>   Input fa File to sort\n"
		"\t\t-o     <str>   Output sort fa File,otherwise[STDOUT]\n"
		"\n"
		"\t\t-s     <str>   Rank the Seq by [name/length]\n"
		"\t\t               Default [name]\n"
		"\t\t-d             descending sort. default is ascending sort\n"
		//"\t\t-reverse         reverse the result,no seq\n"
		"\t\t-h             Show this help\n"
		"\n";
	return 1;
}


int parse_Acmd_FA18(int argc, char **argv ,In3str1v * para_FA18 )
{
	if (argc <=2 ) {usage_sort();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InPut" || flag  == "i" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_FA18->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut"   || flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_FA18->InStr2=argv[i];
		}       
		else if (flag  == "Sort" ||  flag  == "s")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			string temp=(argv[i]);
			if (temp == "length" || temp == "Length" )
			{
				para_FA18->TF=false;
			}
		}
		else if (flag  == "reverse"  ||  flag  == "descend" ||  flag  == "d")
		{
			para_FA18->TF2=false;
		}
		else if (flag  == "help" || flag  == "h")
		{
			usage_sort();return 0;
		}

		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ( (para_FA18->InStr1).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	return 1 ;
}


class FaStrut
{
	public:
		string ID;
		long Length;
		long Start ;
		long End ;
};


//int main(int argc, char *argv[])
int FA_Sort_main(int argc, char *argv[])
{
	In3str1v * para_FA18 = new In3str1v;

	if( parse_Acmd_FA18(argc, argv, para_FA18 )==0)
	{
		delete  para_FA18 ;
		return 1;
	}

	//char buf[65536];
	//setbuf(stdout, buf);

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);
	char * buffer;

	/*//
	  FILE * pFile;
	  long FileSize;
	  pFile = fopen ((para_FA18->InStr1).c_str(), "rb" );  
	  if (pFile==NULL)  
	  {  
	  fputs ("File error",stderr);  
	  exit (1);  
	  }  

	  fseek (pFile , 0 , SEEK_END);  
	  FileSize = ftell (pFile);  

	  rewind (pFile); 
	  buffer = (char*) malloc (sizeof(char)*FileSize);  
	  if (buffer == NULL)  
	  {  
	  fputs ("Memory error",stderr);   
	  exit (2);  
	  }  

	  size_t result;  
	  result = fread (buffer,1,FileSize,pFile);  
	  if (result != FileSize)  
	  {  
	  fputs ("Reading error",stderr);  
	  exit (3);  
	  }  

	  fclose (pFile);  
	  *///

	string	path=para_FA18->InStr1;
	long FileSize;
	int fd;
	string extA =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);
	if (extA!="gz")
	{
		fd = open((para_FA18->InStr1).c_str(),O_RDONLY);
		FileSize = lseek(fd,0,SEEK_END);
		buffer = (char *) mmap(NULL,FileSize,PROT_READ,MAP_PRIVATE,fd,0);
	}
	else
	{
		igzstream  INID ((para_FA18->InStr1).c_str(),ifstream::in);
		if(!INID.good())
		{
			cerr << "open InputFile error: "<<(para_FA18->InStr1)<<endl;
			delete para_FA18 ; return 1;
		}
		string chrALL;
		getline(INID, chrALL ,'%');
		FileSize=chrALL.length();
		buffer=(char*) malloc (sizeof(char)*FileSize);
		strcpy(buffer,chrALL.c_str());
		INID.close();
		chrALL="";
	}

	list <FaStrut>  ChrSeq;
	//vector <FaStrut>  ChrSeq;
	FaStrut  This ;
	This.Start=0;
	This.ID="";
	long ii=1;

	while(buffer[ii]!='\n')
	{
		if (buffer[ii]!=' ')
		{
			This.ID=This.ID+buffer[ii];
		}
		ii++;
	}

	char *EE = strchr(buffer+ii,'>');

	while(EE!=NULL)
	{
		ii=EE-buffer;
		This.End=ii;
		This.Length=ii-This.Start;
		ChrSeq.push_back(This);
		This.Start=This.End;
		ii++;
		This.ID="";
		while(buffer[ii]!='\n')
		{
			if (buffer[ii]!=' ')
			{
				This.ID=This.ID+buffer[ii];
			}
			ii++;
		}
		EE = strchr(buffer+ii,'>');
	}

	This.End=FileSize;
	This.Length=FileSize-This.Start ;
	ChrSeq.push_back(This);

	if ((para_FA18->TF) && (para_FA18->TF2))
	{
		ChrSeq.sort([](FaStrut AA,FaStrut BB){return AA.ID>BB.ID;});
	}
	else if (((para_FA18->TF)) && (!(para_FA18->TF2)))
	{
		ChrSeq.sort([](FaStrut AA,FaStrut BB){return AA.ID<BB.ID;});
	}
	else if ((!(para_FA18->TF)) && ((para_FA18->TF2)))
	{
		ChrSeq.sort([](FaStrut AA,FaStrut BB){return AA.Length>BB.Length;});
	}
	else
	{
		ChrSeq.sort([](FaStrut AA,FaStrut BB){return AA.Length<BB.Length;});
	}

	//vector <FaStrut>::iterator it;

	list <FaStrut>::iterator it;

	if  ((para_FA18->InStr2).empty())
	{
		for(it = ChrSeq.begin();it!=ChrSeq.end();it++)
		{
			fwrite(buffer+(it->Start),1,it->Length,stdout);
			//	fprintf(stdout,">%s\n%s",(it->line).c_str(),(Seq[(it->Num)]).c_str());
		}	
	}
	///*/////
	else
	{
		path=para_FA18->InStr2;
		string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);
		if  (ext!="gz")
		{

			int fdV;
			fdV=open((para_FA18->InStr2).c_str(),O_RDWR,0766);
			if(fdV==-1)
			{
				fdV=open((para_FA18->InStr2).c_str(),O_RDWR|O_CREAT,0766);
				if(fdV==-1)
				{ 
					printf("can not open outFile\n");
					exit(-1);
				}
			}

			ftruncate(fdV,FileSize);

			char* p=(char *)(mmap(NULL,FileSize,PROT_READ|PROT_WRITE,MAP_SHARED,fdV,0)); 
			if(p==NULL || MAP_FAILED==p)
			{
				printf("can Not Wrtie\n");
				exit(-1);
			}

			long long shift =0 ;
			for(it = ChrSeq.begin();it!=ChrSeq.end();it++)
			{
				for (long jj=(it->Start); jj<(it->End); jj++)
				{
					p[shift]=buffer[jj];
					shift++;
				}
			}
			munmap(p,FileSize);
			close(fdV);
		}
		else		
		{


			char *cmd;
			cmd = (char*)malloc((para_FA18->InStr2).length() + 20);
			sprintf(cmd, " gzip >  %s ", (para_FA18->InStr2).c_str());
			FILE *out = popen(cmd, "w");
			free(cmd);

			for(it = ChrSeq.begin();it!=ChrSeq.end();it++)
			{
				fwrite(buffer+(it->Start),1,it->Length,out);
			}
			pclose(out);

			/*////
			  gzFile OUTGZ;
			  OUTGZ = gzopen ((para_FA18->InStr2).c_str(), "wb");
			  for(it = ChrSeq.begin();it!=ChrSeq.end();it++)
			  {
			  gzwrite(OUTGZ,buffer+(it->Start),it->Length);
			  }
			  gzclose(OUTGZ);
			//*/
			/*
			   ogzstream OUT ((para_FA18->InStr2).c_str());
			   for(it = ChrSeq.begin();it!=ChrSeq.end();it++)
			   {
			   for (long jj=(it->Start); jj<(it->End); jj++)
			   {
			   OUT.put(buffer[jj]);
			   }
			   }
			   OUT.close();
			//*/
		}
	}


	if (extA=="gz")
	{
		free (buffer);
	}
	else
	{
		close(fd);
	}
	delete para_FA18 ;
	return 0;
}

///////// swimming in the sky and flying in the sea ////////////
