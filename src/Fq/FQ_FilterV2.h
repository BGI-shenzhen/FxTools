#ifndef FQ_FilterV1_H_
#define FQ_FilterV1_H_ 


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <cstdlib>
#include <algorithm>

#include <gzstream.h>
#include "../ALL/comm.h"
#include "FQ_Filter.h"

typedef long long  llong ;
using namespace std;
//using namespace __gnu_cxx ;


int  print_usage_A_B24 ()
{
	cout <<""
		"\n"
		"Usage:filter  -i A1.fq  A2.fq -o B1.fq  B2.fq [options]\n"
		"\n"
		"\t\t-i      <str>     File name of InFq1/InFq2 Input\n"
		"\t\t-o      <str>     File name of OutFq1/OutFq2.gz output\n"
		"\n"
		"\t\t-n      <float>   Cut N reach X of read length[0.1]\n"
		"\t\t-Q      <int>     defined the Low Quality base [5]\n"
		"\t\t-d      <float>   dele LowQ base reach X of read length[0.5]\n"
		"\t\t-f                Cut off the index\n"
		"\t\t-s                Trim the Bad terminal continuous Base[NA]\n"
		"\t\t-l      <int>     Rm too short Read length after trim [30]\n"
		"\t\t-h                show this help\n" 
		"\n";
	return 1;
}
//
int parse_cmd_A_B24(int argc, char **argv , Para_A24 * para )
{
	if (argc <=4  ) {print_usage_A_B24();return 0;}

	int err_flag = 0;
	for(int i = 1; i < argc || err_flag; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InFq1" || flag  == "i")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->InFq1=argv[i];
			if ( (i+1) < argc  && (argv[i+1][0]!='-'))
			{
				i++;
				para->InFq2=argv[i];
			}
		}
		else if (flag  ==  "InFq2")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->InFq2=argv[i];
		}
		else if (flag  ==  "OutFq1"  || flag  == "o")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->OutFq1=argv[i];
			if ( (i+1) < argc  && (argv[i+1][0]!='-'))
			{
				i++;
				para->OutFq2=argv[i];
			}
		}
		else if (flag  ==  "OutFq2")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ;return 0;}
			i++;
			para->OutFq2=argv[i];
		}
		else if (flag  ==  "IndexCut" ||  flag  == "f" )
		{
			para->indexcut=1;
		}
		else if (flag  ==  "trim" || flag  ==  "s" )
		{
			para->adapter1="trim";
		}
		else if (flag  ==  "MinLeng" ||  flag  ==  "l" )
		{
			if(i + 1 == argc) {  LogLackArg( flag ) ; return 0;}
			i++;
			para->minLeng=atoi(argv[i]);
		}
		else if (flag  ==  "OffLowQ" ||  flag  ==  "d")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->OffLowQ=atof(argv[i]);
		}
		else if (flag  ==  "LowQ"  || flag  ==  "Q" )
		{
			if(i + 1 == argc) {  LogLackArg( flag ) ; return 0;}
			i++;
			para->LowQint=atoi(argv[i]);
		}
		else if (flag  == "help"   ||  flag  ==  "h" )
		{
			print_usage_A_B24();return 0;
		}
		else if (flag  == "OffN"  ||  flag  == "n" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ;return 0;}
			i++;
			para->offN=atof(argv[i]);
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if ( (para->OutFq1).empty() || (para->InFq1).empty() )
	{
		cerr<< "-InFq1 -OutFq1 lack argument for the must"<<endl;
		return 0;
	}
	(para->OutFq1)=add_Asuffix(para->OutFq1);

	if  ((para->OutFq2).empty() &&  (para->InFq2).empty())
	{
		return 1 ;
	}
	else if ( (para->OutFq2).empty() || (para->InFq2).empty() )
	{
		cerr<< "-InFq2 -OutFq2 lack argument for the must"<<endl;
		return 0;
	}
	else
	{
		(para->OutFq2)=add_Asuffix(para->OutFq2);
		return  2 ;
	}
}


int No_trim ( string & Seq , string & Quli ,  Para_A24 * para )
{
	string::size_type N_seq_length= 0 ;
	int LOWQ_base_Number=0;
	for(string::size_type ix=0 ; ix<(para->maxLeng) ; ix++)
	{
		if ( Seq[ix]==(para->N_seq) )
		{
			N_seq_length++;
		}
		if (Quli[ix]<=(para->LowQint))
		{
			LOWQ_base_Number++;
		}
	}

	if ( ((LOWQ_base_Number*1.0)/(para->maxLeng)) >= (para->OffLowQ) )
	{
		return  1 ;
	}
	if ((N_seq_length*1.0/(para->maxLeng))>=(para->offN))
	{
		return 1 ;
	}

	return  0 ;         
}

////////////////////////////

int FQSE_FilterNotrim( Para_A24 * para )
{
	int if_run_index=0 ;

	map <string,int> Good_index  ;
	vector <pair <string,int> > Adapter_one ;

	if ((para->indexcut)==1)
	{
		if_run_index=whether_index(para->indexcut,para->InFq1,Good_index);
	}

	DA *BaseD =new DA[(para->maxLeng)*2];

	igzstream IN_1 ((para->InFq1).c_str(),ifstream::in); // ifstream  + gz 
	ogzstream OUT_1 ((para->OutFq1).c_str());

	if(!IN_1.good())
	{
		cerr << "open IN File error: "<<para->InFq1<<endl;
		return 1;
	}
	if(!OUT_1.good())
	{
		cerr << "open OUT File error: "<<para->OutFq1<<endl;
		return 1;
	} 


	int UC[128]={4};
	UC['A']=0; UC['C']=1;UC['G']=3;UC['T']=2 ; UC['N']=4;
	UC['a']=0; UC['c']=1;UC['g']=3;UC['t']=2 ; UC['n']=4;


	string ID_1 ,seq_1,temp_1,Quly_1 ;
	long long read_Fisrt=0 ,Base_Fisrt=0 , read_Second=0 , Base_Second=0;
	bool No_Trim=true ;
	if ((para->adapter1)=="trim")
	{
		No_Trim=false ;
	}



	if (if_run_index !=0 )
	{
		string::size_type pos_two,pos_one ,leng_cut ; 
		string index_sub ;
		while(!IN_1.eof())
		{
			getline(IN_1,ID_1);
			getline(IN_1,seq_1);
			getline(IN_1,temp_1);
			getline(IN_1,Quly_1);

			if (ID_1.empty())
			{
				continue ;
			}
			pos_one=ID_1.rfind('#')+1;
			pos_two=ID_1.rfind('/');
			leng_cut=pos_two-pos_one;
			index_sub=ID_1.substr(pos_one,leng_cut);

			if (!Good_index.count(index_sub))
			{
				continue ;
			}
			read_Fisrt+=1 ; Base_Fisrt+=seq_1.size();
			if (No_Trim)
			{
				if (No_trim ( seq_1,Quly_1,  para )==1)
				{
					continue ;
				}
			}
			else
			{
				if (trim ( seq_1,Quly_1,  para )==1)
				{
					continue ;
				}
			}
			int seq_length=seq_1.size();
			read_Second+=1 ; Base_Second+=seq_length; 

			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
			if  (seq_length==(para->maxLeng))
			{
				StortD ( seq_1 , BaseD ,para , UC );
			}
		}
	}
	else
	{
		while(!IN_1.eof())
		{
			getline(IN_1,ID_1);
			getline(IN_1,seq_1);
			getline(IN_1,temp_1);
			getline(IN_1,Quly_1);

			if (ID_1.empty())
			{
				continue ;
			}

			read_Fisrt+=1 ; Base_Fisrt+=seq_1.size();

			if (No_Trim)
			{
				if (No_trim ( seq_1,Quly_1,  para )==1)
				{
					continue ;
				}
			}
			else
			{
				if (trim ( seq_1,Quly_1,  para )==1)
				{
					continue ;
				}
			}

			int seq_length=seq_1.size();
			read_Second+=1 ; Base_Second+=seq_length; 
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
			if  (seq_length==(para->maxLeng))
			{
				StortD ( seq_1 ,  BaseD ,para ,UC);
			}
		}

	}


	IN_1.close();
	OUT_1.close();

	read_Fisrt=(read_Fisrt);
	read_Second=(read_Second);
	cout<<"##Final Stat OUT##"<<endl;
	cout<<"#Original Read\t"<<read_Fisrt<<"\t#After Filter Read\t"<<read_Second<<endl;
	cout<<"#Original Base\t"<<Base_Fisrt<<"\t#After Filter Base\t"<<Base_Second<<endl;
	int ALentmp=((para->maxLeng));
	for (int ii=0 ; ii<ALentmp; ii++)
	{
		int Sum=(BaseD[ii]).get(0)+(BaseD[ii]).get(2)+(BaseD[ii]).get(4)+(BaseD[ii]).get(1)+(BaseD[ii]).get(3);

		if (Sum<1)
		{
			continue ;
		}
		cout<<ii+1<<"\t"<<(BaseD[ii]).get(0)<<"\t"<<(BaseD[ii]).get(1)<<"\t"<<(BaseD[ii]).get(2)<<"\t"<<(BaseD[ii]).get(3)<<"\t"<<(BaseD[ii]).get(4)<<"\t"<<"\n";
	}
	delete [] BaseD ;
}
//programme entry
///////// swimming in the sky and flying in the sea ////////////
int FQ_FilterNotrim_main(int argc, char **argv)
{
	Para_A24 * para = new Para_A24;
	int Flag_para=parse_cmd_A_B24(argc, argv ,para ) ;
	if ( Flag_para ==0)
	{
		delete  para ; 
		return 1 ;
	}
	int tmpShiftQ=GetShiftQ((para->InFq1)) ;
	(para->LowQint)	=  (para->LowQint)+tmpShiftQ ;

	igzstream INFQ ((para->InFq1).c_str(),ifstream::in); // ifstream  + gz
	if(!INFQ.good())
	{
		cerr << "open IN File error: "<<para->InFq1<<endl;
		return 1;
	}
	para->maxLeng=0 ;
	int BB=0 ;
	while( (!INFQ.eof()) && ( BB< 16888 ))
	{
		string  line ;
		getline(INFQ,line);
		if (line.length()<=0)  { continue  ; }
		getline(INFQ,line);
		int tmp=line.length();
		if (tmp> (para->maxLeng) )
		{
			(para->maxLeng)=tmp;
		}
		getline(INFQ,line);
		getline(INFQ,line);
		BB++;
	}
	INFQ.close();


	if ( Flag_para ==1)
	{
		FQSE_FilterNotrim (para );
		delete para;
		return  2;
	}
	int if_run_index=0 ;
	map <string,int> Good_index  ;

	if ((para->indexcut)==1)
	{
		if_run_index=whether_index(para->indexcut,para->InFq1,Good_index);
	}


	DA *BaseD =new DA[(para->maxLeng)*2];

	igzstream IN_1 ((para->InFq1).c_str(),ifstream::in); // ifstream  + gz 
	igzstream IN_2 ((para->InFq2).c_str(),ifstream::in); // ifstream  + gz 
	ogzstream OUT_1 ((para->OutFq1).c_str());
	ogzstream OUT_2 ((para->OutFq2).c_str());

	if(!IN_1.good())
	{
		cerr << "open IN File error: "<<para->InFq1<<endl;
		return 1;
	}
	if(!IN_2.good())
	{
		cerr << "open IN File error: "<<para->InFq2<<endl;
		return 1;
	}
	if(!OUT_1.good())
	{
		cerr << "open OUT File error: "<<para->OutFq1<<endl;
		return 1;
	} 
	if(!OUT_2.good())
	{
		cerr << "open OUT File error: "<<para->OutFq2<<endl;
		return 1;
	}



	int UC[128]={4};
	UC['A']=0; UC['C']=1;UC['G']=3;UC['T']=2 ; UC['N']=4;
	UC['a']=0; UC['c']=1;UC['g']=3;UC['t']=2 ; UC['n']=4;




	string ID_1 ,seq_1,temp_1,Quly_1 ;
	string ID_2 ,seq_2,temp_2,Quly_2 ;
	long long read_Fisrt=0 , read_Second=0 ,Base_Fisrt=0 ,Base_Second =0;

	bool No_Trim=true ;
	if ((para->adapter1)=="trim")
	{
		No_Trim=false ;
	}

	if (if_run_index !=0 )
	{
		string::size_type pos_two,pos_one ,leng_cut ; 
		string index_sub ;
		while(!IN_1.eof())
		{
			getline(IN_1,ID_1);
			getline(IN_1,seq_1);
			getline(IN_1,temp_1);
			getline(IN_1,Quly_1);
			getline(IN_2,ID_2);
			getline(IN_2,seq_2);
			getline(IN_2,temp_2);
			getline(IN_2,Quly_2);
			pos_one=ID_1.rfind('#')+1;
			pos_two=ID_1.rfind('/');
			leng_cut=pos_two-pos_one;
			index_sub=ID_1.substr(pos_one,leng_cut);
			if (! Good_index.count(index_sub))
			{
				continue ;
			}            
			read_Fisrt+=1 ; Base_Fisrt+=seq_1.size();  Base_Fisrt+=seq_2.size();
			if (No_Trim)
			{
				if (No_trim (seq_1,Quly_1,  para  )==1)
				{
					continue ;
				}
				if (No_trim (seq_2,Quly_2, para  )==1)
				{
					continue ;            
				}
			}
			else
			{
				if (trim (seq_1,Quly_1,  para  )==1)
				{
					continue ;
				}
				if (trim (seq_2,Quly_2, para  )==1)
				{
					continue ;            
				}
			}
			int Read_1_length=seq_1.size(); 
			int Read_2_length=seq_2.size(); 

			read_Second+=1 ; Base_Second+=Read_1_length; Base_Second+=Read_2_length;
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
			OUT_2<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n" ;
			if  ( Read_2_length==(para->maxLeng)  &&  Read_1_length==(para->maxLeng) )
			{
				StortD ( seq_1 ,  seq_2 , BaseD ,para , UC);
			}
		}
	}
	else
	{
		while(!IN_1.eof())
		{
			getline(IN_1,ID_1);
			getline(IN_1,seq_1);
			getline(IN_1,temp_1);
			getline(IN_1,Quly_1);

			getline(IN_2,ID_2);
			getline(IN_2,seq_2);
			getline(IN_2,temp_2);
			getline(IN_2,Quly_2);
			if (ID_1.empty())
			{
				continue ;
			}

			read_Fisrt+=1 ; Base_Fisrt+=seq_1.size();  Base_Fisrt+=seq_2.size();


			if (No_Trim)
			{
				if (No_trim (seq_1,Quly_1,  para  )==1)
				{
					continue ;
				}
				if (No_trim (seq_2,Quly_2, para  )==1)
				{
					continue ;            
				}
			}
			else
			{
				if (trim (seq_1,Quly_1,  para  )==1)
				{
					continue ;
				}
				if (trim (seq_2,Quly_2, para  )==1)
				{
					continue ;            
				}
			}
			int Read_1_length=seq_1.size(); 
			int Read_2_length=seq_2.size(); 

			read_Second+=1 ; Base_Second+=Read_1_length; Base_Second+=Read_2_length;
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
			OUT_2<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n" ;
			if  ( Read_2_length==(para->maxLeng)  &&  Read_1_length==(para->maxLeng) )
			{
				StortD ( seq_1 ,  seq_2 , BaseD ,para , UC );
			}

		}
	}
	IN_2.close();
	OUT_2.close();
	IN_1.close();
	OUT_1.close();

	read_Fisrt=(read_Fisrt)*2 ;
	read_Second=(read_Second)*2;
	cout<<"##Final Stat OUT##"<<endl;
	cout<<"#Original Read\t"<<read_Fisrt<<"\t#After Filter Read\t"<<read_Second<<endl;
	cout<<"#Original Base\t"<<Base_Fisrt<<"\t#After Filter Base\t"<<Base_Second<<endl;
	int ALentmp=(2*(para->maxLeng));
	for (int ii=0 ; ii<ALentmp; ii++)
	{
		int Sum=(BaseD[ii]).get(0)+(BaseD[ii]).get(2)+(BaseD[ii]).get(4)+(BaseD[ii]).get(1)+(BaseD[ii]).get(3);
		if (Sum<1)
		{
			continue ;
		}
		cout<<ii+1<<"\t"<<(BaseD[ii]).get(0)<<"\t"<<(BaseD[ii]).get(1)<<"\t"<<(BaseD[ii]).get(2)<<"\t"<<(BaseD[ii]).get(3)<<"\t"<<(BaseD[ii]).get(4)<<"\t"<<"\n";
	}
	delete para ;
	delete [] BaseD ;
	return 0 ;
}

#endif  // FQ_FilterV1_H_


