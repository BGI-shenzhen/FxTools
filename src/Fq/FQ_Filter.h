#ifndef FQ_Filter_H_
#define FQ_Filter_H_ 

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
#include "../ALL/DataClass.h"

typedef long long  llong ;
using namespace std;
//using namespace __gnu_cxx ;

class Para_A24 {
	public:
		string InFq1;
		string InFq2;
		string OutFq1;
		string OutFq2;

		string adapter1 ;
		string adapter2 ;
		char minBaseQ;
		char  N_seq;
		float offN;
		float OffLowQ ;

		int indexcut ;
		int minLeng ;
		int LowQint ;
		unsigned int maxLeng ;

		Para_A24()
		{
			InFq1="";
			InFq2="";
			OutFq1="";
			OutFq2="";

			adapter1="NA" ;
			adapter2="NA" ;
			minBaseQ='B' ;
			N_seq='N';
			LowQint=5;

			offN=0.1 ;
			OffLowQ=0.5;
			indexcut=0 ;
			minLeng=30 ;
			maxLeng=100;
		}
} ;
 
int  print_usage_A24 ()
{
	cout <<""
		"\n"
		"Usage:filterV1  -InFq1 A1.fq -InFq2 A2.fq -OutFq1 B1.fq -OutFq2 B2.fq [options]\n"
		"\n"
		"\t\t-InFq1        <str>   File name of InFq1 Input\n"
		"\t\t-OutFq1       <str>   File name of OutFq1.gz output\n"
		"\t\t-InFq2        <str>   File name of InFq2 Input\n"
		"\t\t-OutFq2       <str>   File name of OutFq2.gz output\n"
		"\n"
		"\t\t-MinLeng      <int>   The min  Read length [30]\n"
		"\t\t-MaxLeng      <int>   The max  Read length [100]\n"
		"\t\t-IndexCut             Cut off the index \n"
		"\t\t-OffN       <float>   Cut N reach X of read length[0.1]\n"
		"\t\t-help                 show this help\n" 
		"\n";
	return 1;
}
//
int parse_cmd_A24(int argc, char **argv , Para_A24 * para )
{
	if (argc <=4  ) {print_usage_A24();return 0;}

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

		if (flag  == "InFq1" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->InFq1=argv[i];
		}
		else if (flag  ==  "InFq2")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->InFq2=argv[i];
		}
		else if (flag  ==  "OutFq1")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->OutFq1=argv[i];
		}
		else if (flag  ==  "OutFq2")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ;return 0;}
			i++;
			para->OutFq2=argv[i];
		}
		else if (flag  ==  "Adapter1")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->adapter1=argv[i];
		}
		else if (flag  ==  "Adapter2")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ;return 0;}
			i++;
			para->adapter2=argv[i];
		}
		else if (flag  ==  "IndexCut")
		{
			para->indexcut=1;
		}

		else if (flag  ==  "MaxLeng")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			para->maxLeng=atoi(argv[i]);
		}
		else if (flag  ==  "MinLeng" )
		{
			if(i + 1 == argc) {  LogLackArg( flag ) ; return 0;}
			i++;
			para->minLeng=atoi(argv[i]);
		}
		else if (flag  == "help")
		{
			print_usage_A24();return 0;
		}
		else if (flag  == "OffN")
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


int trim ( string & Seq ,  string & Quli ,  Para_A24 * para  )
{
	string A=Seq ;
	string B=Quli ;
	string::size_type i_end=A.size()-1;        
	string::size_type s_start=0 ,s_end =i_end ;
	for(string::size_type ix=0 ; ix<=i_end ; ix++)
	{
		if ( A[ix]==(para->N_seq) || B[ix]==(para->minBaseQ) )
		{
			s_start=ix;
			continue ;
		}
		else
		{
			s_start=ix;
			break  ;
		}
	}


	for(string::size_type ix=s_end ; ix>s_start ; --ix)
	{
		if ( A[ix]==(para->N_seq) || B[ix]==(para->minBaseQ) )
		{
			s_end=ix;
			continue ;
		}
		else
		{
			s_end=ix;
			break ;
		}
	}
	int  cut_length=s_end-s_start+1 ;
	if (cut_length<(para->minLeng))
	{
		return  1 ; 
	}
	string::size_type N_seq_length= 0 ;
	for(string::size_type ix=s_start; ix<=s_end ; ix++)
	{
		if ( A[ix]==(para->N_seq) )
		{
			N_seq_length++;
		}
	}
	if ( (N_seq_length*1.0/cut_length)>=(para->offN))
	{
		return 1 ;
	}
	else
	{
		Seq = A.substr(s_start,cut_length); 
		Quli =B.substr(s_start,cut_length);

	}
	return  0 ;         
}


int StortD (string  Seq1 ,  string Seq2 , DA  BB[] ,  Para_A24 * para ,int  UC[]  )
{    
		for(string::size_type ix=0; ix<(para->maxLeng) ; ix++)
		{
			int ey=UC[Seq1[ix]];
			(BB[ix]).ADD(ey);
		
			ey=UC[Seq2[ix]] ;
			(BB[ix+(para->maxLeng)]).ADD(ey);

		}
		
		return 1 ;
}


int StortD (string  Seq1 ,  DA  BB[] ,  Para_A24 * para ,int   UC[]  )
{    
		for(string::size_type ix=0; ix<(para->maxLeng) ; ix++)
		{
			//int ey=((int(Seq1[ix]))%10)%5;
			int ey=UC[Seq1[ix]];
			(BB[ix]).ADD(ey);
		}
		return 1 ;
}




////////////////////////////




int cmp(const pair<string,int>& x,const pair<string,int>& y)   
{   
	return x.second<y.second;   
}   

int cmp_first(const pair<string,int>& x,const pair<string,int>& y)   
{   
	return x.first<y.first;   
}   



void sortMapByValue(map<string,int> & tMap, vector <pair <string,int> > & tVector)   
{   
	for(map<string,int>::iterator curr=tMap.begin();curr!=tMap.end();curr++)   
	{   
		tVector.push_back(make_pair(curr->first,curr->second));   
	}   
	sort(tVector.begin(),tVector.end(),cmp);   
}  



int whether_index( const int & cut_ind ,  const string &  Fq_input_1 , map <string,int> & MyMap )
{
	if (cut_ind==0 )
	{
		return 0 ;
	}
	else if (cut_ind<0.01)
	{
		return 0 ;
	}

	string ID ,index_sub ,temp ;
	map <string,int> All_index ;
	int All_read_SE =1 ;
	igzstream Temp (Fq_input_1.c_str(),ifstream::in); 
	if(!Temp.good())
	{
		cerr<<"can't open the file " <<Fq_input_1<<endl ;
		return  0 ;
	}

	// @BB00FCABXX:6:1:1528:2000#0/1
	// @A812M3ABXX:1:1101:1307:2206#CAGATCAT/1 
	getline(Temp,ID);
	string::size_type pos_one=ID.rfind('#')+1;
	string::size_type pos_two=ID.rfind('/');
	string::size_type leng_cut=pos_two-pos_one;
	if ( leng_cut< 3 )
	{
		Temp.close();
		Temp.clear();
		cerr<<"The file " <<Fq_input_1<<" is not index file "<<endl ;
		cerr<<"so -indexcut maybe no use for this"<<endl ;
		return 0 ;
	}
	index_sub=ID.substr(pos_one,leng_cut);
	All_index[index_sub]++;    getline(Temp,temp);  getline(Temp,temp);  getline(Temp,temp);
	while(!Temp.eof() )
	{
		getline(Temp,ID);
		pos_one=ID.rfind('#')+1;
		pos_two=ID.rfind('/');
		leng_cut=pos_two-pos_one;
		index_sub=ID.substr(pos_one,leng_cut);
		All_index[index_sub]++;
		All_read_SE++;
		getline(Temp,temp);  getline(Temp,temp);  getline(Temp,temp);
	}
	Temp.close();
	Temp.clear();
	vector <pair <string,int> > tVector;   
	sortMapByValue(All_index,tVector);

	int Temp_cout_off = All_read_SE ;
	map  <string,int> ::const_iterator map_it=All_index.begin();

	for(int i=tVector.size()-1;  i>=0 ; i--)
	{   
		if ((Temp_cout_off*1.0/All_read_SE)>0.9)
		{
			MyMap[tVector[i].first]=tVector[i].second ;
			cout<<" Good  Index "<<tVector[i].first<<" : "<<tVector[i].second<<endl;
		}
		else
		{
			if (tVector[i].first =="") { Temp_cout_off=Temp_cout_off-tVector[i].second ; ;continue ;}
			cout<<"Remove Index "<<tVector[i].first<<" : "<<tVector[i].second<<endl;
		}        
		Temp_cout_off=Temp_cout_off-tVector[i].second ;

	}   

	return 1 ;
}

int Adapter_Map(  const string & Adapter_input ,  vector <pair <string,int> > & tVector  )
{
	igzstream Temp (Adapter_input.c_str(),ifstream::in); 
	if(!Temp.good())
	{
		cerr<<"can't open the file " <<Adapter_input<<endl ;
		return  0 ;
	}
	string line  ,ID_seq;
	int temp ,start ;
	getline(Temp,line);
	vector <pair <string,int> > Vector_temp ;
	while(!Temp.eof())
	{
		getline(Temp,line);
		istringstream isone (line,istringstream::in);
		isone>>ID_seq>>temp>>start ;
		ID_seq="@"+ID_seq;
		Vector_temp.push_back(make_pair(ID_seq,start));
	}
	Temp.close();
	Temp.clear();

	sort(Vector_temp.begin(),Vector_temp.end(),cmp_first);
	ID_seq=(Vector_temp[0]).first ;
	start=(Vector_temp[0]).second ;    
	int Vector_temp_size=Vector_temp.size();
	for( int i = 1 ; i < Vector_temp_size; i++)
	{
		line=(Vector_temp[i]).first ;
		temp=(Vector_temp[i]).second ;
		if (line!=ID_seq)
		{
			tVector.push_back(make_pair(ID_seq,start));
			ID_seq=line ; start=temp ;
		}
		else if (start > temp)
		{
			start=temp ;
		}
	}
	tVector.push_back(make_pair(ID_seq,start));
	return 1 ;
}


int BinSearch( const vector <pair <string,int> >  tVector  , string  & ID , int low ,int & hight )
{
	//   int low=0 , hight=tVector.size()-1 ,;       
	int   mid ;
	while(low<hight)
	{
		mid=int((low+hight)/2);
		if( (tVector[low]).first == ID )
		{         
			return  low ;
		}
		else if (   (tVector[low]).first > ID )
		{
			low=mid+1 ;
		}
		else
		{
			hight=mid-1 ;
		}
	}
	return  (-1) ;
}

//int Remove_Adapter(  string & ID , string & Seq ,  string & Quli , map <string,int>  AdapterMapA   )
//int Remove_Adapter(  string & ID , string & Seq ,  string & Quli ,std::hash_map <string,int>  AdapterMapA   )    
int Remove_Adapter  ( string & ID , string & Seq ,  string & Quli ,  const vector <pair <string,int> >  tVector     ,int  &  end_of_Search , Para_A24 * para  )
	//int Remove_Adapter(  string & ID , string & Seq ,  string & Quli ,hash_map <string,int>  AdapterMapA   )    
{
	int key=BinSearch (tVector , ID , 0 ,  end_of_Search ) ;
	int value ;
	if(key==-1)
	{
		return 1 ;
	}
	else
	{
		value=(tVector[key]).second ;
		end_of_Search=key ; 
	}

	if ( value < (para->minLeng) )
	{
		return 0 ;
	}
	else 
	{
		Seq=Seq.substr(0,value);
		Quli=Quli.substr(0,value);
		return 1 ;
	}
}

int FQSE_Filter( Para_A24 * para )
{
	int if_run_index=0 ;
	int if_runAdapter_one = 0 ;

	map <string,int> Good_index  ;
	vector <pair <string,int> > Adapter_one ;

	if ((para->indexcut)==1)
	{
		if_run_index=whether_index(para->indexcut,para->InFq1,Good_index);
	}

	if ( (para->adapter1)!="NA" )
	{
		if_runAdapter_one= Adapter_Map(para->adapter1 ,Adapter_one) ;
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

	if (if_run_index !=0  && if_runAdapter_one==0 )
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
			if (trim (seq_1,Quly_1,  para )==1)
			{
				continue ;
			}
			int seq_length=seq_1.size();
			read_Second+=1 ; Base_Second+=seq_length ;
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
			if  (seq_length==(para->maxLeng))
			{
			StortD ( seq_1 , BaseD ,para ,UC );
			}
		}
	}
	else if (if_run_index ==0 && if_runAdapter_one==0  )
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

			if (trim (seq_1,Quly_1, para  )==1)
			{
				continue ;
			}
			
			int seq_length=seq_1.size();
			read_Second+=1 ; Base_Second+=seq_length ;
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
			if  (seq_length==(para->maxLeng))
			{
			StortD ( seq_1 , BaseD ,para ,UC );
			}


		}
	}
	else if (if_run_index !=0 && if_runAdapter_one!=0 )
	{

		string::size_type pos_two,pos_one ,leng_cut ; 
		string index_sub ;
		int end_of_Adapter_one = Adapter_one.size()-1 ;
		while(!IN_1.eof())
		{
			getline(IN_1,ID_1);
			getline(IN_1,seq_1);
			getline(IN_1,temp_1);
			getline(IN_1,Quly_1);
			if (Remove_Adapter(ID_1,seq_1,Quly_1,Adapter_one,end_of_Adapter_one, para )==0)
			{
				continue ;           
			}

			pos_one=ID_1.rfind('#')+1;
			pos_two=ID_1.rfind('/');
			leng_cut=pos_two-pos_one;
			index_sub=ID_1.substr(pos_one,leng_cut);
			if (! Good_index.count(index_sub))
			{
				continue ;
			}
			read_Fisrt+=1 ; Base_Fisrt+=seq_1.size(); 
			if (trim (seq_1,Quly_1, para  )==1)
			{
				continue ;
			}
			read_Second+=1 ; Base_Second+=seq_1.size();
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
		}
	}
	else if (if_run_index ==0 && if_runAdapter_one!=0  )
	{
		int end_of_Adapter_one = Adapter_one.size()-1 ;
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
			if (Remove_Adapter(ID_1,seq_1,Quly_1,Adapter_one,end_of_Adapter_one,para)==0)
			{
				continue ;           
			}
			if (trim (seq_1,Quly_1, para  )==1)
			{
				continue ;
			}
			read_Second+=1 ; Base_Second+=seq_1.size();
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
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
int FQ_Filter_main(int argc, char **argv)
{
	Para_A24 * para = new Para_A24;        
	int Flag_para=parse_cmd_A24(argc, argv ,para ) ;
	if ( Flag_para ==0)
	{
		delete  para ; 
		return 1;
	}
	int tmpShiftQ=GetShiftQ(para->InFq1);
	if (tmpShiftQ==33)
	{
		para->minBaseQ='#';
	}

	else if ( Flag_para ==1)
	{
		FQSE_Filter (para );
		delete para;
		return  2;
	}
	int if_run_index=0 ;
	int if_runAdapter_one = 0 ;
	int if_runAdapter_two = 0 ;
	map <string,int> Good_index  ;

	//    map <string,int> Adapter_one ;
	//    map <string,int> Adapter_two ;
	vector <pair <string,int> > Adapter_one ;
	vector <pair <string,int> > Adapter_two ;

	//  hash_map <const char *,int> Adapter_one ;
	//  hash_map <const char *,int> Adapter_two ;

	if ((para->indexcut)==1)
	{
		if_run_index=whether_index(para->indexcut,para->InFq1,Good_index);
	}

	if ( (para->adapter1)!="NA" )
	{
		if_runAdapter_one= Adapter_Map(para->adapter1 ,Adapter_one) ;
	}

	if ( (para->adapter2)!="NA" )
	{
		if_runAdapter_two= Adapter_Map(para->adapter1 ,Adapter_two) ;
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


	int UC[256]={4};
    UC['A']=0; UC['C']=1;UC['G']=3;UC['T']=2 ; UC['N']=4;
    UC['a']=0; UC['c']=1;UC['g']=3;UC['t']=2 ; UC['n']=4;




	string ID_1 ,seq_1,temp_1,Quly_1 ;
	string ID_2 ,seq_2,temp_2,Quly_2 ;
	long long read_Fisrt=0 , read_Second=0 ,Base_Fisrt=0 ,Base_Second =0;

	if (if_run_index !=0  && if_runAdapter_one==0 && if_runAdapter_two==0 )
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
			if (trim (seq_1,Quly_1,  para   )==1)
			{
				continue ;
			}
			if (trim (seq_2,Quly_2, para   )==1)
			{
				continue ;            
			}
			
			 int Read_1_length=seq_1.size();
			 int Read_2_length=seq_2.size();

			read_Second+=1 ; Base_Second+=Read_1_length; Base_Second+=Read_2_length;
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
			OUT_2<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n" ;
			if (   (Read_1_length  == (para->maxLeng))  && ( Read_2_length  == (para->maxLeng) )  )
			{
			StortD ( seq_1 ,  seq_2 , BaseD ,para ,UC );
			}
		}
	}
	else if (if_run_index ==0 && if_runAdapter_one==0 && if_runAdapter_two==0 )
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

			if (trim (seq_1,Quly_1, para  )==1)
			{
				continue ;
			}
			if (trim (seq_2,Quly_2, para  )==1)
			{
				continue ;            
			}
		

			int Read_1_length=seq_1.size();
			int Read_2_length=seq_2.size();

			read_Second+=1 ; Base_Second+=Read_1_length; Base_Second+=Read_2_length;
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
			OUT_2<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n" ;
			if (   (Read_1_length  == (para->maxLeng))  && ( Read_2_length  == (para->maxLeng) )  )
			{
			StortD ( seq_1 ,  seq_2 , BaseD ,para ,UC );
			}

		}
	}
	else if (if_run_index !=0 && if_runAdapter_one!=0 && if_runAdapter_two!=0 )
	{

		string::size_type pos_two,pos_one ,leng_cut ; 
		string index_sub ;
		int end_of_Adapter_one = Adapter_one.size()-1 ;
		int end_of_Adapter_two = Adapter_two.size()-1 ;
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
			//int Remove_Adapter(  string & ID , string & Seq ,  string & Quli ,const  map <string,int> & AdapterMap  )
			if (Remove_Adapter(ID_1,seq_1,Quly_1,Adapter_one,end_of_Adapter_one, para )==0)
			{
				continue ;           
			}
			if (Remove_Adapter(ID_2,seq_2,Quly_2,Adapter_two,end_of_Adapter_two, para )==0)
			{
				continue ;           
			}

			pos_one=ID_1.rfind('#')+1;
			pos_two=ID_1.rfind('/');
			leng_cut=pos_two-pos_one;
			index_sub=ID_1.substr(pos_one,leng_cut);
			if (! Good_index.count(index_sub))
			{
				continue ;
			}
			read_Fisrt+=1 ; Base_Fisrt+=seq_1.size();  Base_Fisrt+=seq_2.size();
			if (trim (seq_1,Quly_1, para  )==1)
			{
				continue ;
			}
			if (trim (seq_2,Quly_2, para  )==1)
			{
				continue ;            
			}
			read_Second+=1 ; Base_Second+=seq_1.size(); Base_Second+=seq_2.size();
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
			OUT_2<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n" ;
		}
	}


	else if (if_run_index ==0 && if_runAdapter_one!=0 && if_runAdapter_two!=0 )
	{
		int end_of_Adapter_one = Adapter_one.size()-1 ;
		int end_of_Adapter_two = Adapter_two.size()-1 ;
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
			if (Remove_Adapter(ID_1,seq_1,Quly_1,Adapter_one,end_of_Adapter_one,para)==0)
			{
				continue ;           
			}
			if (Remove_Adapter(ID_2,seq_2,Quly_2,Adapter_two,end_of_Adapter_two,para)==0)
			{
				continue ;           
			}
			if (trim (seq_1,Quly_1, para  )==1)
			{
				continue ;
			}
			if (trim (seq_2,Quly_2, para  )==1)
			{
				continue ;            
			}
			read_Second+=1 ; Base_Second+=seq_1.size(); Base_Second+=seq_2.size();
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
			OUT_2<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n" ;
		}
	}

	else if (if_run_index ==0 && if_runAdapter_one ==0 && if_runAdapter_two!=0 )
	{
		int end_of_Adapter_two = Adapter_two.size()-1 ;

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
			if (ID_1.empty())   {   continue ;  }

			read_Fisrt+=1 ; Base_Fisrt+=seq_1.size();  Base_Fisrt+=seq_2.size();
			if (Remove_Adapter(ID_2,seq_2,Quly_2,Adapter_two,end_of_Adapter_two,para)==0)
			{
				continue ;           
			}
			if (trim (seq_1,Quly_1, para  )==1)
			{
				continue ;
			}
			if (trim (seq_2,Quly_2, para  )==1)
			{
				continue ;            
			}
			read_Second+=1 ; Base_Second+=seq_1.size(); Base_Second+=seq_2.size();
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
			OUT_2<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n" ;
		}
	}

	else if (if_run_index ==0 && if_runAdapter_one!=0 && if_runAdapter_two==0 )
	{

		int end_of_Adapter_one = Adapter_one.size()-1 ;
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
			if (ID_1.empty())   {   continue ;  }
			read_Fisrt+=1 ; Base_Fisrt+=seq_1.size();  Base_Fisrt+=seq_2.size();
			if (Remove_Adapter(ID_1,seq_1,Quly_1,Adapter_one,end_of_Adapter_one,para)==0)
			{
				continue ;           
			}
			if (trim (seq_1,Quly_1, para  )==1)
			{
				continue ;
			}
			if (trim (seq_2,Quly_2, para  )==1)
			{
				continue ;            
			}
			read_Second+=1 ; Base_Second+=seq_1.size(); Base_Second+=seq_2.size();
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
			OUT_2<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n" ;
		}
	}


	else if (if_run_index !=0 && if_runAdapter_one==0 && if_runAdapter_two!=0 )
	{

		string::size_type pos_two,pos_one ,leng_cut ; 
		string index_sub ;
		int end_of_Adapter_two = Adapter_two.size()-1 ;

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
			if (Remove_Adapter(ID_2,seq_2,Quly_2,Adapter_two,end_of_Adapter_two,para)==0)
			{
				continue ;           
			}

			pos_one=ID_1.rfind('#')+1;
			pos_two=ID_1.rfind('/');
			leng_cut=pos_two-pos_one;
			index_sub=ID_1.substr(pos_one,leng_cut);
			if (! Good_index.count(index_sub))
			{
				continue ;
			}
			read_Fisrt+=1 ; Base_Fisrt+=seq_1.size();  Base_Fisrt+=seq_2.size();
			if (trim (seq_1,Quly_1, para  )==1)
			{
				continue ;
			}
			if (trim (seq_2,Quly_2, para  )==1)
			{
				continue ;            
			}
			read_Second+=1 ; Base_Second+=seq_1.size(); Base_Second+=seq_2.size();
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
			OUT_2<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n" ;
		}
	}


	else if (if_run_index !=0 && if_runAdapter_one!=0 && if_runAdapter_two==0 )
	{

		string::size_type pos_two,pos_one ,leng_cut ; 
		string index_sub ;
		int end_of_Adapter_one = Adapter_one.size()-1 ;

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
			if (Remove_Adapter(ID_1,seq_1,Quly_1,Adapter_one,end_of_Adapter_one,para)==0)
			{
				continue ;           
			}            
			pos_one=ID_1.rfind('#')+1;
			pos_two=ID_1.rfind('/');
			leng_cut=pos_two-pos_one;
			index_sub=ID_1.substr(pos_one,leng_cut);
			if (! Good_index.count(index_sub))
			{
				continue ;
			}
			read_Fisrt+=1 ; Base_Fisrt+=seq_1.size();  Base_Fisrt+=seq_2.size();
			if (trim (seq_1,Quly_1, para  )==1)
			{
				continue ;
			}
			if (trim (seq_2,Quly_2, para  )==1)
			{
				continue ;            
			}
			read_Second+=1 ; Base_Second+=seq_1.size(); Base_Second+=seq_2.size();
			OUT_1<<ID_1<<"\n"<<seq_1<<"\n"<<temp_1<<"\n"<<Quly_1<<"\n" ;
			OUT_2<<ID_2<<"\n"<<seq_2<<"\n"<<temp_2<<"\n"<<Quly_2<<"\n" ;
		}
	}

	//    Base_Fisrt=Base_Fisrt-seq_1.size(); 
	//    Base_Fisrt=Base_Fisrt-seq_2.size(); 

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
		cout<<ii+1<<"\t"<<(BaseD[ii]).get(0)<<"\t"<<(BaseD[ii]).get(1)<<"\t"<<(BaseD[ii]).get(2)<<"\t"<<(BaseD[ii]).get(3)<<"\t"<<(BaseD[ii]).get(4)<<"\t"<<endl;
	}

	delete para ;
	delete [] BaseD ;
	return 0 ;
}

#endif  // FQ_Filter_H_

///////// swimming in the sky and flying in the sea ////////////
