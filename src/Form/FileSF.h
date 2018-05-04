#ifndef FileSF_H_
#define FileSF_H_
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

int  print_File01()
{
	cout <<""
		"\n"
		"\tUsage: SF  -InFile1 <A.in> -InFile2 <B.in> \n"
		"\n"
		"\t\t-a    <str>    Input File1 for same or diff\n"
		"\t\t-b    <str>    Input File2 for same or diff\n"
		"\t\t\n"
		"\t\t-ID1  <int>    The ID in the File1 row 1,2 [2]\n"
		"\t\t-ID2  <int>    The ID in the File2 row 1,2 [2]\n"
		"\t\t-o    <str>    OutPut File, or STDOUT\n"
		"\t\t-s    <int>    1: same in file2 ;2 : diff in file1 ; [1]\n"
		"\t\t               3: diff in file2 ;4 : same in file1 ; [1]\n"
		"\t\t               6: same in file1 & file2  7,8:same repeat \n"
		"\t\t-h             show this help\n"
		"\n";
	return 1;
}

int parse_File01(int argc, char **argv , In3str1v * para_File01 )
{
	if (argc <=2 ) {print_File01();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InFile1"  || flag=="a")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_File01->InStr1=argv[i];
		}
		else if (flag  ==  "InFile2"  || flag=="b")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_File01->InStr2=argv[i];
		}
		else if (flag  ==  "OutPut"  || flag=="o")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_File01->InStr3=argv[i];
		}
		else if (flag  == "Swich" || flag=="s")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			(para_File01->InInt)=atoi(argv[i]);
		}
		else if (flag  == "ID1")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			(para_File01->List)[0]=argv[i];
		}
		else if (flag  == "ID2" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			(para_File01->List)[1]=argv[i];
		}
		else if (flag  == "help" || flag=="h")
		{
			print_File01();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	if  ((para_File01->InStr1).empty() || ((para_File01->InStr2).empty()))
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}
	if (!(para_File01->InStr3).empty())
	{
		(para_File01->InStr3)=add_Asuffix(para_File01->InStr3) ;
	}
	return 1 ;
}
////////////////////////////

int File_SF_main(int argc,char *argv[])
	//int main(int argc,char *argv[])
{    
	In3str1v * para_File01 = new In3str1v ;
	(para_File01->InInt)=1;
	(para_File01->List).push_back("2");
	(para_File01->List).push_back("2");
	if (parse_File01( argc, argv, para_File01)==0)
	{
		delete para_File01 ;
		return 1;
	}
	igzstream IN_F1 ((para_File01->InStr1).c_str(),ifstream::in);
	igzstream IN_F2 ((para_File01->InStr2).c_str(),ifstream::in);
	if (IN_F1.fail())
	{
		cerr<<"Can't open file "<<(para_File01->InStr1)<<endl;
		return 1;
	}
	if (IN_F2.fail())
	{
		cerr<<"Can't open file "<<(para_File01->InStr2)<<endl;
		return 1;
	}

	vector <string> inf;
	vector <int> ID1;
	split((para_File01->List)[0],inf,",");
	for (unsigned int ii=0 ; ii<inf.size(); ii++)
	{
		ID1.push_back(atoi(inf[ii].c_str())-1);
	}
	int lengthID1=ID1.size();
	inf.clear();
	vector <int> ID2;
	split((para_File01->List)[1],inf,",");
	for (unsigned int ii=0 ; ii<inf.size(); ii++)
	{
		ID2.push_back(atoi(inf[ii].c_str())-1);
	}

	int lengthID2=ID2.size();

	if ((para_File01->InStr3).empty())
	{
		map <string,bool> Hash ;
		if ((para_File01->InInt)==1)
		{
			string line ;
			while(getline(IN_F1,line))
			{
				vector <string> Tmp ;
				split(line,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID1 ; ii++)
				{
					key+=Tmp[ID1[ii]];
				}
				Hash.insert(map <string, bool> :: value_type(key,true));
			}
			while(!IN_F2.eof())
			{
				string  line2 ;
				getline(IN_F2,line2) ;
				if (line2.length()<1){continue  ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				if (Hash.find(key)!=Hash.end())
				{
					cout<<line2<<"\n";
				}
			}
		}
		else if ((para_File01->InInt)==3)
		{
			string line ;
			while(getline(IN_F1,line))
			{
				vector <string> Tmp ;
				split(line,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID1 ; ii++)
				{
					key+=Tmp[ID1[ii]];
				}
				Hash.insert(map <string, bool> :: value_type(key,true));
			}
			while(!IN_F2.eof())
			{
				string  line2 ;
				getline(IN_F2,line2) ;
				if (line2.length()<1){continue  ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				if (Hash.find(key)==Hash.end())
				{
					cout<<line2<<"\n";
				}
			}
		}
		else if ((para_File01->InInt)==4)
		{
			string line ;
			while(getline(IN_F2,line))
			{
				vector <string> Tmp ;
				split(line,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				Hash.insert(map <string, bool> :: value_type(key,true));
			}
			while(!IN_F1.eof())
			{
				string  line2 ;
				getline(IN_F1,line2) ;
				if (line2.length()<1){continue  ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID1 ; ii++)
				{
					key+=Tmp[ID1[ii]];
				}
				if (Hash.find(key)!=Hash.end())
				{
					cout<<line2<<"\n";
				}
			}
		}
		else if ((para_File01->InInt)==2)
		{
			string line ;
			while(getline(IN_F2,line))
			{
				vector <string> Tmp ;
				split(line,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				Hash.insert(map <string, bool> :: value_type(key,true));
			}
			while(!IN_F1.eof())
			{
				string  line2 ;
				getline(IN_F1,line2) ;
				if (line2.length()<1){continue  ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID1 ; ii++)
				{
					key+=Tmp[ID1[ii]];
				}
				if (Hash.find(key)==Hash.end())
				{
					cout<<line2<<"\n";
				}
			}
		}
		else if ((para_File01->InInt)==6)
		{
			string line ;
			map <string,string> Hashstr ;
			while(getline(IN_F2,line))
			{
				vector <string> Tmp ;
				split(line,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				Hashstr.insert(map <string, string> :: value_type(key,line));
			}
			while(!IN_F1.eof())
			{
				string  line2 ;
				getline(IN_F1,line2) ;
				if (line2.length()<1){continue  ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID1 ; ii++)
				{
					key+=Tmp[ID1[ii]];
				}
				map <string,string> ::  iterator it=Hashstr.find(key);
				if (it!=Hashstr.end())
				{                    
					cout<<line2<<"\n";
					cout<<it->second<<"\n";
				}
			}
		}
		else if ((para_File01->InInt)==5)
		{
			string line ;
			map <string,string> Hashstr ;
			while(getline(IN_F2,line))
			{
				vector <string> Tmp ;
				split(line,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				Hashstr.insert(map <string, string> :: value_type(key,line));
			}
			cout<<"<<--Diff in File1-->>"<<"\n";
			while(!IN_F1.eof())
			{
				string  line2 ;
				getline(IN_F1,line2) ;
				if (line2.length()<1){continue ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID1 ; ii++)
				{
					key+=Tmp[ID1[ii]];
				}
				map <string,string> ::  iterator it=Hashstr.find(key);
				if (it==Hashstr.end())
				{                    
					cout<<line2<<"\n";
				}
				else
				{
					Hashstr.erase(it);
				}
			}
			cout<<"<<--Diff in File2-->>"<<"\n";
			map <string,string> ::  iterator it=Hashstr.begin();
			for( ; it!=Hashstr.end() ; it++)
			{
				cout<<it->second<<"\n";
			}
		}
		else if ((para_File01->InInt)==8)
		{
			string line ;
			while(getline(IN_F2,line))
			{
				vector <string> Tmp ;
				split(line,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				Hash.insert(map <string, bool> :: value_type(key,false));
			}
			cout<<"<<--Same in File1-->>"<<"\n";

			while(!IN_F1.eof())
			{
				string  line2 ;
				getline(IN_F1,line2) ;
				if (line2.length()<1){continue  ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				map <string,bool> ::  iterator it=Hash.find(key);
				if (it!=Hash.end())
				{
					cout<<line2<<"\n";
					it->second=true;
				}
			}

			cout<<"<<--Same in File2-->>"<<"\n";

			igzstream IN_F2_S ((para_File01->InStr2).c_str(),ifstream::in);
			if (IN_F2_S.fail())
			{
				cerr<<"Can't open file "<<(para_File01->InStr2)<<endl;
				return 1;
			}

			while(!IN_F2_S.eof())
			{
				string  line2 ;
				getline(IN_F2_S,line2) ;
				if (line2.length()<1){continue  ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				map <string,bool> ::  iterator it=Hash.find(key);
				if ((it->second))
				{
					cout<<line2<<"\n";
				}
			}
		}
		else if ((para_File01->InInt)==7)
		{
			map <string,vector <string> > Hashstr_2 ;
			string  line2 ;
			while(getline(IN_F2,line2))
			{
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				map <string,vector <string> > ::iterator it=Hashstr_2.find(key);
				if (it==Hashstr_2.end())
				{
					vector  <string> tmpA (1,line2);

					Hashstr_2.insert(map  <string,vector <string> > :: value_type(key,tmpA));
				}
				else
				{
					(it->second).push_back(line2);
				}
			}

			map <string,vector <string> > Hashstr_1 ;
			string  line1 ;
			while(getline(IN_F1,line1))
			{
				vector <string> Tmp ;
				split(line1,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID1 ; ii++)
				{
					key+=Tmp[ID1[ii]];
				}
				map <string,vector <string> > ::iterator it=Hashstr_1.find(key);
				if (it==Hashstr_1.end())
				{
					vector  <string> tmpB(1,line1);
					Hashstr_1.insert(map  <string,vector <string> > :: value_type(key,tmpB));
				}
				else
				{
					(it->second).push_back(line1);
				}
			}


			map <string,vector <string> > ::const_iterator map_it=Hashstr_1.begin();
			while(map_it!=Hashstr_1.end())
			{
				string key=map_it->first;
				map <string,vector <string> > ::iterator it=Hashstr_2.find(key);
				if (it!=Hashstr_2.end())
				{
					for (int A=0;A<(map_it->second).size(); A++)
					{
						cout<<(map_it->second)[A]<<"\n";
					}
					for (int A=0;A<(it->second).size(); A++)
					{
						cout<<(it->second)[A]<<"\n";
					}
				}
				map_it++;
			}
		}

	}
	else
	{
		ogzstream OUT ((para_File01->InStr3).c_str());
		if (OUT.fail())
		{
			cerr<<"can't open OutPut File"<<(para_File01->InStr3)<<endl;
			return 1;
		}

		map <string,bool> Hash ;
		if ((para_File01->InInt)==1)
		{
			string line ;
			while(getline(IN_F1,line))
			{
				vector <string> Tmp ;
				split(line,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID1 ; ii++)
				{
					key+=Tmp[ID1[ii]];
				}
				Hash.insert(map <string, bool> :: value_type(key,true));
			}
			while(!IN_F2.eof())
			{
				string  line2 ;
				getline(IN_F2,line2) ;
				if (line2.length()<1){continue  ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				if (Hash.find(key)!=Hash.end())
				{
					OUT<<line2<<"\n";
				}
			}
		}
		else if ((para_File01->InInt)==3)
		{
			string line ;
			while(getline(IN_F1,line))
			{
				vector <string> Tmp ;
				split(line,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID1 ; ii++)
				{
					key+=Tmp[ID1[ii]];
				}
				Hash.insert(map <string, bool> :: value_type(key,true));
			}
			while(!IN_F2.eof())
			{
				string  line2 ;
				getline(IN_F2,line2) ;
				if (line2.length()<1){continue  ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				if (Hash.find(key)==Hash.end())
				{
					OUT<<line2<<"\n";
				}
			}
		}
		else if ((para_File01->InInt)==4)
		{
			string line ;
			while(getline(IN_F2,line))
			{
				vector <string> Tmp ;
				split(line,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				Hash.insert(map <string, bool> :: value_type(key,true));
			}
			while(!IN_F1.eof())
			{
				string  line2 ;
				getline(IN_F1,line2) ;
				if (line2.length()<1){continue  ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID1 ; ii++)
				{
					key+=Tmp[ID1[ii]];
				}
				if (Hash.find(key)!=Hash.end())
				{
					OUT<<line2<<"\n";
				}
			}
		}
		else if ((para_File01->InInt)==2)
		{
			string line ;
			while(getline(IN_F2,line))
			{
				vector <string> Tmp ;
				split(line,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				Hash.insert(map <string, bool> :: value_type(key,true));
			}
			while(!IN_F1.eof())
			{
				string  line2 ;
				getline(IN_F1,line2) ;
				if (line2.length()<1){continue  ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID1 ; ii++)
				{
					key+=Tmp[ID1[ii]];
				}
				if (Hash.find(key)==Hash.end())
				{
					OUT<<line2<<"\n";
				}
			}
		}
		else if ((para_File01->InInt)==6)
		{
			string line ;
			map <string,string> Hashstr ;
			while(getline(IN_F2,line))
			{
				vector <string> Tmp ;
				split(line,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				Hashstr.insert(map <string, string> :: value_type(key,line));
			}
			while(!IN_F1.eof())
			{
				string  line2 ;
				getline(IN_F1,line2) ;
				if (line2.length()<1){continue  ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID1 ; ii++)
				{
					key+=Tmp[ID1[ii]];
				}
				map <string,string> ::  iterator it=Hashstr.find(key);
				if (it!=Hashstr.end())
				{                    
					OUT<<line2<<"\n";
					OUT<<it->second<<"\n";
				}
			}
		}
		else if ((para_File01->InInt)==5)
		{
			string line ;
			map <string,string> Hashstr ;
			while(getline(IN_F2,line))
			{
				vector <string> Tmp ;
				split(line,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				Hashstr.insert(map <string, string> :: value_type(key,line));
			}
			OUT<<"<<--Diff in File1-->>"<<"\n";
			while(!IN_F1.eof())
			{
				string  line2 ;
				getline(IN_F1,line2) ;
				if (line2.length()<1){continue ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID1 ; ii++)
				{
					key+=Tmp[ID1[ii]];
				}
				map <string,string> ::  iterator it=Hashstr.find(key);
				if (it==Hashstr.end())
				{                    
					OUT<<line2<<"\n";
				}
				else
				{
					Hashstr.erase(it);
				}
			}
			OUT<<"<<--Diff in File2-->>"<<"\n";
			map <string,string> ::  iterator it=Hashstr.begin();
			for( ; it!=Hashstr.end() ; it++)
			{
				OUT<<it->second<<"\n";
			}
		}
		else if ((para_File01->InInt)==8)
		{
			string line ;
			while(getline(IN_F2,line))
			{
				vector <string> Tmp ;
				split(line,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				Hash.insert(map <string, bool> :: value_type(key,false));
			}
			OUT<<"<<--Same in File1-->>"<<"\n";

			while(!IN_F1.eof())
			{
				string  line2 ;
				getline(IN_F1,line2) ;
				if (line2.length()<1){continue  ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				map <string,bool> ::  iterator it=Hash.find(key);
				if (it!=Hash.end())
				{
					OUT<<line2<<"\n";
					it->second=true;
				}
			}

			OUT<<"<<--Same in File2-->>"<<"\n";

			igzstream IN_F2_S ((para_File01->InStr2).c_str(),ifstream::in);
			if (IN_F2_S.fail())
			{
				cerr<<"Can't open file "<<(para_File01->InStr2)<<endl;
				return 1;
			}

			while(!IN_F2_S.eof())
			{
				string  line2 ;
				getline(IN_F2_S,line2) ;
				if (line2.length()<1){continue  ;}
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				map <string,bool> ::  iterator it=Hash.find(key);
				if ((it->second))
				{
					OUT<<line2<<"\n";
				}
			}
		}
		else if ((para_File01->InInt)==7)
		{
			map <string,vector <string> > Hashstr_2 ;
			string  line2 ;

			while(getline(IN_F2,line2))
			{
				vector <string> Tmp ;
				split(line2,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID2 ; ii++)
				{
					key+=Tmp[ID2[ii]];
				}
				map <string,vector <string> > ::iterator it=Hashstr_2.find(key);
				if (it==Hashstr_2.end())
				{
					vector  <string> tmpA (1,line2);
					Hashstr_2.insert(map  <string,vector <string> > :: value_type(key,tmpA));

				}
				else
				{
					(it->second).push_back(line2);
				}
			}

			map <string,vector <string> > Hashstr_1 ;
			string  line1 ;
			while(getline(IN_F1,line1))
			{
				vector <string> Tmp ;
				split(line1,Tmp," \t");
				string key;
				for (int ii=0 ; ii<lengthID1 ; ii++)
				{
					key+=Tmp[ID1[ii]];
				}
				map <string,vector <string> > ::iterator it=Hashstr_1.find(key);
				if (it==Hashstr_1.end())
				{
					vector  <string> tmpB(1,line1);
					Hashstr_1.insert(map  <string,vector <string> > :: value_type(key,tmpB));
				}
				else
				{
					(it->second).push_back(line1);
				}
			}


			map <string,vector <string> > ::const_iterator map_it=Hashstr_1.begin();
			while(map_it!=Hashstr_1.end())
			{
				string key=map_it->first;
				map <string,vector <string> > ::iterator it=Hashstr_2.find(key);
				if (it!=Hashstr_2.end())
				{
					for (int A=0;A<(map_it->second).size(); A++)
					{
						OUT<<(map_it->second)[A]<<"\n";
					}
					for (int A=0;A<(it->second).size(); A++)
					{
						OUT<<(it->second)[A]<<"\n";
					}
				}
				map_it++;
			}
		}

		OUT.close();
	}


	IN_F1.close();
	IN_F2.close();
	delete para_File01 ;
	return 0 ;
}

#endif // FileSF_H_

////////////////////////swimming in the sea & flying in the sky //////////////////


