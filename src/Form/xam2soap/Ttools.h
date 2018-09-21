/*
 ******************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Fan Zhang
 *CREATE DATE : 2010-7-14
 *CLASS NAME: 
 *FUNCTION : some useful function to the SOAPsnp
 *FILE NAME : tool.h
 *UPDATE DATE : 2011-9-23
 *UPDATE BY : Fan Zhang
 *******************************************************************************
 */


#ifndef Ttools_H_
#define Ttools_H_

#include <string>
#include <iostream>
//#include <strstream>
#include <vector>
#include <string>
#include <sstream>


using namespace std ;


#define NOUSE_ALIGNMENT "alignment is no use"
#define HIT_COUNT_X0		"X0:i:"
#define HIT_COUNT_H0		"H0:i:"
#define HIT_COUNT_XA		"XA:Z:"
#define HIT_COUNT_XTU		"XT:A:U"
#define HIT_COUNT_XTR		"XT:A:R"

/* integer to string*/
//std::string Int2Str(int num);
string Int2Str (int A ) ;
string  modif_cigarsize( string  cigar  );
// A function to spilt string s into vector vec according to char splitchar
void TStringSplit(std::string s, char splitchar, std::vector<std::string>& vec);
/* count the indel length from the cigar string. return the indel length and the position*/
int Tcount_indel_len(const std::string cigar, int &pos, std::string &seq, std::string &qual);
bool Tcount_indel_len_filter(const std::string cigar, int & pos,std::string &seq,std::string &qual);
/* format the sam text to the soap text*/
std::string Talignment_format(const std::string &sam_ali);
std::string Talignment_format_Filter(const std::string &sam_ali,int MappingQ , bool ID);
/* count soft clip number at the beginning*/
int Tcount_soft_clip(const std::string cigar, int &last_S);

/*
 ******************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Bill Tang
 *CREATE DATE : 2010-7-14
 *CLASS NAME: 
 *FUNCTION : definition of the useful functions
 *FILE NAME : tool.cpp
 *UPDATE DATE : 2011-9-23
 *UPDATE BY : Fan Zhang
 *******************************************************************************
 */

//#include "Ttools.h"
//using namespace std;
/* integer to string*/
/*
   std::string Int2Str(int num) {
   char num_str[256] = {0};
   sprintf(num_str, "%d", num);
   return num_str;
   }
   */
// A function to spilt string s into vector vec according to char splitchar
void TStringSplit(std::string s, char splitchar, std::vector<std::string>& vec) {
	// promise the vec is empty
	if(vec.size()>0) 
		vec.clear();
	// chomp the string
	while(s.size() > 0 && s[0] == splitchar)
		s = s.substr(1, s.length() - 1);
	while(s.size() > 0 && s[s.length() - 1] == splitchar)
		s = s.substr(0, s.length() - 1);
	int length = s.length();
	int start=0;
	for(int i=0;i<length;i++)
	{
		if(s[i] == splitchar)
		{
			vec.push_back(s.substr(start,i - start));
			while(s[i + 1] == splitchar)
				i ++;
			start = i+1;
		}
		else if(i == length-1)// attach last
		{
			vec.push_back(s.substr(start,i+1 - start));
		}
	}
}

string  modif_cigarsize( string  cigar  )
{
	int cigarsize=cigar.size();
	const char *tmp = cigar.c_str();
	size_t Final_End=cigarsize-1;
	if (tmp[Final_End] == 'H' )
	{
		while(1)
		{
			Final_End--;
			if ( tmp[Final_End] < 48  || tmp[Final_End] > 57 )
			{
				break ;
			}
		}
	}
	size_t Final_Start=-1;
	while(1)
	{
		Final_Start++;
		if ( tmp[Final_Start] < 48  || tmp[Final_Start] > 57 )
		{
			break;
		}
	}

	if  ( tmp[Final_Start] == 'H' )
	{
		Final_Start++;
	}
	else
	{
		Final_Start=0;
	}

	string H=cigar.substr(Final_Start, Final_End-Final_Start+1);
	return H;

}





//hewm  stat read_length
int Tcount_Read_len (const std::string cigar_v1 ) 
{
	string  cigar=modif_cigarsize(cigar_v1);
	const char *tmp = cigar.c_str();
	int begin = 0;
	int end;
	int len = 0;
	int Fianle_Read_length=0;
	int cigarsize=cigar.size();
	for (end = 0; end < cigarsize; end ++)
	{
		switch (tmp[end])
		{
			case 'M':
			case 'D':
			case 'N':
				{

					len= atoi(cigar.substr(begin, end - begin).c_str());
					begin = end + 1;
					Fianle_Read_length+=len;
					break;
				}
			case 'I':
			case 'S':
			case 'P':
			case 'H':
				{
					begin = end + 1;
					break;
				}
			default : break;
		}
	}
	
	return  Fianle_Read_length;
}





/* count the indel length from the cigar string. return the indel length and the position*/
int Tcount_indel_len (const std::string cigar_v1, int & pos,std::string &seq,std::string &qual) {
	string  cigar=modif_cigarsize(cigar_v1);
	const char *tmp = cigar.c_str();
	int begin = 0;
	int end;
	int len = 0;
	pos = 0;
	int add=0;
	std::string tmp_seq,tmp_qual;
	int cigarsize=cigar.size();
	for (end = 0; end < cigarsize; end ++) {
		switch (tmp[end]) {
			case 'M':
				{

					len= atoi(cigar.substr(begin, end - begin).c_str());
					tmp_seq+=seq.substr(pos,len);
					tmp_qual+=qual.substr(pos,len);
					pos +=len;
					begin = end + 1;
					break;
				}
			case 'I':
				{
					len = atoi(cigar.substr(begin, end - begin).c_str()); 
					pos +=len;
					begin = end + 1;
					break;
				}
			case 'D':
			case 'N':
				{
					len = atoi(cigar.substr(begin, end - begin).c_str()); 
					for(int i=0;i!=len;++i)
					{
						tmp_seq+="N";
						tmp_qual+="B";
					}
					begin = end +1;
					add+=len;
					break;
				}
			case 'S':
			case 'P':
				{
					len= atoi(cigar.substr(begin, end - begin).c_str());
					begin = end + 1;
					pos+=len;
					break;
				}
			case 'H':
				{
					len= atoi(cigar.substr(begin, end - begin).c_str());
					begin = end + 1;
					pos+=len;
					add-=len;
					break;
				}
			default : break;
		}
	}
	seq=tmp_seq;
	qual=tmp_qual;
	return add;
}

/* format the sam text to the soap text*/
std::string Talignment_format(const std::string &sam_ali) {
	if(sam_ali.empty()) {
		return NOUSE_ALIGNMENT;
	}
	std::vector<std::string> vec;
	TStringSplit(sam_ali, '\t', vec);
	if ((vec.size() < 11) || (vec[5] == "*") )
		return NOUSE_ALIGNMENT;
	std::string format;
	int pos = 0;
	int best_hit = 1;
	int flag = atoi(vec[1].c_str());
	int last_clip_num;
	int clip_num = Tcount_soft_clip(vec[5], last_clip_num);
	int vec9size=vec[9].size();
	if (clip_num + last_clip_num > vec9size || ((flag>>10) & 1))//change here to remove dup 20110716
		return NOUSE_ALIGNMENT;

	if (flag & (0x1 << 6)) {
		format = vec[0] + "/1" + "\t"; // add query name
	} else if (flag & (0x1 << 7)) {
		format = vec[0] + "/2" + "\t"; // add query name
	} else {
		format = vec[0] + "\t";  // add query name
	}
	//changed 20110921
	std::string tmp_seq=vec[9];//vec[9].substr(clip_num, vec[9].size() - last_clip_num - clip_num);
	std::string tmp_qual=vec[10];//vec[10].substr(clip_num, vec[10].size() - last_clip_num - clip_num);
	Tcount_indel_len(vec[5],pos,tmp_seq,tmp_qual);
	//int length_change=Tcount_indel_len(vec[5],pos,tmp_seq,tmp_qual);
	//tmp_seq=tmp_seq.substr(clip_num, tmp_seq.size()- last_clip_num - clip_num);
	//tmp_qual=tmp_qual.substr(clip_num, tmp_qual.size()- last_clip_num - clip_num);
	format += tmp_seq + "\t"; // add sequence
	format += tmp_qual + "\t"; // add quality
	int vecsize=vec.size();
	for (int i = 11; i < vecsize; i++) {
		if (vec[i].find(HIT_COUNT_XA) != std::string::npos)			
		{
			vector <string> tmp ;
			split( vec[i] ,  tmp , ";");
			best_hit=1+tmp.size(); 
			break;
		}
		else if (vec[i].find(HIT_COUNT_H0) != std::string::npos) {
			best_hit = atoi(vec[i].substr(5).c_str());
			break;
		} else if (vec[i].find(HIT_COUNT_X0) != std::string::npos) {
			best_hit = atoi(vec[i].substr(5).c_str());
			break;
		}
		else if (vec[i].find(HIT_COUNT_XTU)!= std::string::npos)
		{
			best_hit =1;
			break;
		}
		else if (vec[i].find(HIT_COUNT_XTR)!= std::string::npos)
		{
			best_hit =2;
			break;
		}
	}
	format += Int2Str(best_hit)  + "\t"; // add number of best hit.

	format += ((flag & (0x1 << 6)) ? "a" : "b"); // add a/b
	format += "\t";
	format += Int2Str(tmp_seq.size());//vec[9].size() - clip_num - last_clip_num+length_change); // add length
	format += "\t";
	format += ((flag & (0x1 << 4)) ? "-" : "+");  // add strand
	format += "\t";
	format += vec[2] + "\t";  // add chr
	//int coor_change=atoi(vec[3].c_str())+clip_num;//change pos by adding softclip length
	//char real_coord[1024];
	//sprintf(real_coord,"%d",coor_change);
	//string new_coord(real_coord);
	format += vec[3] + "\t";  // add location
	//format += Int2Str(Tcount_indel_len(vec[5], pos));  // add number of mismatch
	//correction for Bill Tang 20110711
	format += "0";//ignore mismatch setting constant to deal with indel in cigar
	format += "\t";
	//format += Int2Str(pos - clip_num);  // add indel position
	format += vec[5];
	vecsize=vec.size();
	for (int i=10;i!=vecsize;++i)
	{
		if(vec[i].find("MD:Z:")!=vec[i].npos)
		{
			std::string tmp=vec[i].substr(vec[i].find("MD:Z:")+5,vec[i].size()-5);
			format +="\t";
			format +=tmp;
		}
	}
	//if(length_change>10){exit(0);}
	return format;
}

/**
 * DATE: 2010-9-9
 * FUNCTION: count the soft clip number at the read's beginning.
 * PARAMETER: cigar: the cigar string. last_S: last soft clip's length.
 * RETURN:	soft clip number.
 */
int Tcount_soft_clip(const std::string cigar, int &last_S)
{
	const char *tmp = cigar.c_str();
	last_S = 0;
	int CigarSize= cigar.size() ;
	// count the last soft clip number.
	if (tmp[CigarSize - 1] == 'S')
	{
		int j = CigarSize - 2;
		while (tmp[j] >= '0' && tmp[j] <= '9')
		{
			j--;
		}
		last_S = atoi(cigar.substr(j + 1, CigarSize - j - 1).c_str());
	}

	int i = 0;
	for (; i < CigarSize ; ++i)
	{
		if (tmp[i] >= '0' && tmp[i] <= '9')
		{
			// jump to the next position.
			continue;
		}
		else if (tmp[i] == 'S')
		{
			// return the soft clip number.
			int num = atoi(cigar.substr(0, i).c_str());
			return num;
		}
		else
		{
			// no soft clip infront of the reads.
			break;
		}
	}
	return 0;
}













// hewm//

bool Tcount_indel_len_filter(const std::string cigar_v1, int & pos,std::string &seq,std::string &qual) {
	string  cigar=modif_cigarsize(cigar_v1);
	const char *tmp = cigar.c_str();
	int begin = 0;
	int end;
	int len = 0;
	pos = 0;
	int add=0;
	bool NoGapInder=false;
	std::string tmp_seq,tmp_qual;
	int cigarsize=cigar.size();
	for (end = 0; end < cigarsize; end ++) {
		switch (tmp[end]) {
			case 'M':
				{

					len= atoi(cigar.substr(begin, end - begin).c_str());
					tmp_seq+=seq.substr(pos,len);
					tmp_qual+=qual.substr(pos,len);
					pos +=len;
					begin = end + 1;
					break;
				}
			case 'I':
				{
					len = atoi(cigar.substr(begin, end - begin).c_str()); 
					pos +=len;
					begin = end + 1;
					NoGapInder=true ; //hewm
					break;
				}
			case 'N':
			case 'D':
				{
					len = atoi(cigar.substr(begin, end - begin).c_str()); 
					for(int i=0;i!=len;++i)
					{
						tmp_seq+="N";
						tmp_qual+="B";
					}
					begin = end +1;
					add+=len;
					NoGapInder=true ; //hewm
					break;
				}
			case 'S':
			case 'P':
				{
					len= atoi(cigar.substr(begin, end - begin).c_str());
					begin = end + 1;
					pos+=len;
					break;
				}
			case 'H':
				{
					len= atoi(cigar.substr(begin, end - begin).c_str());
					begin = end + 1;
					pos+=len;
					add-=len;
					break;
				}
			default : break;
		}

	}
	seq=tmp_seq;
	qual=tmp_qual;
	return NoGapInder ;
	//	return add;
}











/* format the sam text to the soap text*/
std::string Talignment_format_filter(const std::string &sam_ali,int MappingQ , bool FiGap ) {
	if(sam_ali.empty()) {
		return NOUSE_ALIGNMENT;
	}
	std::vector<std::string> vec;
	TStringSplit(sam_ali, '\t', vec);
	int NowmappingQ=atoi(vec[4].c_str());

	int vec_size=vec.size();
	if ((vec_size < 11) || (vec[5] == "*") || (vec[2] == "*") || (NowmappingQ<MappingQ))
	{
		return NOUSE_ALIGNMENT;
	}

	std::string format;
	int pos = 0;
	int best_hit = 1;
	int flag = atoi(vec[1].c_str());
	int last_clip_num;
	int clip_num = Tcount_soft_clip(vec[5], last_clip_num);
	int vec9size=vec[9].size();
	if (clip_num + last_clip_num > vec9size || ((flag>>10) & 1))//change here to remove dup 20110716
		return NOUSE_ALIGNMENT;

	if (flag & (0x1 << 6)) {
		format = vec[0] + "/1" + "\t"; // add query name
	} else if (flag & (0x1 << 7)) {
		format = vec[0] + "/2" + "\t"; // add query name
	} else {
		format = vec[0] + "\t";  // add query name
	}
	//changed 20110921
	std::string tmp_seq=vec[9];//vec[9].substr(clip_num, vec[9].size() - last_clip_num - clip_num);
	std::string tmp_qual=vec[10];//vec[10].substr(clip_num, vec[10].size() - last_clip_num - clip_num);
	bool Gap= Tcount_indel_len_filter(vec[5],pos,tmp_seq,tmp_qual);
	if ( (FiGap) &&  (Gap) )
	{
		return NOUSE_ALIGNMENT ;
	}
	//    ol Tcount_indel_len_filter(const std::string cigar, int & pos,std::string &seq,std::string &qual)

	//int length_change=Tcount_indel_len(vec[5],pos,tmp_seq,tmp_qual);
	//tmp_seq=tmp_seq.substr(clip_num, tmp_seq.size()- last_clip_num - clip_num);
	//tmp_qual=tmp_qual.substr(clip_num, tmp_qual.size()- last_clip_num - clip_num);
	format += tmp_seq + "\t"; // add sequence
	format += tmp_qual + "\t"; // add quality
	int vecsize=vec.size();
	for (int i = 11; i < vecsize; i++) {
		if (vec[i].find(HIT_COUNT_H0) != std::string::npos) {
			best_hit = atoi(vec[i].substr(5).c_str());
			break;
		} else if (vec[i].find(HIT_COUNT_X0) != std::string::npos) {
			best_hit = atoi(vec[i].substr(5).c_str());
			break;
		}
	}
	format += Int2Str(best_hit)  + "\t"; // add number of best hit.

	format += ((flag & (0x1 << 6)) ? "a" : "b"); // add a/b
	format += "\t";
	format += Int2Str(tmp_seq.size());//vec[9].size() - clip_num - last_clip_num+length_change); // add length
	format += "\t";
	format += ((flag & (0x1 << 4)) ? "-" : "+");  // add strand
	format += "\t";
	format += vec[2] + "\t";  // add chr
	//int coor_change=atoi(vec[3].c_str())+clip_num;//change pos by adding softclip length
	//char real_coord[1024];
	//sprintf(real_coord,"%d",coor_change);
	//string new_coord(real_coord);
	format += vec[3] + "\t";  // add location
	//format += Int2Str(Tcount_indel_len(vec[5], pos));  // add number of mismatch
	//correction for Bill Tang 20110711
	format += "0";//ignore mismatch setting constant to deal with indel in cigar
	format += "\t";
	//format += Int2Str(pos - clip_num);  // add indel position
	format += vec[5];
	vecsize=vec.size();
	for (int i=10;i!=vecsize;++i)
	{
		if(vec[i].find("MD:Z:")!=vec[i].npos)
		{
			std::string tmp=vec[i].substr(vec[i].find("MD:Z:")+5,vec[i].size()-5);
			format +="\t";
			format +=tmp;
		}
	}
	//if(length_change>10){exit(0);}
	return format;
}


#endif  // //


