
////////////////////////swimming in the sea & flying in the sky //////////////////


/*
 * DataClass.h
 *
 *  Created on: 2011-11-21
 *      Author: hewm@genomics.org.cn
 */

#ifndef DataClass_H_
#define DataClass_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <cstdlib>

using namespace std;

typedef char Base;
typedef unsigned long long   ubit64_t;
typedef long double LDBL ;
typedef long long  llong ;

string   Int2Base[4][4]=
{
	{"A", "M" , "W" , "R" },

	{"M", "C" , "Y" , "S" },

	{"W", "Y" , "T" , "K" },

	{"R", "S" , "K" , "G" },
} ;

Base ArryBase[5]={'A','C','T','G','N'};
string BaseNum[5]={"A","C","T","G","N"};



LDBL get_quality_error ( char Q , int minQ )
{
	int Q_int=Q-minQ ;
	//cerr<<"Qming:"<<(para->minQ)<<endl;
	LDBL error=pow(10,-0.1*Q_int);
	return error ;
} ;


double  get_het( int A_length ,int B_length , double HET )
{
	double het_pr=0.25-0.05*((A_length+B_length)/5);
	double AB_het=max(HET,het_pr);
	if  ( B_length==0 || A_length==0 )
	{
		AB_het=HET;
	} 
	return AB_het ;
}



///////////////// q_seq  for site ////////
class Site {

	string q_base[5];
	int sum_hit ;

	public:

	Site ( void )
	{
		q_base[0]=q_base[1]=q_base[2]=q_base[3]=q_base[4]="";
		sum_hit=0;
	}

	int  Add_Quality (Base A ,Base QA,int hit )
	{        
		A=toupper(A);
		if (A == 'N') {return 0 ;}
		else
		{
			int alleleA_Num=((A>>1)&3);            
			q_base[alleleA_Num]= q_base[alleleA_Num]+QA ;
			sum_hit+=hit ;
			return 1 ;
		}
	}

	void Destory_Site()
	{
		q_base[0]=q_base[1]=q_base[2]=q_base[3]=q_base[4]="";
		sum_hit=0;
		//    return 1 ;
	}

	string get_q(int bin_base) 
	{
		return q_base[bin_base];
	}
	int get_hit() 
	{
		return sum_hit ;
	}


	int  CopySame ( Site B )
	{
		q_base[0]=B.get_q(0);
		q_base[1]=B.get_q(1);
		q_base[2]=B.get_q(2);
		q_base[3]=B.get_q(3);
		sum_hit=B.get_hit();
		return 1; 
		//    return *this ;
	}

	string get_q(Base A)
	{
		int bin_base=((A>>1)&3);
		return q_base[bin_base];
	}
	int get_q_length(int A) 
	{
		int length=(q_base[A]).length();
		return length ;
	}

	int get_q_length(Base A) 
	{
		int bin_base=((A>>1)&3);
		int length=(q_base[bin_base]).length();
		return length ;
	}

	int get_Two_allele(int  Ref , Base & B , int & Ref_length , int & B_length )
	{
		//        int RefAllele=((Ref>>1)&3) ;
		Ref_length=q_base[Ref].length();
		int tmpB=4;
		B_length=0;
		for (int i=0 ; i<4 ; i++)
		{
			if (Ref^i)
			{
				int tmpL=q_base[i].length();
				if (B_length<tmpL)
				{
					tmpB=i; 
					B_length=tmpL ;
				}
			}
		}
		B=ArryBase[tmpB];
		return tmpB ;
	}

	double get_mean_hit( int & SumDepth )
	{
		if (sum_hit==0)
		{
			SumDepth=0;
			return 25.0 ;
		}
		else
		{
			SumDepth=((q_base[0]).length())+((q_base[1]).length())+((q_base[2]).length())+((q_base[3]).length()) ;
			double mean_hit=(sum_hit*1.0)/(SumDepth) ;
			return mean_hit ;
		}
	}

	string  Genotype_Likelihood ( int A , int  B ,int A_length, int B_length ,int minQ )
	{
		LDBL A_LLHD=1 ;
		LDBL B_LLHD=1 ;
		A_length = (A_length<125)?A_length:125;
		B_length = (B_length<125)?B_length:125;

		for(int ix=0 ; ix<A_length  ; ix++)
		{
			LDBL  A_ix_error=get_quality_error (q_base[A][ix] , minQ );
			A_LLHD=A_LLHD*(1-A_ix_error);
			B_LLHD=B_LLHD*A_ix_error ;
		}

		for(int ix=0 ; ix<B_length  ; ix++)
		{
			LDBL  B_ix_error=get_quality_error (q_base[B][ix],  minQ);
			B_LLHD=B_LLHD*(1-B_ix_error);
			A_LLHD=A_LLHD*B_ix_error ;
		}
		string genotype="N";

		if (A_LLHD>=B_LLHD )
		{
			genotype=BaseNum[A];
		}
		else
		{
			genotype=BaseNum[B];
		}

		return genotype  ; 
	}

	string Genotype_Likelihood(int A ,int B,int A_length, int B_length , int minQ , double HET )
	{
		LDBL AA_LLHD=1.0 ;
		LDBL AB_LLHD=1.0 ;
		LDBL BB_LLHD=1.0 ;
		A_length = (A_length<125)?A_length:125;
		B_length = (B_length<125)?B_length:125;

		for(int ix=0 ; ix<A_length ; ix++)
		{
			LDBL  A_ix_error=get_quality_error ((q_base[A])[ix] , minQ); 
			AA_LLHD=AA_LLHD*(1-A_ix_error);
			BB_LLHD=BB_LLHD*A_ix_error ;
			AB_LLHD=AB_LLHD*0.5;
		}
		for(int ix=0 ; ix<B_length ; ix++)
		{
			LDBL  B_ix_error=get_quality_error ((q_base[B])[ix] ,  minQ);
			BB_LLHD=BB_LLHD*(1-B_ix_error);
			AA_LLHD=AA_LLHD*B_ix_error ;
			AB_LLHD=AB_LLHD*0.5;
		}

		double ABHet= get_het(  A_length , B_length, HET ) ;
		double Bland=(1-ABHet)/2.0 ;
		LDBL P_AA=Bland*AA_LLHD;  LDBL P_BB=Bland*BB_LLHD;  LDBL P_AB=ABHet*AB_LLHD;
		string genotype="N";

		if (P_AA>=P_BB && P_AA>=P_AB )
		{
			genotype=Int2Base[A][A];
		}
		else if (P_BB>P_AB)
		{
			genotype=Int2Base[B][B];
		}
		else
		{
			genotype=Int2Base[A][B];
		}

		return genotype  ; 
	}


	string Genotype( int &  B_length , int MinQ )
	{

		int A=0 , A_length=q_base[0].length();
		for (int ii=1 ; ii<4 ; ii++)
		{
			int Len=q_base[ii].length() ;
			if (A_length<Len)
			{
				A=ii ; 
				A_length=Len ;
			}
		}
		if (A_length==0) {return "N"  ;}

		int B=4 ;
		B_length=0 ;
		for (int i=0 ; i<4 ; i++)
		{
			if (A^i)
			{
				int tmpL=q_base[i].length();
				if (B_length<tmpL)
				{
					B=i; 
					B_length=tmpL ;
				}
			}
		}

		LDBL AA_LLHD=1.0 ;
		LDBL AB_LLHD=1.0 ;
		LDBL BB_LLHD=1.0 ;
		A_length = (A_length<125)?A_length:125;
		B_length = (B_length<125)?B_length:125;

		for(int ix=0 ; ix<A_length ; ix++)
		{
			LDBL  A_ix_error=get_quality_error ((q_base[A])[ix] , MinQ ); 
			AA_LLHD=AA_LLHD*(1-A_ix_error);
			BB_LLHD=BB_LLHD*A_ix_error ;
			AB_LLHD=AB_LLHD*0.5;
		}
		for(int ix=0 ; ix<B_length ; ix++)
		{
			LDBL  B_ix_error=get_quality_error ((q_base[B])[ix] , MinQ );
			BB_LLHD=BB_LLHD*(1-B_ix_error);
			AA_LLHD=AA_LLHD*B_ix_error ;
			AB_LLHD=AB_LLHD*0.5;
		}

		double ABHet= get_het(  A_length , B_length, 0.001 ) ;
		double Bland=(1-ABHet)/2.0 ;
		LDBL P_AA=Bland*AA_LLHD;  LDBL P_BB=Bland*BB_LLHD;  LDBL P_AB=ABHet*AB_LLHD;
		string genotype="N";

		if (P_AA>=P_BB && P_AA>=P_AB )
		{
			genotype=Int2Base[A][A];
		}
		else if (P_BB>P_AB)
		{
			genotype=Int2Base[B][B];
		}
		else
		{
			genotype=Int2Base[A][B];
		}

		return genotype  ; 
	}




	string Genotype( int &  B_length ,int MinQ , char Abase )
	{
		int A=((Abase>>1)&3), A_length=q_base[A].length();
		int B=4; B_length=0 ;
		for (int i=0 ; i<4 ; i++)
		{
			if (A^i)
			{
				int tmpL=q_base[i].length();
				if (B_length<tmpL)
				{
					B=i; 
					B_length=tmpL ;
				}
			}
		}
		LDBL AA_LLHD=1.0 ;
		LDBL AB_LLHD=1.0 ;
		LDBL BB_LLHD=1.0 ;
		A_length = (A_length<125)?A_length:125;
		B_length = (B_length<125)?B_length:125;

		for(int ix=0 ; ix<A_length ; ix++)
		{
			LDBL  A_ix_error=get_quality_error ((q_base[A])[ix] , MinQ ); 
			AA_LLHD=AA_LLHD*(1-A_ix_error);
			BB_LLHD=BB_LLHD*A_ix_error ;
			AB_LLHD=AB_LLHD*0.5;    
		}

		for(int ix=0 ; ix<B_length ; ix++)
		{
			LDBL  B_ix_error=get_quality_error ((q_base[B])[ix] ,  MinQ );
			BB_LLHD=BB_LLHD*(1-B_ix_error);
			AA_LLHD=AA_LLHD*B_ix_error ;
			AB_LLHD=AB_LLHD*0.5;
		}

		double ABHet= get_het(  A_length , B_length, 0.001 ) ;
		double Bland=(1-ABHet)/2.0 ;
		LDBL P_AA=Bland*AA_LLHD;  LDBL P_BB=Bland*BB_LLHD;  LDBL P_AB=ABHet*AB_LLHD;
		string genotype="N";

		if (P_AA>=P_BB && P_AA>=P_AB )
		{
			genotype=Int2Base[A][A];
		}
		else if (P_BB>P_AB)
		{
			genotype=Int2Base[B][B];
		}
		else
		{
			genotype=Int2Base[A][B];
		}

		return genotype  ; 
	}

	////////swimming in the sky and flying in the sea *///////////
} ;


class In3str1v {
	public:
		string InStr1 ;
		string InStr2 ;
		string InStr3 ;
		vector <string> List ;
		bool  TF ;
		int   InInt ;
		bool TF2 ;
		In3str1v()
		{
			InStr1="";
			InStr2="";
			InStr3="";
			TF=true ;
			TF2=true ;
			InInt=0 ;
		}
};


class ParaClass {
	public:
		string InPut1;
		string InPut2;
		string OutPut1;
		string OutPut2;
		string InStr1 ;
		string InStr2 ;
		string InStr3 ;

		char  Inchar1;
		char  Inchar2;
		float Infloat1;
		float Infloat2 ;

		int InInt1 ;
		int InInt2 ;
		int InInt3 ;
		unsigned int InUnInt1 ;

		llong  Inllong1 ;

		bool  TF ;
		double Indouble;
		ParaClass()
		{
			InPut1="";
			InPut2="";
			OutPut1="";
			OutPut2="";
			InStr1="";
			InStr2="";
			InStr3="";

			Inchar1='@';
			Inchar2='N';
			Infloat1=0.5;
			Infloat2=0.0;

			InInt1=0;
			InInt2=0;
			InInt3=0;
			InUnInt1=1;

			Inllong1=0;
			Indouble=1.0;
			TF=true;
		}
};


class  Region {
	public:
		llong Start ;
		llong End ;
		Region()
		{
			Start=0;
			End=0;
		}
};


class DA
{
	public :
		llong D[5];
		DA()
		{
			D[0]=0; D[1]=0;
			D[2]=0; D[3]=0;
			D[4]=0;
		}
		int ADD (int base)
		{
			D[base]++;
			return  D[base] ;
		}
		int get(int base )
		{
			return  D[base] ;
		}
};


class DepthDis
{
	public:
		ubit64_t DA;
		ubit64_t DB;
		DepthDis()
		{
			DA=0;
			DB=0;
		}
};




typedef struct {
	double r;       // a fraction between 0 and 1
	double g;       // a fraction between 0 and 1
	double b;       // a fraction between 0 and 1
} rgb;

typedef struct {
	double h;       // angle in degrees
	double s;       // a fraction between 0 and 1
	double v;       // a fraction between 0 and 1
} hsv;





//////////////// swimming in the sky and flying in the sea ////////////////


///////////////// q_seq  for site ////////


class Para_Formt01 {
	public:
		string input_file ;
		string OutSamFile ;
		string OutBamFile ;
		int PE ;
		int sort ;
		int shiftQ ;
		string Dict ;
		Para_Formt01()
		{
			input_file="" ;
			OutSamFile="" ;
			OutBamFile="" ;
			Dict="";
			PE=0 ;
			sort=0;
			shiftQ=0;
		}
};

class SamLine
{
	public:
		string RID;
		int Flag ;
		string seq;
		string Qseq;
		string cigar;
		string chr;
		llong position ;
		int mapQ ;
		string NM_i ;
		string MD ;
		bool IF ;
		int isize ;
		llong coor ;
		string XorD ;
		void copy( SamLine *  B )
		{
			RID=B->RID;   Flag=B->Flag;    seq=B->seq ;
			Qseq=B->Qseq; cigar=B->cigar; chr=B->chr;
			position=B->position ; mapQ=B->mapQ;  IF=B->IF;
			NM_i=B->NM_i ;MD=B->MD ; isize=B->isize;
			coor=B->coor; XorD=B->XorD;
		}
		SamLine()
		{
			RID="";  seq=""; Qseq=""; cigar=""; chr ="";
			Flag=0 ; position=0; mapQ =0 ; IF=false;
			NM_i=""; MD="";
			isize=0; coor=0; XorD="*";
		}
		void rm ()
		{
			RID="";  seq=""; Qseq=""; cigar=""; chr ="";
			Flag=0 ; position=0; mapQ =0 ; IF=false;
			NM_i=""; MD="";
			isize=0; coor=0; XorD="*";
			//return 1 ;
		}
		void Print(ogzstream & OUT )
		{
			OUT<<RID<<"\t"<<Flag<<"\t"<<chr<<"\t"<<position<<"\t"<<mapQ<<"\t"<<cigar<<"\t"<<XorD<<"\t"<<coor<<"\t"<<isize<<"\t"<<seq<<"\t"<<Qseq<<"\t"<<NM_i<<"\t"<<MD<<endl;
		}
		void OUT2str(string & OutStr )
		{
			OutStr=RID+"\t"+Int2Str(Flag)+"\t"+chr+"\t"+Int2Str(position)+"\t"+Int2Str(mapQ)+"\t"+cigar+"\t"+XorD+"\t"+Int2Str(coor)+"\t"+Int2Str(isize)+"\t"+seq+"\t"+Qseq+"\t"+NM_i+"\t"+MD;
		}
};


	////////swimming in the sky and flying in the sea *///////////





#endif /* DataClass_H_ */

//////////////// swimming in the sky and flying in the sea ////////////////
