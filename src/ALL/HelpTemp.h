#ifndef HelpLess_H_
#define HelpLess_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <cstdlib>
#include <stdio.h>
#include <gzstream.h>
#include "../ALL/comm.h"
#include <ncurses.h>

#define TV_MAX_GOTO  40
#define Max_Input_Line 524286
typedef long long llong ;
using namespace std;
///////////////////

class Less
{
    public:
        igzstream  File ;
        int mrow ;
        int mcol ;
        int row_shift ;
        int col_shift ;
        WINDOW *whelp , *command;
        string PringLine[Max_Input_Line];
        int Now_Pring_Line ;
        int pos ;
        llong Before ;
        bool NoEndFile ;
        string findstr ;
        int FindFlag ;
};

int Less_init ( Less  * tv , int argc , char *argv[] ) ;

int readFile( Less  * tv  , int Max )
{
    if (tv->NoEndFile)
    {
        while( (!(tv->File).eof()) && ((tv->Now_Pring_Line)<Max))
        {
            string temp;
            getline(tv->File,temp);
            (tv->PringLine)[(tv->Now_Pring_Line)]=temp;
            (tv->Now_Pring_Line)++;
        }
        if ((tv->File).eof())
        {
            tv->NoEndFile=false ;
        }
        return 1;
    }
    return 0;
}

int findstrpos(Less  * tv )
{
    if ((tv->findstr).empty())
    {
        return 0;
    }
    else
    {
        for (int kk=(tv->pos)+1 ; kk<(tv->Now_Pring_Line); kk++)
        {
            if ((tv->PringLine[kk]).find(tv->findstr)!=string::npos)
            {
                (tv->pos)=kk;
                return 1;
            }
        }
        while(tv->NoEndFile)
        {
            tv->Before+=(tv->pos);
            tv->pos=0;
            readFile( tv ,Max_Input_Line );
            for (int kk=(tv->pos)+1 ; kk<(tv->Now_Pring_Line); kk++)
            {
                if ((tv->PringLine[kk]).find(tv->findstr)!=string::npos)
                {
                    (tv->pos)=kk;
                    return 1;
                }
            }
        }
        return 0;
    }
}



static void less_tv_win_command(Less *tv)
{
    char str[256], *p;
    int i, l = 0;
    wborder(tv->command, '|', '|', '-', '-', '+', '+', '+', '+');
    mvwprintw(tv->command, 1, 2, "find:str     find the string of str");
    mvwprintw(tv->command, 2, 2, "=Num         goto the pos of Num");
    mvwprintw(tv->command, 3, 2, "Enter command:");
    for (;;) {
        int c = wgetch(tv->command);
        wrefresh(tv->command);
        if (c == KEY_BACKSPACE || c == '\010' || c == '\177') {
            if (l>0){ --l;}
        }
        else if (c == KEY_ENTER || c == '\012' || c == '\015')
        {
            int  _beg;
            if (str[0] == '=') {
                _beg = strtol(str+1, &p, 10) - 1;
                if (_beg > -1) {
                    tv->pos = _beg+1;
                    return ;
                }
            } else {
                int ii=0 ;
                string A;
                for ( ii=0 ; ii<256 ; ii++)
                {
                    if (str[ii] == ':')
                    {
                        break ;
                    }
                    A=A+str[ii];
                }
                if (A=="find")
                {
                    string B=str;
                    tv->findstr=B.substr(ii+1);
                    tv->FindFlag=2;
                    (tv->FindFlag)+=findstrpos(tv);
                }
                return;
            }
        }
        else if (isgraph(c))
        {
            if (l < TV_MAX_GOTO) str[l++] = c;
        }
        else if (c == '\027') {l = 0;}
        else if (c == '\033') { return; }
        str[l] = '\0';
        for (i = 0; i < TV_MAX_GOTO; ++i) mvwaddch(tv->command, 3 , 16 + i, ' ');
        mvwprintw(tv->command, 3, 16, "%s", str);

    }
}

int tv_pl_func( Less  * tv )
{
    clear();
    int End=tv->mrow-1;
    for (int ii=0 ; ii<End ; ii++ )
    {
        int kk=ii+(tv->pos);
        int attr = 0 ,x=2;
        string tmp=(tv->PringLine)[kk];
        int B=tmp.length();
        if (tmp[0] == '#')
        {
            x=1;
        }
        else if (tmp.find("sage")!=string::npos )
        {
            x=8;
        }
        else if (tmp.find("document")!=string::npos )
        {
            x=10;
            attr|=A_UNDERLINE;
        }
        else if (tmp.find("tools")!=string::npos)
        {
            x=3;
        }
        else if (tmp.find("hewm")!=string::npos)
        {
            x=1;
            attr|=A_BLINK;
        }
        else if (tmp.find("-")!=string::npos )
        {
            x=6;
        }
        else if (tmp.find("help")!=string::npos  ||  tmp.find("Help")!=string::npos ||  tmp.rfind("elp")!=string::npos)
        {
            x=3;
        }
        else 
        {
            x=2;
        }

        string Num=Int2Str((tv->Before)+kk);
        if (!(tv->findstr).empty() &&  tmp.find(tv->findstr)!=string::npos)
        {
            string::size_type lagPos=0;
            map <string::size_type,bool> FindPos ;
            string::size_type length=(tv->findstr).length();
            while(tmp.find(tv->findstr,lagPos)!=string::npos)
            {
                lagPos=tmp.find(tv->findstr,lagPos);
                FindPos[lagPos]=true;
                lagPos+=length;
            }
            int attr_two=attr ;
            attr_two |= COLOR_PAIR(10);
            map <string::size_type,bool>:: iterator it=FindPos.begin();
            int head=0;
            attr |= COLOR_PAIR(x);
            for( ; it!=FindPos.end(); it++)
            {
                int sublength=(it->first)-head;
                string B=tmp.substr(head,sublength);
                attron(attr);
                mvaddstr(ii , (tv->col_shift+head),B.c_str()); 
                attroff(attr);
                attron(attr_two);
                mvaddstr(ii ,(tv->col_shift)+(it->first) , (tv->findstr).c_str());
                attroff(attr_two);
                head=(it->first)+length;
            }
            string B=tmp.substr(head);
            attron(attr);
            mvaddstr(ii , (tv->col_shift+head),B.c_str());
            mvaddstr(ii,0,Num.c_str());
            attroff(attr);
        }
        else
        {
            attr |= COLOR_PAIR(x);
            attron(attr);        
            if (tv->col_shift>0)
            {
                mvaddstr(ii,tv->col_shift,tmp.c_str());
            }
            else
            {
                int A=0-(tv->col_shift);
                A=(A>B)?B:A;
                string BB=tmp.substr(A);
                mvaddstr(ii,0,BB.c_str());
            }        
            mvaddstr(ii,0,Num.c_str());
            attroff(attr);
        }
    }
    string AA ;
    if (tv->FindFlag==0)
    {
        AA="? for Help ; q for leave ; / for command model";
    }
    else if (tv->FindFlag==2)
    {
        AA="No find : "+(tv->findstr);
    }
    else if ((tv->FindFlag)==3)
    {
        AA="Find : "+(tv->findstr);
    }

    else if (tv->FindFlag==4)
    {
         AA="The End";
    }
    move(End/2,tv->col_shift);
    mvaddstr(End,tv->col_shift,AA.c_str());
    refresh();
    tv->FindFlag=0;
    return  0;
}

int Less_init ( Less  * tv , int argc , char *argv[] )
{
    initscr();keypad(stdscr, TRUE); 
    clear();
    noecho();cbreak();
    tv->mrow = 24; tv->mcol = 80;
    tv->findstr="hewm";
    tv->FindFlag=0;
    tv->Before=0;
    getmaxyx(stdscr, tv->mrow, tv->mcol) ;
    tv->row_shift=0  ; tv->col_shift=tv->mcol/20 ;
    tv->Now_Pring_Line=0 ;
    tv->NoEndFile=true ;
    tv->pos=0 ;
    tv->whelp = newwin(15, 40, 5, 5);
    tv->command = newwin(5, 50, 10, 5);
    start_color();
    init_pair(0, COLOR_CYAN, COLOR_GREEN);
    init_pair(1, COLOR_BLUE, COLOR_BLACK);
    init_pair(2, COLOR_GREEN, COLOR_BLACK);
    init_pair(3, COLOR_YELLOW, COLOR_BLACK);
    init_pair(4, COLOR_WHITE, COLOR_BLACK);
    init_pair(5, COLOR_GREEN, COLOR_BLACK);
    init_pair(6, COLOR_CYAN, COLOR_BLACK);
    init_pair(7, COLOR_YELLOW, COLOR_BLACK);
    init_pair(8, COLOR_RED, COLOR_BLACK);
    init_pair(9, COLOR_BLUE, COLOR_BLACK);
    init_pair(10, COLOR_RED, COLOR_CYAN);
    init_pair(11, COLOR_RED, COLOR_GREEN);
    init_pair(12, COLOR_RED, COLOR_BLUE);


    string INPara=argv[0] ;
    string HelpFile=INPara+".Readme" ;
    if ( access(argv[2],0)==0 )
    {
        (tv->File).open(argv[2],ifstream::in);
    }
    else if ( access(HelpFile.c_str(),0)== 0 )
    {
        (tv->File).open(HelpFile.c_str(),ifstream::in);
    }
    else if ( access(argv[1],0)==0 )
    {
        (tv->File).open(argv[2],ifstream::in);
    }
    else
    {
        tv->PringLine[0]="Program: NewFxTools ";
        tv->PringLine[1]="Version: 1.02    hewm2008@gmail.com    2018-02-02";
        tv->PringLine[2]="";
        tv->PringLine[3]="\t\tUsage:";
        tv->PringLine[4]="\t\tFatools        Tools For Fasta";
        tv->PringLine[5]="\t\tFqtools        Tools For Fastq";
        tv->PringLine[6]="\t\tFormtools     Tools For Form convert";
        tv->PringLine[7]="\t\tFiletools     Tools For Specified File";
        tv->PringLine[8]="";
        tv->PringLine[9]="\t\tHelp          Show this help";
        tv->PringLine[10]="";
        string A="Can't Find the help document: "+HelpFile;
        tv->PringLine[11]=A;
        tv->PringLine[12]="\t\t";
        tv->PringLine[13]="please contact hewm2008@qq.com hewm2008@gmail.com to get this doc";
        tv->PringLine[14]="See https://github.com/BGI-shenzhen/Reseqtools  ";
        tv->PringLine[15]="See http://code.google.com/p/reseqtools ";
        tv->PringLine[16]="Join Communication & Discussion QQ Group: 125293663 ";
        tv->NoEndFile=false ;
        tv->Now_Pring_Line=16;
    }
    readFile( tv , 168 );
    tv_pl_func(tv);
    readFile( tv ,Max_Input_Line);
    return 1;
}

int Less_des(Less  * tv )
{
    delwin(tv->whelp);
    delwin(tv->command);
    tv->File.close();
    endwin();
    return 0 ;
}

void Less_win_help(Less *tv) {
    int r = 1;
    WINDOW *win = tv->whelp;
    wborder(win, '|', '|', '-', '-', '+', '+', '+', '+');
    mvwprintw(win, r++, 2, "        -=-    Help    -=- ");
    r++;
    mvwprintw(win, r++, 2, "?          This window");
    mvwprintw(win, r++, 2, "Arrows     Small scroll movement");
    mvwprintw(win, r++, 2, "a,s,w,d    Small scroll movement");
    mvwprintw(win, r++, 2, "A,S,W,D    Large scroll movement");
    mvwprintw(win, r++, 2, "ctrl-H     Scroll 1k left");
    mvwprintw(win, r++, 2, "ctrl-L     Scroll 1k right");
    mvwprintw(win, r++, 2, "space      Scroll one screen");
    mvwprintw(win, r++, 2, "n          next findstr");
    mvwprintw(win, r++, 2, "/          for command model");
    mvwprintw(win, r++, 2, "q          Exit");
    wrefresh(win);
    wgetch(win);
}

void Less_tv_loop(Less *tv,  int argc , char *argv[])
{
    /*
       struct timeval tv; 
       tv.tv_sec = 5;
       tv.tv_usec = 0;
       */
    while(1)
    {
        ////*//
        int c = getch();
//        if (c)
 //       {
            //            int c = getch();
            switch (c)
            {
                case 'h':
                case 'H':
                case '?': Less_win_help(tv); break;
                case '\033':
                case 'Q':
                case 'q': goto end_loop;
                case KEY_LEFT:
                case 'a':
                          --(tv->col_shift); break;
                case KEY_RIGHT:
                case 'd':
                case 'l': ++(tv->col_shift); break;
                case KEY_SLEFT:
                case 'A':(tv->col_shift) -= 20; break;
                case KEY_SRIGHT:
                case 'D':
                case 'L': (tv->col_shift) += 20; break;
                case '\010': (tv->pos) -= 1000; break;
                case '\014': (tv->pos) += 1000; break;
                case KEY_NPAGE:
                case ' ': (tv->pos) += tv->mrow; break;
                case KEY_UP:
                case 'w':
                case 'j': --(tv->pos); break;
                case 'n':tv->FindFlag=2;(tv->FindFlag)+=findstrpos(tv); break;
                case '/':less_tv_win_command(tv); break;
                case KEY_DOWN:
                case 's':
                case 'k': ++(tv->pos); break;
                case 'W': (tv->pos)-=20;break;
                case 'S': (tv->pos)+=20;break;
                case KEY_BACKSPACE:
                case KEY_PPAGE:
                case '\177': (tv->pos) -= tv->mrow ; break;
                case KEY_RESIZE: getmaxyx(stdscr, tv->mrow, tv->mcol); break;
                default: continue;
            }
     //   }
/*///
        else
        {
            (tv->pos)++;
            sleep(1);
        }
        //*///
        ///*///
        //      if (tv->row_shift < 0) {tv->row_shift = 0;}
        if (tv->pos<0){tv->pos=0;}
        else if ((tv->pos)<(tv->Before))
        {
            int _beg=(tv->pos);
            Less_init ( tv , argc , argv );
            tv->pos = _beg;
        }
        else if ( (tv->pos)>Max_Input_Line  && (tv->NoEndFile))
        {
            tv->Before+=(tv->pos);
            tv->pos=0;
            readFile( tv ,Max_Input_Line );
        }
        else if ( (tv->pos)>(tv->Now_Pring_Line)  && (!(tv->NoEndFile)) )
        {
            tv->pos=(tv->Now_Pring_Line);
            tv->FindFlag=4;
        }
        tv_pl_func(tv);
        tv->FindFlag=0;
        if (tv->NoEndFile)
        {

        }
    }
end_loop:
    return;
}

int help_main ( int argc , char *argv[]  )
{
    Less  * tv =new Less ;
    Less_init(tv ,argc,argv);
    Less_tv_loop(tv, argc ,  argv );
    Less_des (tv);
    return 0;
}

#endif //HelpLess_H_
////////////////////////swimming in the sea & flying in the sky //////////////////


