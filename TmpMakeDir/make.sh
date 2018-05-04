#!/bin/sh
#$ -S /bin/sh
#Version1.0	P_bc_rd@genomics.org.cn	2018-02-2
echo Start Time : 
date

echo  g++ --std=c++11 	-g	-O3	-o	FxTools	FxTools.cpp  -lhts	-lncurses	-lm	-lpthread	-lboost_thread	-lboost_system	-lz	-I	/ifshk4/BC_PUB/biosoft/newblc/01.Usr/include/	-L	/ifshk4/BC_PUB/biosoft/newblc/01.Usr/lib -L./src/lib  -static
g++	 --std=c++11 -g	-O3	-o	FxTools	FxTools.cpp  -lhts	-lncurses	-lm	-lpthread	-lboost_thread	-lboost_system	-lz	-I	/ifshk4/BC_PUB/biosoft/newblc/01.Usr/include/	-L	/ifshk4/BC_PUB/biosoft/newblc/01.Usr/lib -L./src/lib   -static
# gcc > 4.9  

echo End Time : 
date
