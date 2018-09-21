#!/usr/bin/perl-w
use strict;
#explanation:this program is edited to 
#edit by HeWeiMing;   Mon May 14 16:38:54 HKT 2018
#Version 1.0    hewm@genomics.org.cn 

die  "Version 1.0\t2018-05-14;\nUsage: $0 <InPut><Out>\n" unless (@ARGV ==1);

#############Befor  Start  , open the files ####################

open (IA,"$ARGV[0]") || die "input file can't open $!";

#open (OA,">$ARGV[1]") || die "output file can't open $!" ;

print "void More_HelpFM()\n";
print "{\n";
################ Do what you want to do #######################
print "cout<<\"\\n\\\n";
	while(<IA>) 
	{ 
		chomp ; 
		if  ($_ eq "")
		{
			print "\\n\\\n";
		}
		else
		{
	    print "\\t\\t$_\\n\\\n"	;
	}
	}
print "\\n\";\n";
#print "return 1 ;\n";
print "}\n";
close IA;
#close OA ;

######################swimming in the sky and flying in the sea ###########################
