#!/usr/bin/perl


use strict;
use warnings;

sub thouse{
  my ($num) = @_;
  $num = reverse $num;     # reverse the number's digits
  $num =~ s/(\d{3})/$1,/g; # insert comma every 3 digits, from beginning
  $num = reverse $num;     # Reverse the result
  $num =~ s/^\,//;         # remove leading comma, if any
  return $num;
}


if($#ARGV == -1){
  die "Usage\ntab2latex.pl [options] fileWithTabSeparatedData\nOptions:\n\t-s\t1,2,3Use thousand separator on fields 1,2,3\n";
}

my $thousSep=0;
my @thousSepFields=0;

for(my $i=0;$i<$#ARGV;$i++){
  if($ARGV[$i] eq "-s"){
    $thousSep=1;
    @thousSepFields=split(",",$ARGV[$i+1]);
  }
}

my $fileToUse=$ARGV[$#ARGV];


print "\\documentclass[a4paper,10pt]{article}\n";
print "\\usepackage[landscape]{geometry}\n";
print "\\usepackage[dvips]{graphicx,epsfig}\n";
print "\\usepackage{color}\n";
print " \n";
print "\n";
print "\\begin{document}\n";



my $firstLine=1;
open(FILE,$fileToUse) or die "Cannot open ".$fileToUse."\n";
while(my $line = <FILE>){
  my @array=split("\t",$line);
  if($firstLine){
    print "\\begin{tabular}{";
    print "l"x($#array+1);
    print "}\n";

    $firstLine=0;
  }

  if($thousSep){
    for(my $i=0;$i<=$#array;$i++){
      foreach my $index (@thousSepFields){
	if(($i+1) == $index){
	  $array[$i]=thouse($array[$i]);
	  last;
	}
      }

    }
  }

  for(my $i=0;$i<=$#array;$i++){
   $array[$i] =~ s/\%/\\%/g;
  }

  my $string= join(" & ",@array);
  chomp($string);
  print ($string." \\\\ \n");
}
close(FILE) or die "Cannot close ".$fileToUse."\n";


print "\\end{tabular}\n";
print "\\end{document}\n";
