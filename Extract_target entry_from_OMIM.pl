#! /usr/bin/perl
use LWP::Simple;
use Getopt::Std;
getopts "f:e:h:";
if ((! defined $opt_f) &&(!defined $opt_e)||((defined $opt_f) && (defined $opt_e))|| (defined Sopt_h)){
die
"*********************************
Two methods can be choosed to extract information from OMIM:
   Usagel:Extract_from_OMIM.pl [-e entry]
   Usage2:Extract_fromOMIM.pl [-f filename];
     -h: help and usage. 
	 -e: single OMIM entry
	 -f:a file contains more than one oMIM entryn
	 *******************\n";
}
my @information=("Animal Model", "Biochemical Features","Clinical Features", 
   "Clinical Management","Cloning","Cytogenetics", "Description","Diagnosis",
   "volution", "Gene Family"."Gene Function", "Gene Structure", "GeneTherapy", 
   "Genetic Variability","Genotype","Genot ype/Phenotype Correlations", 
   "Heterogeneity", "History", "Inheritance" ."Mapping"."Molecular Genetics", 
   "Nomenclature","OthezFeatures", "Pathogenesis","Phenotype"."Population Genetics");
   
if(defined $opt_e) {
   my $content =
   get ("http://api.omim.org/api/entry?mimNumber=$opt_e&include=text&apiKey=NxXfhjnoT1GKYAk2wgHe0g");
   die "Couldn't get OMIM data" unless defined $content;
   open (OUT, ">Entry_$opt_e.txt");
   for $i(@information) {
   if($content=~/\<textSectionTitle\>$i\<\/textSectionTitle\>\n\<textSectionContent\>(.*?)\<\/textSectionContent\>/s) {
   my $match=$1;
   $match=~s/\n//g;
   print OUT $i."\t".$match. "\n";
   }
   else{
   print OUT $i."\t"."NA"."\n";
   }
  }
   close OUT;
}

if (defined $opt_f){
   open ENTRY, $opt_f;
   while (<ENTRY>) {
     chomp;
     my @item=split/\t/;
     push @entrvID. $item[0]
   if($item[0]=~/^\d+$/);
   }
   print "Couldn 't get OMIM data of \"Entry: $ID\""."\n" unless defined 
   Şcontent;
   next unless defined $content;
   open (OUT, ">Entry_$ID.txt");
   for $i(@information){
       if ($content=~/\<textSectionTitle\>$i\<\/textSectionTitle\>\n\<textSectionContent\>(.*?)\<\/textSectionContent\>/s){
         my $match=$1;
         $match=~s/\n//g;
         print OUT $i."\t".$match."\n";
        }
        else{
         print OUT $i."\t"."NA"."\n";
            }
        }
    close OUT;
}


