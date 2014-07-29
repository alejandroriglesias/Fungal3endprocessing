#!/usr/bin/perl 
use strict;
use warnings;

my %hash1;
my %hash2;
my ($protein, $species, $evalue);
my $threshold = $ARGV[0];                               #Introduce the e-value in the command-line directly
my $threshold2 = $ARGV[1];                              #Introduce the other e-value in the command-line (if needed)

open (IN, "</home/alexx/Fungal3endprocessing/BLASTp/Global_BLASTp/BLASTp_parser_results/tabular_outputs/Formatted_tabular_outputs/tabular_results_tab_out" ) || die "failed $!\n";
my $line = <IN>;
my @lines = split "\n", $line;
#chomp $line;
while (<IN>){
    chomp $_;
    if ($_ eq "") {                                     #Jump the blank lines between the lines to avoid compilation errors
        next;
    }
    my @tab = split "\t", $_;   
    $protein=$tab[0];                                   #First element of the tab array (protein_ID)
    $species=$tab[1];                                   #Second element of the tab array (fungal genome species)
    $tab[2] =~ s/^\s+|\s+$//g;                          #Delete the whitespace that appears before the e-value (weird)
    $evalue=$tab[2];                                    #Third element of the tab array (e-value)
    #print $species;
    #print $evalue;
    if (!defined($hash2{$protein})) {
        $hash2{$protein}=$protein; 
    }
    if ($evalue <=$threshold) {
    #}
    
    #if (($evalue >= $threshold) && ($evalue <=$threshold2)) {                               #This is where you can select the e-value range you want(in the shell)
        
        #print $hash2{$protein};   
       $hash1{$species}{$protein}= $evalue;             #Obtain the e-value of the selected elements (protein and species)
    
       
       #print $hash1{$species}{$protein};
    }
    
    
    
}

close IN;
#print "\n";
print "species\t";                                      
foreach my $prot  (sort (keys %hash2)){
    print $hash2{$prot},"\t";
}
print "\n";
foreach my $spec (sort (keys %hash1)){
    print $spec, "\t";
    foreach my $prot (sort (keys %hash2)){
           if (!defined($hash1{$spec}{$prot})) {
            print "NO", "\t";
           }
           else{
            print $hash1{$spec}{$prot}, "\t";
            }
    }
    print "\n";
}