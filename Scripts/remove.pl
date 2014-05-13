#!/usr/bin/perl
use strict;
use warnings;


#open(IN,"</home/alexx/Fungal3endprocessing/InterProScan_results/MGG_00680/MGG_00680_1Aspergillus_fumigatus.fasta_InterPro_out") or die "open read failed: $!";
#open(OUT, ">/home/alexx/Fungal3endprocessing/MGG_00680_1Aspergillus_fumigatus.fasta_InterPro_out2") or die "open write failed: $!";

open(IN, shift(@ARGV)) or die "open read failed: $!";

my $flag = 1;

while (my $line = <IN>){
	#chomp $line;
 	if($line =~ m/##FASTA/){
 	    $flag = 0;
 	}
        if ($flag == 1) {
            print $line;
        }
        
}   
