#!/usr/bin/perl -w
use strict;

my $file=<"/home/alexx/Fungal3endprocessing/BLASTp/Global_BLASTp/BLASTp_parser_results/0.0_results">;
my @lines = split "\n", $file;
foreach my $line(@lines){
            if ($file =~ m/As/i) {
            print STDERR "Yes!", "\n";
        }
}