#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;


      
my @files=<"/home/alexx/Fungal3endprocessing/BLASTp/Global_BLASTp/*out">;       #Go through all the *out files in the folder
    foreach my $file (@files){

        my $input = new Bio::SearchIO(-format => 'blast', 
                                      -file   => $file);


            while( my $result = $input->next_result ) {

                while( my $hit = $result->next_hit ) {

                    while( my $hsp = $hit->next_hsp ) {

                            #if( $hsp->length('total') > 50 )
                            #if ($hsp->evalue == 0.0) {
                            if ( ($hsp->evalue >= 1e-10) && ($hsp->evalue <= 1e-5)) {          
        
print "Query=", $result->query_name, "\t", "Hit=", $hit->name, "\t", "Description=", $hit->description, "\t", " E_value=", $hsp->evalue, "\n";
                                    
                                        open (OUT, ">>/home/alexx/Fungal3endprocessing/BLASTp/Global_BLASTp/BLASTp_parser_results/e-10_e-5_hits_results") or die "Can't open $! \n";                #print to a file                             
                                        print OUT "Query=",  $result->query_name, "\t", " Hit=",  $hit->name, "\t", " Description=",  $hit->description, "\t", " E_value=", $hsp->evalue, "\n","\n";
            }                 
        }
      }
    }  
  }    
    
