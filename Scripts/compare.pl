

#!/usr/bin/perl

# Author: Olivia Ho-Shing. BIO 343: Laboratory Method in Genomics. 20 October 2009. 

# SUMMARY: This program was created to accept two text files containing amino acid sequences in FASTA format. 
# Sequences in the file designated as the query proteome will be compared with sequences in the other file, 
# designated as the "database". Each query sequence is compared pairwise to each database sequence by calling
# bl2seq.exe, a local BlastP alignment. If any sequences from the query proteome have significant alignments 
# (E-value is less than or equal to the given threshold value), the sequence is considered conserved between 
# the two proteomes, and the name of the sequence as it is labeled in the query proteome will be printed to a 
# file for conserved sequences. Query sequences with no significant alignments to any database sequences is 
# considered unique to the query proteome, and the name is printed to a file for unique sequences. 

use strict; 
use warnings; 
use CGI;                        # Loads the package to use CGI interface for BLAST on the web. 
# use POSIX;                    # Loads the package to use POSIX on the web. May not be necessary for this program.  


####################################################################################################################################################################
# USER INPUT goes here. Designate your file names and an E-value threshold for the BlastP alignment within the quotations (no spaces) on the right side of 
# the equal sign. Result file names should be specific identifiers for your project. Plain text files are recommended for the input files containing the proteomes. 
# Input files must have sequences in FASTA format. All alignments with an E-value greater than the designated E-value threshold will not be a signficant hit. 


    my $queryProteome          = "/home/alexx/Fungal3endprocessing/FASTA_fungal_genomes/Multiple_genomes/Ashbya_gossypii.fasta";                     # Input name (case-sensitive) of query file in FASTA format
    my $compareProteome        = "/home/alexx/Fungal3endprocessing/FASTA_fungal_genomes/Multiple_genomes/Magnaporthe_oryzae.fasta";                    # Input name (case-sensitive) of database file in FASTA format
    
    my $uniqueResultsFile      = "/home/alexx/Fungal3endprocessing/FASTA_fungal_genomes/Multiple_genomes/unique.txt";        # Create a specific name of a .doc, .txt, or .xls file to list unique protein names
    my $conservedResultsFile   = "/home/alexx/Fungal3endprocessing/FASTA_fungal_genomes/Multiple_genomes/conserved.txt";     # Create a specific name of a .doc, .txt, or .xls file to list conserved protein names
    
    my $threshold              =   0.0001;                                    # Designate an E-value threshold (recommended <= 0.001)
    

# Please do not alter anything below this line, unless you have reason to validate the program output and debug. If so, follow directions where marked "PROGRAMMER". 
####################################################################################################################################################################



  #--------------------------------------------------------------------------------------------------------------------------------------#
  #--------------------------------------------------------------------------------------------------------------------------------------#
  # PROGRAMMER can uncomment the following four lines, 2 lines in the middle of this code (near line 135), and the last line in this code 
  # (near line 170) in order to create a default text file "testOutput.txt" to see the output of bl2seq BlastP alignment. 
  
  my $testbl2Output; 												  # Creates a variable to see Blast output after program runs.
  my $OUTPUT_HANDLE;												  # Creates file handle to write an output file. 
  $testbl2Output = 'testbl2Output.txt'; 							  # Names file testbl2Output.txt. 
  open ($OUTPUT_HANDLE, ">", $testbl2Output ) or die $!; 			  # Opens output text file for writing, or end program. 
  
  #--------------------------------------------------------------------------------------------------------------------------------------#
  #--------------------------------------------------------------------------------------------------------------------------------------#



  my $format = "%.1e";                                                # Defines the format for the threshold value given - this format matches the format in the bl2seq output.
  my $formatthresh = sprintf($format, $threshold);                    # Formats the user input threshold value. 


## -- Open the proteome files for reading, and open the output files for writing

  open(my $fh, '<', $queryProteome) or die $!;                        # Opens query proteome file or end program. 
  my @qWords = <$fh>;                                                 # Reads the file. 
  my $qProt = join("",@qWords);                                       # Joins the sequences and then split them into an array @myProtein. 
  my @myProtein = split(/>/, $qProt);                                 # Separates the sequences by the FASTA format ">". Removes the > from beginning of array elements. 
  #print "\nScalar myProtein: " , scalar(@myProtein) , "\n";          # PROGRAMMER can uncomment to check the number of array items. 
 
  open(my $fn, '<', $compareProteome) or die $!;                      # Opens database proteome file or end program. 
  my @hWords = <$fn>;                                                 # Reads the file. 
  my $hProt = join("",@hWords);                                       # Joins the sequences and then split them into an array @theirProtein. 
  my @theirProtein = split(/>/, $hProt);                              # Separates the sequences by the FASTA format ">". Removes the >. 
  #print "\nScalar theirProtein: " , scalar(@theirProtein) , "\n";    # PROGRAMMER can uncomment to check the number of array items.

  my $UNIQUE_HANDLE;												  # Creates file handle to write into unique names file. 
  open ($UNIQUE_HANDLE, ">", $uniqueResultsFile ) or die $!;          # Opens unique results text file for writing or end program.  

  my $CONSERVED_HANDLE;												  # Creates file handle to write into conserved names file. 
  open ($CONSERVED_HANDLE, ">", $conservedResultsFile ) or die $!;    # Opens conserved results text file for writing or end program. 


## -- Loop through the query sequences and database sequences, and store in temporary files for pairwise BlastP compare

  my $output; 														  # Declares variable to hold BlastP alignment output 
  
  my $a = 1;                  										  # Declares variable to count array element position of the @myProtein (query) array. The first element position (0) is empty.  
  
  while ( $a < scalar(@myProtein) ) {								  # Begins while loop to loop through the query sequences by the array position. While loop ends after calling the last array element. 
        
        $output = ""; 												  # Ensures the output variable is empty at the beginning of each loop when a new query sequence is to be aligned. 
        my $myName;       									          # Constructs variable to store protein's name. 
        $myProtein[$a] =~ m/(.*)\W(.*)/;                              # Uses regular expression to find in array element a string (the name), a whitespace (the newline) and another string (the sequence).
        															  # Captures these two strings. 
        $myName = $1;  												  # Defines and saves the first captured string as the sequence name.
     	my $mySequence = $2; 										  # Defines and saves the second captured string as the sequence.  
     
     
        my $tempOurProtein;   										  # Creates a temporary file within this while loop for the single captured query sequence. 
		my $TEMP_OURPROTEIN;  										  # Creates file handle for temporary file. 
		$tempOurProtein = 'tempOurProtein.txt';  					  # Calls the file tempOurProtein.txt. 
        open ($TEMP_OURPROTEIN, ">", $tempOurProtein ) or die $!;	  # Opens temporary file for writing or end program.  
		print $TEMP_OURPROTEIN ">$myName\n"; 						  # Writes the sequence name at top of temporary file for checking by user. 
		                                                              # At end of program, should contain only the last sequence from the query file.
        print $TEMP_OURPROTEIN $mySequence;                           # Writes the sequence after a newline into temporary file. 
		close ( $TEMP_OURPROTEIN );   								  # Closes the temporary file for writing. 


        my $d = 1;                 									  # Declares variable to count array element position of the @theirProtein (database) array. First element position (0) is empty.   
   
        while ($d < scalar(@theirProtein) ) {						  # Begins nested while loop to loop through the database sequences by array position. While loop ends after last array element. 
   
              my $theirName;                                          # Constructs variable to store protein's name. 
     		  $theirProtein[$d] =~ m/(.*)\W(.*)/;                     # Uses regular expression to find in array element a string (the name), a whitespace (the newline) and another string (the sequence).
     		 													      # Captures these two strings. 
     		  $theirName = $1;  									  # Defines and saves the first captured string as the sequence name.
     		  my $theirSequence = $2;  								  # Defines and saves the second captured string as the sequence.
     		  
     		  my $tempTheirProtein;    					              # Creates a temporary file within the nested while loop for the single captured database sequence. 
			  my $TEMP_THEIRPROTEIN;  								  # Creates file handle for temporary file. 
			  $tempTheirProtein = 'tempTheirProtein.txt';  			  # Calls the file tempTheirProtein.txt.  
			  open ($TEMP_THEIRPROTEIN, ">", $tempTheirProtein )      # Opens temporary file for writing or end program. 
			    or die $!; 
			  print $TEMP_THEIRPROTEIN ">$theirName\n"; 			  # Writes the sequence name at top of temporary file for checking by user.
			  print $TEMP_THEIRPROTEIN $theirSequence;  			  # Writes the sequence after a newline into temporary file. 
		      close ( $TEMP_THEIRPROTEIN );    						  # Closes the temporary file for writing. 

        ## -- Use bl2seq.exe to call BlastP to BLAST the query sequence to the database sequence

			  $output = $output . `bl2seq -i $tempOurProtein -j $tempTheirProtein -p blastp -e $threshold`;
			  														  # The above line uses bl2seq.exe to perform the pairwise alignment between the two temporary files created.
			  														  # The output is stored in the output variable, and as the database sequences are looped through, all outputs are concatenated. 
			  														  # The entire concatenated output is cleared once the database array has been completely looped through, ending the nested while loop. 
			  														  
			 
			  #--------------------------------------------------------------------------------------------------------------------------------------#
			  #--------------------------------------------------------------------------------------------------------------------------------------#
			  # PROGRAMMER can uncomment the following two lines to print the default test file "testbl2Output.txt". Next, uncomment the last line. 
			  
			  print $OUTPUT_HANDLE $output;						  # Prints the concatenated output to a text file, noting the loop iterations. 
			  print $OUTPUT_HANDLE "\n --------------------- End loop iteration --------------------- \n"; 
			  
			  #--------------------------------------------------------------------------------------------------------------------------------------#
			  #--------------------------------------------------------------------------------------------------------------------------------------#


			  $d++;   												  # Increments to next amino acid sequence in database array
        }       													  # Ends database loop when reach the end of the database array @theirProtein


## -- Extract number of E-value hits from output, and sort into unique protein file and conserved protein file 

        my $regex = 'Number of sequences better than ' . $formatthresh . ': ' . '[^0]';     
                                                                      # Searches for the number of hits above the E-value threshold ($formatthresh). 
                                                                      # Will find this line if the number of hits is not 0
                                                                      
        if ( $output =~ m/$regex/gi ) {                               # If this line does occur anywhere in the output, 
        print $CONSERVED_HANDLE "$myName\n"; }                        # Then the query matches a part of the database, and the protein is conserved. Saves to Conserved Results file
        else {                                                        
        print $UNIQUE_HANDLE "$myName\n"; }                           # Or else it is unique, and saves to Unique Results File

        $a++;                                                         # Increments to next amino acid sequence in query array
  }                                                                   # Ends query loop when it's reached the end of the query array @ourProtein. 


## -- Close files after writing is complete. Files are found in same directory as program 

  close ( $CONSERVED_HANDLE );                                        # Closes the conserved amino acids file for writing. 
  close ( $UNIQUE_HANDLE );                                           # Closes the unique amino acids file for wriring. 

  #-------------------------------------------------------------------#
  #-------------------------------------------------------------------#
  # PROGRAMMER can uncomment the follow line to print the default test. 
  
  close ( $OUTPUT_HANDLE );                
  
  #-------------------------------------------------------------------#
  #-------------------------------------------------------------------#


 
