#!/usr/bin/perl -w
use strict;
use Bio::DB::SeqFeature::Store;
use Bio::SeqFeature::Generic;

#my $ID = <stdin>;
my $ID = $ARGV[0];									#The ID is introduced when running the script in the terminal
my $ID2 =$ARGV[1];
#my $ID2 = $ARGV[1];
 # Open the sequence database
 my $db      = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                                -dsn     => 'dbi:mysql:Fungal3endprocessing',
						-user => 'root',
						-password => '');


my @features = $db->features(-seqid => $ID, -type => "protein_match");                  #Give me only the IPScan features (that's why I specify "protein_match")
my %featurehash;                                                                        #The map built-in function (next line) allows me to build a new list from another list while modifying the elements.
map {$featurehash{$_.$_->start} = $_} @features;                                        #'Map' se usa para ejecutar el código entre {} y que se aplique a todas las posiciones del array, y en este caso se usa para: to solve the triplication problem (OJO! Specify the name+$_->start, to avoid eliminating same features with the same name but different start!!) The idea is: if the start is the same, store it in the same place (that way we avoid repetition of the features)
@features = map {$featurehash{$_}} keys %featurehash;                                   #Meter los identificadores del hash (protein name) en @features usando 'map'. The 'key' function returns an array of all the keys of the %feature hash. This is stored in the @features array.
my $thisfeature = $features[0];                                                         #For the reference protein (MGG_... or S.cerevisiae protein)
                                                                                                                                                                        
foreach my $thisfeature(@features){

    my $name = $thisfeature->name;
    my $source = $thisfeature->source;
    my $start = $thisfeature->start;
    my $end = $thisfeature->end;
}

my @otherfeatures = $db->features(-seqid => $ID2, -type => "protein_match");
my %otherfeaturehash;                                                       
map {$otherfeaturehash{$_.$_->start} = $_} @otherfeatures;                              #To solve the triplication problem (OJO! Specify the name+$_->start, to avoid eliminating same features with the same name but different start!!
@otherfeatures = map {$otherfeaturehash{$_}} keys %otherfeaturehash;
my $thatfeature = $features[1];                                                         #For the top BLAST hit (in this case)


foreach my $thatfeature(@otherfeatures){

    my $name = $thatfeature->name;
    my $source = $thatfeature->source;
    my $start = $thatfeature->start;
    my $end = $thatfeature->end;
       
}


my @sorted_features = sort {$a->start <=> $b->start || $a->end <=> $b->end || $a->name cmp $b->name} @features;             #Sort (ordenar) first by the start, then the end (in case the start is the same?), and then order it also by name (in case the end is the same?)
my @namefeature = "";
foreach (@sorted_features){
    push(@namefeature, $_->name);                                                                                           #store all the names into the @namefeature array
    print "Query ", $_->name, "\t", $_->start, "\t", $_->end, "\n";
}
my @sorted_otherfeatures = sort {$a->start <=> $b->start || $a->end <=> $b->end || $a->name cmp $b->name} @otherfeatures;   #Sort first by the start, then the end, and then order it also by name
my @nameotherfeature = "";
foreach (@sorted_otherfeatures){
    push(@nameotherfeature, $_->name);                                                                                       #store all the names into the @nameotherfeature array
    print "Candidate ", $_->name, "\t", $_->start, "\t", $_->end, "\n";
}



#sub compare {
#    my ($sorted_features, $sorted_otherfeatures)= @_;
#    my @sorted_features = @$sorted_features;
#    my @sorted_otherfeatures = @$sorted_otherfeatures;
#if (@sorted_features ~~ @sorted_otherfeatures) {                                                                                      #compare the number of domains.(scalar gives you the number of elements in the array?)
#    print "Same numbers of domains", "\n";
##
#    if (@namefeature eq @nameotherfeature) {                                                                                   #check the order of the domains 
#        print "Same domain order along the protein when compared (order is the same in both)", "\n";
#    }else{
#        print "Different domain order along the protein when compared (order is not the same)", "\n";
#    }
#}
##
#else {                                                                                          
#    print "Not the samber number of domains", "\t", "Query domains: ", scalar @namefeature, "\t", "Candidate domains: ", scalar @nameotherfeature, "\n";
#}
#}
#compare (\@sorted_features, \@sorted_otherfeatures);


#my $stringnames = join (',', @namefeature);                                                                                  # connect all the elements (names) of @namefeature into a single string    
#my $stringnamesotherfeature = join (',', @nameotherfeature);                                                                 # connect all the elements (names) of @nameotherfeature into a single string
if (@namefeature ~~ @nameotherfeature) {                                                                                      #compare the number of domains.(scalar gives you the number of elements in the array?)
   print "Same numbers of domains", "\n";
#
    if (@namefeature eq @nameotherfeature) {                                                                                  #check the order of the domains 
        print "Same domain order along the protein when compared (order is the same in both)", "\n";
   }else{
       print "Different domain order along the protein when compared (order is not the same)", "\n";
    }
}
#
else {                                                                                          
    print "Not the samber number of domains", "\t", "Query domains: ", scalar @namefeature, "\t", "Candidate domains: ", scalar @nameotherfeature, "\n";
    if (scalar @namefeature < scalar @nameotherfeature) {                                                                   #if the query has less domains than the candidate (scalar context)
##        #(my $subname = $stringnames) =~ s/,/.*/g;                                                                           #get rid of the commas in the name list and substitute it with *
        if (@namefeature eq @nameotherfeature) {                                                                         
           print "Query domain set is a subset of candidate domain set", "\n";                                             #The query domains are also present in the candidate protein, even if the query protein has less domains
       }else {
            print "Query domain set is not a subset of candidate domain set (some candidate domains lacking in the query protein set) ", "\n";  #Some candidate domains are not present in the query domain set
        }
    }
    if (scalar @namefeature > scalar @nameotherfeature) {                                                                   #if the candidate has less domains than the query (scalar context)
##        #(my $subname = $stringnamesotherfeature) =~ s/,/.*/g;
        if (@namefeature eq @nameotherfeature) {    
           print "Candidate domain set is a subset of query domain set", "\n";                                             #The candidate domains are also present in the query protein, even if the candidate protein has less domains
        }else{
            print "Candidate domain set is not a subset of query domain set (some query domains lacking in the candidate domain set)", "\n";    #Some query domains are not present in the candidate domain set
        }
  }   
}