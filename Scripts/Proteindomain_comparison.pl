#!/usr/bin/perl -w
use strict;
use Bio::DB::SeqFeature::Store;
use Bio::SeqFeature::Generic;
use experimental 'smartmatch';

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


my %protein;
$protein{$ID} = \@sorted_features;
$protein{$ID2} = \@sorted_otherfeatures;
#print @sorted_features , "\n";
#print @sorted_otherfeatures, "\n";

if ( scalar(@{$protein{$ID}}) == scalar(@{$protein{$ID2}})) { # @{$protein{$ID}} means "take the list memory address referred to in protein{$ID} and turn it into a list"
    print "same number of domains\n";
} else {
    print "different number of domains\n";
}

my $longest = (scalar(@{$protein{$ID}}) > scalar(@{$protein{$ID2}}))?$ID:$ID2;  # the CONDITION ? X : Y  operator says "condition true, return X, otherwise return Y"
my $shortest = ($longest eq $ID)?$ID2:$ID;
###
my @shortestlist = map {$_->name} @{$protein{$shortest}};
my @longestlist = map {$_->name} @{$protein{$longest}};


foreach my $element(@{$protein{$longest}}){
    
    if ($element->name ~~ @shortestlist) {                                                                           # if the current element is a member of @list2
        my( $position_longest )= grep {$longestlist[$_] eq $element->name} 0..scalar(@longestlist);
       my( $position_shortest )= grep {$shortestlist[$_] eq $element->name} 0..scalar(@shortestlist);
       unless ($position_longest == $position_shortest) {
         print "domain $element is in position $position_longest in $longest protein; in position $position_shortest in $shortest protein\n";
       }
   } else {
       print "domain $element is not part of protein $shortest\n";
   }
}

#foreach my $element(@{$protein{$shortest}}){
#    
#if ($element->name ~~ @longestlist) {                                                                           # if the current element is a member of @list1
#        my( $position_shortest )= grep {$shortestlist[$_] eq $element->name} 0..scalar(@shortestlist);
#       my( $position_longest )= grep {$longestlist[$_] eq $element->name} 0..scalar(@longestlist);
#       unless ($position_shortest == $position_longest) {
#         print "element $element is in position $position_shortest in $shortest protein; in position $position_longest in $longest protein\n";
#       }
#   } else {
#       print "element $element is not part of protein $longest\n";
#   }
#}
