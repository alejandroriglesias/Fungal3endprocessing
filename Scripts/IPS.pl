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
my %features;                                                                 
my %junk = map {$features{$_.$_->start} = $_} @features;                                #To solve the triplication problem (OJO! Specify the name+$_->start, to avoid eliminating same features with the same name but different start!!)
@features = map {$features{$_}} keys %features;
my $thisfeature = $features[0];                                                         #For the reference protein (MGG_... or S.cerevisiae protein)
                                                   
#@features = $db->get_features_by_location($ID);

foreach my $thisfeature(@features){

    my $name = $thisfeature->name;
    my $source = $thisfeature->source;
    my $start = $thisfeature->start;
    my $end = $thisfeature->end;
        #print "Query_domain_name: ", $thisfeature->name, "\t", "Query_Source: ", $thisfeature->source, "\t", "Query_domain_start: ", $thisfeature->start, "\t", "Query_domain_end: ", $thisfeature->end, "\n";
}

my @otherfeatures = $db->features(-seqid => $ID2, -type => "protein_match");
my %otherfeatures;                                                       
my %junk2 = map {$otherfeatures{$_.$_->start} = $_} @otherfeatures;                     #To solve the triplication problem (OJO! Specify the name+$_->start, to avoid eliminating same features with the same name but different start!!
@otherfeatures = map {$otherfeatures{$_}} keys %otherfeatures;
my $thatfeature = $features[1];                                                         #For the top BLAST hit (in this case)

#@features = $db->get_features_by_location($ID);

foreach my $thatfeature(@otherfeatures){

    my $name = $thatfeature->name;
    my $source = $thatfeature->source;
    my $start = $thatfeature->start;
    my $end = $thatfeature->end;
        #print "Candidate_domain_name: ", $thatfeature->name, "\t", "Candidate_Source: ", $thatfeature->source, "\t", "Candidate_domain_start: ", $thatfeature->start, "\t", "Candidate_domain_end: ", $thatfeature->end, "\n";
}


my @sorted_features = sort {$a->start <=> $b->start || $a->end <=> $b->end || $a->name cmp $b->name} @features;             #Sort first by the start, then the end, and then order it also by name
#my @sorted_features2 = sort {$a->end <=> $b->end} @features;
    #my @sorted_features2 = sort $thisfeature->end(@features);
my @namefeature = "";
foreach (@sorted_features){
    push(@namefeature, $_->name);                                                                                           #store all the names into the @namefeature array
    print "Query ", $_->name, "\t", $_->start, "\t", $_->end, "\n";
}
my @sorted_otherfeatures = sort {$a->start <=> $b->start || $a->end <=> $b->end || $a->name cmp $b->name} @otherfeatures;   #Sort first by the start, then the end, and then order it also by name
#my @sorted_otherfeatures2 = sort {$a->end <=> $b->end} @otherfeatures;
my @nameotherfeature = "";
foreach (@sorted_otherfeatures){
    push(@nameotherfeature, $_->name);                                                                                       #store all the names into the @nameotherfeature array
    print "Candidate ", $_->name, "\t", $_->start, "\t", $_->end, "\n";
#        print "Candidate ", (join "\n", @sorted_otherfeatures), "\n";
}

my $stringnames = join (',', @namefeature);                                                                                  # connect all the elements (names) of @namefeature into a single string    
my $stringnamesotherfeature = join (',', @nameotherfeature);                                                                 # connect all the elements (names) of @nameotherfeature into a single string
if (scalar @namefeature == scalar @nameotherfeature) {                                                                       #compare the number of domains
    print "Same numbers of domains", "\n";

    if ($stringnames eq $stringnamesotherfeature) {                                                                          #check the order of the domains 
        print "Same domain order along the protein when compared (order is the same in both)", "\n";
    }else{
        print "Different domain order along the protein when compared (order is not the same)", "\n";
    }
}

else {                                                                                          
    print "Not the samber number of domains", "\t", "Query domains: ", scalar @namefeature, "\t", "Candidate domains: ", scalar @nameotherfeature, "\n";
    if (scalar @namefeature < scalar @nameotherfeature) {                                                                   #if the query has less domains than the candidate
        (my $subname = $stringnames) =~ s/,/.*/g;                                                                           #get rid of the commas in the name list and substitute it with *
        if ($subname eq $stringnamesotherfeature) {                                                                         
            print "Query domain set is a subset of candidate domain set", "\n";                                             #The query domains are also present in the candidate protein, even if the query protein has less domains
        }else {
            print "Query domain set is not a subset of candidate domain set (some candidate domains lacking in the query protein set) ", "\n";
        }
    }
    if (scalar @namefeature > scalar @nameotherfeature) {                                                                   #if the candidate has less domains than the query
        (my $subname = $stringnamesotherfeature) =~ s/,/.*/g;
        if ($stringnames eq $subname) {    
            print "Candidate domain set is a subset of query domain set", "\n";                                             #The candidate domains are also present in the query protein, even if the candidate protein has less domains
        }else{
            print "Candidate domain set is not a subset of query domain set (some query domains lacking in the candidate domain set)", "\n";    #Some query domains are not present in the candidate domain set
        }
    }   
}