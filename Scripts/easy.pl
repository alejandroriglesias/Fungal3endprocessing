use strict;
use warnings;

use Bio::DB::SeqFeature::Store;
use Bio::SeqFeature::Generic;

#my $ID = <stdin>;
my $ID = $ARGV[0];									#The ID is introduced when running the script in the terminal
 # Open the sequence database
 my $db      = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                                -dsn     => 'dbi:mysql:Fungal3endprocessing',
						-user => 'root',
						-password => '');
 

#
#my $iterator = $db->get_seq_stream;
#  while (my $feature = $iterator->next_seq) {
#    print $feature->primary_tag,' ',$feature->display_name,"\n";
#  }
#my @features = $db->get_features_by_type("polypeptide");
## my @features = $db->features(-seq_id => 'CADAFLAT00007143|mRNA');
## my $seq = $db->get_Seq_by_id('CADAFLAT00007143|mRNA');
my @fulllengthfeatures = $db->features(-seqid => $ID, -type => "polypeptide");		#Get the features of the full-length protein
my $FULLLENGTH = shift @fulllengthfeatures;
my $seqlength = $FULLLENGTH->seq->length;

my @features = $db->get_features_by_location($ID);
											#Get the rest of the features
my $thisfeature = $features[0];

#my $seq_id =$thisfeature->ref;
#my $seq = $thisfeature->seq;
#my $id = $thisfeature->id;
#my $start = $thisfeature->start;
#my $end = $thisfeature->end;



my $seqobj = Bio::Seq->new(-id => $thisfeature->load_id, 
			   -display_id => $thisfeature->load_id,
				-seq => $FULLLENGTH->seq->seq);

#my @features = $seq->features;
$ID=0;
foreach my $feature(@features){
    ++$ID;
	my $term = " ";
	if ($feature->attributes('signature_desc')){
		my @terms = $feature->attributes('signature_desc');
		$term = $terms[0];
	}
	if ($feature->type eq 'protein_match:PANTHER'){
		my $name = $feature->name;
		open (IN, "panther.txt") || die "can't open file $!\n";
		while (my $line = <IN>){
			if ($line =~ /^$name/){
				my @elements = split '\t', $line;
				$term = $elements[1];
				last;
			}
		}
	}
        #if ($feature-> type eq 'protein_match:SUPERFAMILY') {
            #my $name = $feature->name;
            #open(IN, "") or die "can't open file $!\n";
            # while (my $line = <IN>){
              #  if ($line = ~/^$name/) {
             #       my @elements = split '\t', $line;
            #        $term = $elements[1];
           #         last;
          #      }
                
         #   }
        
        if ($feature->start > '20') {
            print "Shorter start than original!","\n";
        }
        
        if ($feature->end < '500') {
            print "Shorter end than original!","\n";
        }
        
	print "ID ", $feature->display_name, " Type ", $feature->type, " name " , $feature->name, " start ", $feature->start, "  end ", $feature->end, " DESC: ", $term, "\n";
        
        
        
       
        

	#my $compliant_feature = Bio::SeqFeature::Generic->new(
	#-primary_id => $feature->primary_id,
	#-start => $feature->start,
	#-end => $feature->end,
	#-strand => $feature->strand,
	#-display_name => $feature->display_name,
	#-seq_id => $feature->seq_id,
	#-type => $feature->type.$ID,
	#-source_tag => $feature->source_tag.$ID,
	#-primary_tag => $feature->primary_tag.$ID,
	#);
	#$compliant_feature->{term} = $term;
	

	#$seqobj->add_SeqFeature($compliant_feature);

        #print "Start ", $feature->start, " End ", $feature->end, " Name ", $feature->display_name, " Seq_ID ", $feature->type.$ID, " Type ", $feature->type.$ID, " Source_tag ", $feature->source_tag.$ID, " Primary_tag ", $feature->primary_tag.$ID, "\n";
        
}

