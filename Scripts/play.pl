use strict;
use warnings;

use Bio::DB::SeqFeature::Store;
use Bio::SeqFeature::Generic;


 # Open the sequence database
 my $db      = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                                -dsn     => 'dbi:mysql:3endprocessing',
						-user => 'root',
						-password => '');

#my $sequence_id = 'CBX89940';
#
#my $iterator = $db->get_seq_stream;
#  while (my $feature = $iterator->next_seq) {
#    print $feature->primary_tag,' ',$feature->display_name,"\n";
#  }
my @features = $db->get_features_by_type("polypeptide");
# my @features = $db->features(-seq_id => 'CADAFLAT00007143|mRNA');
# my $seq = $db->get_Seq_by_id('CADAFLAT00007143|mRNA');
my $thisfeature = shift @features;

my $seq = $thisfeature->seq;
my $id = $thisfeature->id;
my $start = $thisfeature->start;
my $end = $thisfeature->end;

@features = $db->get_features_by_location($id, $start, $end);

my $seqobj = Bio::Seq->new(-id => $thisfeature->load_id,
			   -display_id => $thisfeature->load_id,
				-seq => $thisfeature->seq->seq);

#my @features = $seq->features;

foreach my $feature(@features){
	my $compliant_feature = Bio::SeqFeature::Generic->new(
	-primary_id => $feature->primary_id,
	-start => $feature->start,
	-end => $feature->end,
	-strand => $feature->strand,
	-display_name => $feature->display_name,
	-seq_id => $feature->seq_id,
	-type => $feature->type,
	-source_tag => $feature->source_tag,
	-primary_tag => $feature->primary_tag,
	);
	
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
	print "Type ", $feature->type, " name " , $feature->name, " start ", $feature->start, "  end ", $feature->end, " DESC: ", $term, "\n";
	$seqobj->add_SeqFeature($compliant_feature);

}

        use Tk;
        use Bio::SeqIO;
        use Bio::Tk::SeqCanvas;

        Begin();
        MainLoop;

        sub Begin {

        # set up the Tk Windows

        my $MW = MainWindow->new (-title => "Map Of BioSeq Object");
        my $Frame = $MW->Frame()->pack(-side => 'top');
        my $lblSysMess = $MW->Label()->pack(-side => 'bottom', -fill => 'both');

        # Draw the Map
        my $axis_length = 800;  # how large (long axis) I want the final map to be

        my $MapObj = Bio::Tk::SeqCanvas->new(
               $axis_length,
               $Frame,
               $lblSysMess,
               $seqobj,
               -orientation => 'horizontal',
               label => 'primary_tag',
               width => 200,
               );


	}


print ;
print ;
