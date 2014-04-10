use strict;
use warnings;

use Bio::DB::SeqFeature::Store;

 # Open the sequence database
 my $db      = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                                -dsn     => 'dbi:mysql:3EndProcessing',
						-user => 'uroot',
						#-password => 'LetMeIn');

my $sequence_id = 'CADAFLAT00007143|mRNA';
# my @features = $db->get_features_by_type('protein_match:Gene3D');
# my @features = $db->features(-seq_id => 'CADAFLAT00007143|mRNA');
# my $seq = $db->get_Seq_by_id('CADAFLAT00007143|mRNA');
my $seq = $db->segment($sequence_id);
my $seqobj = Bio::Seq->new( -display_id => $sequence_id,
                             -seq => $seq->seq->seq);
my @features = $seq->features;
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

