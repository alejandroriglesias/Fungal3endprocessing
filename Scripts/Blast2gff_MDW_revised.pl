#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           | 
# Blast2Gff.pl                                              |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# EDITED AND DEBUGGED BY:  Mark Wilkinson, April, 2014      |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 04/17/2007                                       |
# UPDATED: 04/18/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Converts BLAST output to GFF format. This the the GFF    |
#  format that is used with the Apollo Genome Annotation    |
#  curation program.                                        |
#  Currently this only works with m8 blast output.          |
#                                                           |
# VERSION:                                                  |
# $Id::                                                  $: |
#                                                           |
#-----------------------------------------------------------+

=head1 NAME

Blast2Gff.pl - Convert BLAST output to GFF format.

=cut

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;                   # Keeps thing running clean
use Getopt::Std;              # Get options from command line
use Path::Iterator::Rule;


=head1 VARIABLES
=cut
#-----------------------------+
# VARIABLES                   |
#-----------------------------+
my $GffAppend;                 # BOOLEAN. Append to GFF file
my $InFolder;                    # Full path to the input blast output file
my $OutFile;                   # Full path to the gff formatted output file
my $AlignFormat;               # Alignment format of the blast output file
                               # ie. -m = 0,8, or 9
my $PrintHelp;                 # Boolean, print the Usage statement
my $BlastDb;                   # Blast database 
my $BlastProg;                 # Blast program (ie. blastn, blastx

my $Usage = "USAGE:\n".
    "Blast2Gff.pl -i InFile.Fasta -o OutFile.gff -d BlastDb\n".
    " -p blastprogram -m AligFormat -s SeqName -a\n\n".
    " -i Full path to the folder with subfolders containing blast outputs" .
    " -o Full path for the GFF formated output file [STRING]\n".
    " -m Format of the algnment outout from blast [INTEGER]\n".
    "    Default value is 8. Valid values are 0,8,9".
    " -p Blast program used [STRING]\n".
    "    Default is blastn\n".
    " -a Append results to the gff file [BOOLEAN]\n".
    "    Default is to overwrite any exiting file.\n".
    " -h Print this help statement [BOOLEAN]\n";

=head1 COMMAND LINE VARIABLES
=cut
#-----------------------------+
# COMMAND LINE VARIABLES      |
#-----------------------------+
my %Options;
getopts('d:i:o:m:p:s:h:a', \%Options);

$PrintHelp = $Options{h};
if ($PrintHelp)
{
    print $Usage;
    exit;
}

$GffAppend = $Options{a};
$InFolder = $Options{i} ||
    die "\n\n\n\aERROR: An input file must be specified.\n\n$Usage\n\n\n";
$OutFile = $Options{o} ||
    die "\n\n\n\aERROR: An output filename must be specified.\n\n$Usage\n\n\n";
# Default output is the full path of the input file with the gff extension
$BlastProg = $Options{p} ||
    "blastn";
$AlignFormat = $Options{m} || 
    "8";                        # Default format is tab delim

my $GFFOUT;
if ($GffAppend)
{
    open ($GFFOUT, ">>".$OutFile) ||
        die "Can not open GFF ouput file.$OutFile.\n";
} else {
    open ($GFFOUT, ">".$OutFile) ||
        die "Can not open GFF ouput file.$OutFile.\n";
}    

 
my $rule = Path::Iterator::Rule->new;
 
my $it = $rule->iter( $InFolder );
while ( my $InFile = $it->() ) {
    print "processing: $InFile\n";
    next if (-d $InFile);

    # If append was selected, just append gff data to the
    # output file
    
    &TabBlast2Gff ($InFile, $GFFOUT, $BlastDb, $BlastProg);
}

close GFFOUT;
print "\n\n\nDONE!\n\n\n";
exit 1;


#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub TabBlast2Gff 
{
    my $In = $_[0];       # Path to blast intput file
    my $Out = $_[1];      # Path to gff output file
    my $Db = $_[2];       # The BLAST databas the hits are derived from
    my $Prog = $_[3];     # BLAST program used
    #my $GStart;           # GFF Start position 
    #my $GEnd;             # GFF End position

#    my $Format = $_[4];   # Format of the blast file
#                          # 8,9, 0 etc
#    my $UseScore = $_[5]; # Score format to use
    
    my $HitNum = "0";
    #-----------------------------+
    # FILE I/O                    |
    #-----------------------------+
    open (BLASTIN, "<".$In) ||
        die "\n\nCan not open BLAST input file.$In.\n\n\n";
    

    while (<BLASTIN>)
    {
        if (substr($_, 0, 1) eq '#')
        {
            # print "$_\n";
            next;
        }

        
        $HitNum++;
   
        my ($QryID, $SubID, $PID, $Len, 
            $MisMatch, $GapOpen, 
            $QStart,$QEnd, $SStart, $SEnd,
            $EVal, $BitScore) = split(/\t/);
        
        my $SStrand;
        my $QStrand;
        my $Frame = ".";
        
        # Set the start to be less then the end
        # This info can be used to dedeuct the strand

        my ($GQstart, $GQend);
        if ($QStart < $QEnd)
        {
            $GQstart = $QStart;
            $GQend = $QEnd;
            $QStrand = "+";
        } elsif ($QStart > $QEnd) {
            $GQstart = $QStart;
            $GQend = $QEnd;
            $QStrand = "-";
        } else {
            die "Unexpected Query Start and End:\n\tS:$QStart\n\tE:$QEnd";
        }

        my ($GSstart, $GSend);
        if ($SStart < $SEnd)
        {
            $GSstart = $SStart;
            $GSend = $SEnd;
            $SStrand = "+";
        } elsif ($QStart > $QEnd) {
            $GSstart = $SStart;
            $GSend = $SEnd;
            $SStrand = "-";
        } else {
            die "Unexpected Subject Start and End:\n\tS:$SStart\n\tE:$SEnd";
        }
        
        # Trim leading white space from Bit score
        $BitScore =~ s/^\s*(.*?)\s*$/$1/;
        
        my $QGFF = 
            $QryID."\t".
            $Prog."\t".    
            $Prog."_HIT_$SubID\t".           # Feature
            $GQstart."\t".          # Start
            $GQend."\t".            # End
            $BitScore."\t".        # Score
            $QStrand."\t".          # Strand
            $Frame."\t".           # Frame
            "Name=$Prog"."_HIT_".$SubID.";ID=$QryID|$SubID|$GQstart|$GQend;Parent=$QryID" .                 # Attribute
            "\n";

        my $SGFF = 
            $SubID."\t".
            $Prog."\t".    
            $Prog."_HIT_$QryID\t".           # Feature
            $GSstart."\t".          # Start
            $GSend."\t".            # End
            $BitScore."\t".        # Score
            $SStrand."\t".          # Strand
            $Frame."\t".           # Frame
            "Name=$Prog"."_HIT_".$QryID.";ID=$SubID|$QryID|$GSstart|$GSend;Parent=$SubID" .                 # Attribute
            "\n";

        print $Out $QGFF, $SGFF;
        print $QGFF, $SGFF;
        
    } # END OF WHILE BLASTIN
    
} # END OF Blast2Gff Subfunction