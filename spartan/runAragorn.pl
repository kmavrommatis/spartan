#!/jgi/tools/bin/perl
use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/../..";
use SetEnvironment;
SetEnvironment::setEnv();
use Aragorn;

my ($inFile, $outFile) = @ARGV;


my $topology = "circular";
#my $topology = "linear";
#my $intronSize;
my $intronSize = 1000;
#my $intronSize = "50,500";
my $transTable;
#my $transTable = 11;


print "Predicting RNAs via Aragorn...\n";
my $start = `date`;
chomp($start);
my $aragornObj = Aragorn->new(	"-inputFile"  => $inFile,
								"-outputFile" => $outFile
							);

$aragornObj->setSequenceTopology( $topology ) if defined($topology);
$aragornObj->setIntronSize( $intronSize ) if defined($intronSize);
$aragornObj->setTranslationTable( $transTable ) if defined($transTable);
$aragornObj->runPrediction();
my $end = `date`;
chomp($end);

print "Started at $start\nEnded at $end\n\nAdios!\n\n";