#!/jgi/tools/bin/perl
use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/../..";



use SetEnvironment;
SetEnvironment::setEnv();
use Infernal;
use File::Basename;

my ($inFile,$outFile) = @ARGV;
my $wfh;
my $workDir=$ENV{ SCRATCH };
if(!$inFile or !-e $inFile){die "Please provide a valid input fasta file\n";}
if(defined($outFile) ){ 
	open $wfh, ">".$outFile;
	$workDir = dirname($outFile)."/";
}
else{$wfh = \*STDOUT;}


#print "Predicting RNAs via Rfam...\n";
my $start = `date`;
chomp($start);
my $rRFAMObj = Infernal->new(	"-inputFile"  => $inFile,
								"-outputFh" => $wfh,
								"-path"       => $workDir
							);
$rRFAMObj->excludeModels( "db" );
#$rRFAMObj->constrainModelsTo( "trnas" );
$rRFAMObj->runPrediction();
my $end = `date`;
chomp($end);
close ($wfh) if defined($outFile);

#print "Started at $start\nEnded at $end\n\nAdios!\n\n";
