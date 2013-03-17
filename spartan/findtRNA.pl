#!/usr/bin/env perl
use strict;
use warnings;

use Log::Log4perl;
use FindBin;
use lib "$FindBin::Bin";
use SetEnvironment;
SetEnvironment::setEnv();
use findtRNAs;
use Getopt::Long;

$Getopt::Long::ignorecase = 0;

my ($inFasta,$outData,$domain,$prType, $logFile, $loglevel,$help);
my($notrnascan,$noinfernal,$yesblast);
$loglevel="INFO";
GetOptions(     
	"input:s"          		=> \$inFasta,
	"gffout:s"        		=> \$outData,
	"domain:s"        		=> \$domain,
	"genometype:s"    		=> \$prType,
	"notrnascan!"			=> \$notrnascan,
	"noinfernal!"			=>\$noinfernal,
	"blast!"				=>\$yesblast,
	"loglevel:s"			=> \$loglevel,
	"log:s"					=>\$logFile,
	"help!"=>\$help
);

if(defined($help)){printUsage();exit(0);}
if(!defined($domain)){$domain="BAE";}
if(!defined($inFasta)){
	printUsage();exit(1);
}
if(!defined($logFile)){$logFile="findrRNA.log";}
my $logger=setupLog();
# Check for tRNAs
my $tRnaFn=$outData;


my $domainCopy = $domain;
if (lc($prType) eq "metagenome" or $domain eq "BAE" ) { $domainCopy = "B,A,E"; }
if (lc($domainCopy) eq "ba") { $domainCopy = "B,A"; }
require findtRNAs;
printInfo();
my $tRNAObj = findtRNAs->new(	"-inputFile"  => $inFasta,
								"-outputFile" => $tRnaFn,
								"-domain"     => $domainCopy
							);
if (defined($prType) and lc($prType) eq 'metagenome') { $tRNAObj->metagenome(); }
else { $tRNAObj->isolate(); }


my($trnascanFlag,$blastFlag,$infernalFlag)=("on","off","on");
if(defined($noinfernal)){$infernalFlag="off";}
if(defined($notrnascan)){$trnascanFlag="off";}
if(defined($yesblast)){$blastFlag="on";}
$tRNAObj->runPrediction($trnascanFlag, $blastFlag,$infernalFlag,$infernalFlag); # No BLAST (only tRNAscan and infernal)


sub printInfo{
	$logger->info("findtRNA.pl v". toolversion());
	$logger->info( "Input      : $inFasta");
	$logger->info( "Output(gff): $outData");
	$logger->info( "Domain     : $domainCopy");

}


sub setupLog{
	my $logConf=qq{
	    log4perl.rootLogger          = $loglevel, Logfile, Screen
	    log4perl.appender.Logfile          = Log::Log4perl::Appender::File
	    log4perl.appender.Logfile.filename = $logFile
	    log4perl.appender.Logfile.layout   = Log::Log4perl::Layout::PatternLayout
	    log4perl.appender.Logfile.layout.ConversionPattern = [%p : %c - %d] - %m{chomp}%n
	    log4perl.appender.Screen         = Log::Log4perl::Appender::Screen
	    log4perl.appender.Screen.stderr  = 0
	    log4perl.appender.Screen.layout = Log::Log4perl::Layout::PatternLayout
	    log4perl.appender.Screen.layout.ConversionPattern = [%p - %d] - %m{chomp}%n
	};
	
	
	Log::Log4perl->init(\$logConf);
	my $logger=Log::Log4perl->get_logger("findrRNA");

	return $logger;
}


sub printUsage{
	print 
	q{
Program: findtRNA.pl by Konstantinos Mavrommatis (KMavrommatis@lbl.gov)
Use this program to predict tRNA genes in bacterial, archaeal and eukaryotic genomes
Arguments:
	-input <input file>. Input file is in fasta format
	-gffout <output file>. The output file will be in gff format
	-domain <domain of organism>. The domain can be any combination of B, A, E for Bacteria, Archaea, Eukaryota respectively. Default BAE i.e. all
	-genometype <type of sequence>. The genome type can be any of metagenome, isolate. Default isolate.
	-log <logfile>
	-loglevel <[INFO], DEBUB, TRACE>
};
}



sub toolversion{
	my $version = '1.0';
	
	return $version;
}
