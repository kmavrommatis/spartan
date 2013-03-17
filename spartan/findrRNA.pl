#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use Log::Log4perl;
use lib "$FindBin::Bin";



use SetEnvironment;
SetEnvironment::setEnv();
use findrRNAs;
#use createGenomeFile;
#use File::Basename;
use Getopt::Long;

$Getopt::Long::ignorecase = 0;

my ($inFasta,$outData, $help,$rnaType, $domain,$logFile,$partial);

my $loglevel="INFO";
GetOptions(     
	"input:s"          		=> \$inFasta,
	"gffout:s"        		=> \$outData,
	"rnatype:s"			=> \$rnaType,
	"domain:s"        		=> \$domain,
	"partial!"    		=> \$partial,
	"loglevel:s"			=> \$loglevel,
	"log:s"					=>\$logFile,
	"help!"=>\$help
);


if(defined($help)){
	printUsage();
	exit(0);
}

if(!defined($logFile)){$logFile="findrRNA.log";}
my $logger=setupLog();
if(!defined($domain)){$domain="BAE";}
if(!defined($rnaType)){$rnaType='tsu,ssu,lsu';}
if(!defined($inFasta)){
	printUsage();
	exit(1);
}
my %h=('tsu'=>0,'ssu'=>0,'lsu'=>0);
if(lc($rnaType) =~/tsu/){ $h{tsu}=1;}
if(lc($rnaType) =~/ssu/){ $h{ssu}=1;}
if(lc($rnaType) =~/lsu/){ $h{lsu}=1;}
# Check for tRNAs

my @t=split(//, uc($domain));
my $domainCopy = join(",", @t);
if ( $domain eq "BAE" ) { $domainCopy = "B,A,E"; }
if (lc($domainCopy) eq "ba") { $domainCopy = "B,A"; }
require findrRNAs;


printInfo();



my $rRNAObj = findrRNAs->new(	"-inputFile"  => $inFasta,
								"-outputFile" => $outData,
								"-domain"     => $domainCopy,
								"-tsu" => $h{tsu},
								"-ssu" => $h{ssu},
								"-lsu" => $h{lsu}
							);
if (defined($partial)) { $rRNAObj->partial(); }
$rRNAObj->runPrediction();



sub printInfo{
	$logger->info("findrRNA.pl v". toolversion());
	$logger->info( "Input      : $inFasta");
	$logger->info( "Output(gff): $outData");
	$logger->info( "Domain     : $domainCopy");
	$logger->info( "5S         : ".  $h{ tsu });
	$logger->info( "16S        : ".  $h{ ssu });
	$logger->info( "23S        : ".  $h{ lsu });
	$logger->info( "accepting partial genes") if defined($partial);
	$logger->info( "accepting only full size genes") if !defined($partial);
}

sub printUsage{
	print 
	q{
Program: findrRNA.pl by Konstantinos Mavrommatis (KMavrommatis@lbl.gov)
Use this program to predict rRNA (23S, 16S, 5S) genes in bacterial, archaeal and eukaryotic genomes
Arguments:
	-input <input file>. Input file is in fasta format
	-gffout <output file>. The output file will be in gff format
	-rnatype <RNA type to search>. The RNA type can be a comma list of tsu, ssu, lsu for 5S, 16S, 23S respectively. Default 'tsu,ssu,lsu' i.e. all
	-domain <domain of organism>. The domain can be any combination of B, A, E for Bacteria, Archaea, Eukaryota respectively. Default BAE i.e. all
	-partial. To allow for partial gene prediction, i.e. genes located at the ends of contigs. [ default is look for full size genes ].
	-log <logfile>
	-loglevel <[INFO], DEBUG, TRACE>
};
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




sub toolversion{
	my $version = '1.0';
	
	return $version;
}