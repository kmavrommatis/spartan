package findrRNAs;

use strict;
use warnings;
use Log::Log4perl;
use FindBin;
use lib "$FindBin::Bin";
use rRNA;
use CommonFunc;
use geneCallerVersions;


sub new{
	my $self={};
        my ($class,@arguments)=@_;
        bless($self,$class);
        $self->{ logger } =  Log::Log4perl->get_logger('findrRNAm');
        for (my $i=0;$i<scalar(@arguments);$i++){
        	if(	$arguments[$i] eq '-inputFile' or
        		$arguments[$i] eq '-outputFile' or
        		$arguments[$i] eq '-domain' or
        		$arguments[$i] eq '-ssu' or 
        		$arguments[$i] eq '-lsu' or
        		$arguments[$i] eq '-tsu'){
        			my $a= substr( $arguments[$i] , 1 );
        			$self->{ $a  }=$arguments[++$i];
        			$self->{logger}->debug("Argument $a : $arguments[$i]");
        		}
			
          	if($arguments[$i] eq '-metagenome'){$self->metagenome();}
        	if($arguments[$i] eq '-isolate'){$self->isolate();}
        } 
        
        # make sure that mandatory information is given
       	if(!defined($self->{ inputFile })){$self->{logger}->logdie( "findrRNAs: No input file was provided");} 
       	if(!defined($self->{ outputFile})){$self->{outputFile} = $self->{inputFile}. ".rRNA";}
        $self->getVersion();
        $self->{partial}=0;
        return $self;
}


sub getVersion{
	my ($self)=@_;
	if( !defined( $self->{version} ) ){
		$self->{version} = geneCallerVersions::getVersion("hmmer");
	}
	
	return $self->{version};
}


# If set to metagenome domains turn to B,A,E
# Not actually doing anything at the moment
# More for future-proofness
sub partial{
	my ($self)=@_;
	$self->{partial}=1;
}



sub runPrediction{
	my ($self,  $trnaScanFlag, $trnaBlastFlag) = @_;
	
	my $dataArray = [];
	my @domains=split(",", $self->{domain});
	my $rrnaCounter = 0;
	my $trnaCounter = 0;
	
	# delete outfile if exists
	unlink($self->{outputFile}) if -f $self->{outputFile};
			
	# search
    my %dataArraysPartialGenes = ('tsu'=>[], 'ssu'=>[], 'lsu'=>[]);

	# depending on the type of gene to process call the appropriate function
	if(defined($self->{ssu}) and $self->{ssu}==1){
		$self->{logger}->debug("Predicting 16S rRNA genes");
		rRNA::findRRNA($self->{inputFile}, \@domains, 'ssu', $dataArray, $self->{partial});
	}
	if(defined($self->{tsu}) and $self->{tsu}==1){
		$self->{logger}->debug("Predicting 5S rRNA genes");
		rRNA::findRRNA($self->{inputFile}, \@domains, 'tsu', $dataArray, $self->{partial});
	}
	if(defined($self->{lsu}) and $self->{lsu}==1){
		$self->{logger}->debug("Predicting 23S rRNA genes");
		rRNA::findRRNA($self->{inputFile}, \@domains, 'lsu', $dataArray, $self->{partial});
	}
	

	# Make sure that there are no overlapping genes
	my $removedID = {};
	$self->cleanRedundant($dataArray, $removedID);
	my $counter = 0;
	my $newData = [];
	foreach my $d( @$dataArray ) {
		my $dHash = &parseAnnotation($d->[8]);
		if( $removedID->{ $dHash->{ID} } or
		($dHash->{Parent} and $removedID->{ $dHash->{Parent} }) ) {
			next;
		} 
		push @$newData, $d;
	}
	
	@$dataArray = (); # Empty old array
	
    # Write to output file
    $self->storeGFF( $newData);
}



# Resolve problems with overlapping hits
sub cleanRedundant{
	my ($self,$dataArray, $removedID) = @_;

	# Sort hit array by seq, start, length.
	@$dataArray = sort{ $a->[0] cmp $b->[0] || CommonFunc::cleanCoord($a->[3]) <=> CommonFunc::cleanCoord($b->[3]) || 
						abs(CommonFunc::cleanCoord($b->[3])-CommonFunc::cleanCoord($b->[4])) <=> 
						abs(CommonFunc::cleanCoord($a->[3])-CommonFunc::cleanCoord($a->[4])) } @$dataArray;
	
	# Score cutoffs for the different method.
	# Besides tRNAscan all values are per nucleotide.
	my $scoreCutoffs = {};
	$scoreCutoffs->{'5S'} = 0.4;
	$scoreCutoffs->{'16S'} = 0.85;
	$scoreCutoffs->{'18S'} = 0.85;
	$scoreCutoffs->{'23S'} = 0.8;
	$scoreCutoffs->{'28S'} = 0.8;
	
	# each array has the following data
	#   contig, method,"tRNA/rRNA/exon", start, end, score, strand, ".", groupData
	
	for (my $i=0; $i<scalar(@$dataArray); $i++) {
		my $query = $dataArray->[$i];
		$self->{logger}->trace("The array contains $i: ",join("\t", @$query));
		
	}
	for (my $i=0; $i<scalar(@$dataArray); $i++) {
		my $query = $dataArray->[$i];
		my $queryHash = &parseAnnotation($query->[8]);
		if ($query->[2] eq 'exon' ) { next; }
		if ($removedID->{ $queryHash->{ID} }) { next; }
		$self->{logger}->trace("Comparing $i: ",join("\t", @$query));
		for(my $j=$i-5;$j<=$i+5 ;$j++){
			next if $j < 0;
			next if $j >= scalar(@$dataArray);
			next if $i == $j;
			
			my $hit = $dataArray->[$j];
			if ($hit->[2] eq 'exon') { next; }
			
			my $hitHash = &parseAnnotation($hit->[8]);
			if ($query->[0] ne $hit->[0]) { next; }
			if ($removedID->{ $hitHash->{ID} }) { next; }
			my $overlap = CommonFunc::isOverlapping(
					CommonFunc::cleanCoord($query->[3]), 
					CommonFunc::cleanCoord($query->[4]), 
					CommonFunc::cleanCoord($hit->[3]), 
					CommonFunc::cleanCoord($hit->[4]));
			my $overlapFlag = 'no';
			my $queryGeneLength = abs(CommonFunc::cleanCoord($query->[4]) - CommonFunc::cleanCoord($query->[3])) + 1;
			my $hitGeneLength = abs(CommonFunc::cleanCoord($hit->[4]) - CommonFunc::cleanCoord($hit->[3])) + 1;
			if($overlap > 30 or $overlap > 0.2 * $queryGeneLength or $overlap > 0.2 * $hitGeneLength) {
				$overlapFlag = 'yes';
			}
			$self->{logger}->trace("       to $j: ",join("\t", @$hit));
			
			# Get the score for each gene
			my $queryScore = $query->[5] / $queryGeneLength;
			my $hitScore = $hit->[5] / $hitGeneLength;
			
			# Mark low score rRNAs
			
			my ($qType) = split " ", $queryHash->{Type};
			my ($hType) = split " ", $hitHash->{Type};
			if( $queryScore < $scoreCutoffs->{$qType} ) { 
				$dataArray->[$i]->[8] .= "; LowScore=$queryScore" if $dataArray->[$i]->[8] !~/LowScore/;
			}
			if( $hitScore < $scoreCutoffs->{$hType} ) { 
				$dataArray->[$j]->[8] .= "; LowScore=$hitScore" if $dataArray->[$i]->[8] !~/LowScore/;
			}
			
			if($overlapFlag eq 'no') { next; }
			
			# Figure out which gene to keep
			#different methods
			if($query->[1] ne $hit->[1]){
				if( $query->[1] =~/INFERNAL/ ){
					$removedID->{ $hitHash->{ID} } = 1;
				}elsif( $hit->[1] =~/INFERNAL/ ){
					$removedID->{ $queryHash->{ID} } = 1;
				}
			}
			else{# Same method
				$self->{logger}->trace( "cleanRedundant: we are comparing results of the same method\n" );
				$self->{logger}->trace("query score is $queryScore and hit score is $hitScore");
				if ( $queryScore > $hitScore ) {
					$removedID->{ $hitHash->{ID} } = 1;
					$self->{logger}->trace("We are removing $hitHash->{ID} due to score comparison");
				}
				elsif ( $queryScore == $hitScore ) {
					# Take longer one
					if ($queryGeneLength >= $hitGeneLength) {
						$removedID->{ $hitHash->{ID} } = 1;
						$self->{logger}->trace("We are removing $hitHash->{ID} due to size");
					}
					else {
						$removedID->{ $queryHash->{ID} } = 1;
						$self->{logger}->trace("We are removing $queryHash->{ID} due to size");
					}
				}
				else {
					$removedID->{ $queryHash->{ID} } = 1;
					$self->{logger}->trace("We are removing $queryHash->{ID} due to score ");
				}
			}
		}
	}
}


# Parse the annotation field
sub parseAnnotation {
	my ($string) = @_;
	my $returnHash = {};
	my @annot = split(";", $string);
	foreach my $a (@annot) {
		my ($tag, $value) = split("=", $a);
		$tag =~ s/^\s+|\s+$//g; # Remove leading and trailing whitespaces
		if (defined($value)) { 
			$value =~ s/\"//g;
			$value =~ s/^\s+|\s+$//g;
		}
		else { $value = ""; }
		$returnHash->{$tag} = $value;
	}
	return $returnHash;
}


# Store the output in the gff file
sub storeGFF{
	my ($self, $dataArray) = @_;
	
	open (my $wfh, ">".$self->{outputFile}) or $self->{logger}->logdie( "Unable to open $self->{outputFile}: $!");
	
	my $geneCount = 0;
	foreach my $d( @$dataArray ) {
		print $wfh join("\t", @$d)."\n" ;
		if ($d->[2] ne "exon") { $geneCount++ };
	}
	close $wfh;
	
	$self->{logger}->info(  "After consolidation of overlapping hits $geneCount RNA genes were stored in the file");
}
1;
