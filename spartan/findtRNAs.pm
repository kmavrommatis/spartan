package findtRNAs;

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin";
use tRNA;
use CommonFunc;
use geneCallerVersions;
use Log::Log4perl;



sub new{
	my $self={};
        my ($class,@arguments)=@_;
        bless($self,$class);
        $self->{logger} =  Log::Log4perl->get_logger('tRNAm');
        for (my $i=0;$i<scalar(@arguments);$i++){
        	if(	$arguments[$i] eq '-inputFile' or
        		$arguments[$i] eq '-outputFile' or
        		$arguments[$i] eq '-domain'){
        			my $a= substr( $arguments[$i] , 1 );
        			$self->{ $a  }=$arguments[++$i];
        		}

          	if($arguments[$i] eq '-metagenome'){$self->metagenome();}
        	if($arguments[$i] eq '-isolate'){$self->isolate();}
        } 
        
        # make sure that mandatory information is given
       	if(!defined($self->{ inputFile })){$self->{logger}->logdie( "findtRNAs: No input file was provided\n");} 
       	if(!defined($self->{ outputFile})){$self->{outputFile} = $self->{inputFile}. ".tRNA";}
        $self->getVersion();
        return $self;
}


sub getVersion{
	my ($self)=@_;
	if( !defined( $self->{trnaBlastVersion} ) ) {
		$self->{trnaBlastVersion} = geneCallerVersions::getVersion("trnablast");
		$self->{methodOrder}->{$self->{trnaBlastVersion}  } =15; # we trust tBLAST more than anything
	}
	if( !defined( $self->{tRNAscanVersion} ) ) {
		$self->{tRNAscanVersion} = geneCallerVersions::getVersion("trnascan");
		$self->{methodOrder}->{$self->{tRNAscanVersion}  } =12; # we trust tRNA scan
	}
	if( !defined( $self->{aragornVersion} ) ) {
		$self->{aragornVersion} = geneCallerVersions::getVersion("aragorn");
		$self->{methodOrder}->{$self->{aragornVersion}  } =10; # aragorn is better than infernal (because it finds introns)
	}
	if( !defined( $self->{cmsearchVersion} ) ) {
		$self->{cmsearchVersion} = geneCallerVersions::getVersion("cmsearch11");;
		$self->{methodOrder}->{$self->{cmsearchVersion}  } =5; # infernal is the last in order.
	}
	
	
	
	return "tRNAscan: ". $self->{tRNAscanVersion} .
			", BLAST: ". $self->{trnaBlastVersion} .
			", Aragorn: ". $self->{aragornVersion}.
			", infernal: ". $self->{cmsearchVersion};
			
	
}


# If set to metagenome domains turn to B,A,E
# Not actually doing anything at the moment
# More for future-proofness
sub metagenome{
	my ($self)=@_;
	$self->{isolate}=undef;
	$self->{metagenome}=1;
}

sub isolate{
	my ($self)=@_;
	$self->{isolate}=1;
	$self->{metagenome}=undef;
}


sub runPrediction{
	my ($self, $trnaScanFlag, $trnaBlastFlag, 
					$aragornFlag, $infernalFlag,$transTable, $circular) = @_;
	
	my $dataArray = [];
	my $trnaCounter = 0;
	
	# delete outfile if exists
	unlink($self->{outputFile}) if -f $self->{outputFile};
			
	# search
	my @domains = split(",", $self->{domain});

	$self->{logger}->info("Checking for tRNA genes.");
	tRNA::findTRNA($self->{inputFile}, \@domains, $dataArray, $trnaScanFlag, 
						$trnaBlastFlag, $aragornFlag, $infernalFlag, $transTable, $circular);
	
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
	
	$dataArray = []; # Empty old array
	
    # Write to output file
    $self->storeGFF($newData);
}


# Resolve problems with overlapping hits
# the algorithm is as follows:

# all trna predictions are kept and considered correct
# ARAGORN and INFERNAL that overlaps trnaScan are removed

# if a prediction is done by both INFERNAL and ARAGORN
# we keep the ARAGORN prediction. This ensures a correct
# detection of introns

# if a prediction is done by INFERNAL only is kept

# if a prediction is done by ARAGORN only it is discarged
# since tests have shown that ARAGORN predicts spurious genes.


sub cleanRedundant{
	my ($self, $dataArray, $removedID) = @_;

	# Sort hit array by seq, start, length.
	@$dataArray = sort{ $a->[0] cmp $b->[0] || $a->[3] <=> $b->[3] || 
						abs($b->[3]-$b->[4]) <=> abs($a->[3]-$a->[4]) } @$dataArray;
	
	# Score cutoffs for the different method.
	# Besides tRNAscan all values are per nucleotide.
	$self->{scoreCutoffs} = {};
	$self->{scoreCutoffs}->{ $self->{trnaBlastVersion} } = 0.8;
	$self->{scoreCutoffs}->{ $self->{tRNAscanVersion} } = 20;
	$self->{scoreCutoffs}->{ $self->{cmsearchVersion} } = 15;
	$self->{scoreCutoffs}->{ $self->{aragornVersion} } = 50;
	# each array has the following data
	#   contig, method,"tRNA/rRNA/exon", start, end, score, strand, ".", groupData
	
	for (my $i=0; $i<scalar(@$dataArray); $i++) {
		my $query = $dataArray->[$i];
		my $queryHash = &parseAnnotation($query->[8]);
		if ($query->[2] eq 'exon' ) { next; }
		if ($removedID->{ $queryHash->{ID} }) { next; }
		$self->{logger}->trace("Comparing $i ", join("\t",@{$query}[0..7] ) );
		for(my $j=$i-7;$j<=$i+7 ;$j++) {
			next if $j < 0;
			next if $j >= scalar(@$dataArray);
			next if $i == $j;
			
	 		
			my $hit = $dataArray->[$j];

			
			if ($hit->[2] eq 'exon') { next; }
			
			$self->{logger}->trace("       to $j ", join("\t",@{$hit}[0..7]   ) );

			my $hitHash = &parseAnnotation($hit->[8]);
			if ($query->[0] ne $hit->[0]) { 
				$self->{logger}->trace("These are coming from $query->[0] and $hit->[0] and will be ignored");
				next; 
			}
			if ($removedID->{ $hitHash->{ID} }) { 
				$self->{logger}->trace("This hit ($hitHash->{ID}) is already marked to be removed and will be ignored ");
				next; }
			
			my $overlap = CommonFunc::isOverlapping($query->[3], $query->[4], $hit->[3], $hit->[4]);
			my $overlapFlag = 'no';
			my $queryGeneLength = abs($query->[4] - $query->[3]) + 1;
			my $hitGeneLength = abs($hit->[4] - $hit->[3]) + 1;
			if($overlap > 30 or $overlap > 0.2 * $queryGeneLength or $overlap > 0.2 * $hitGeneLength) {
				$overlapFlag = 'yes';
				$self->{logger}->trace("The two genes are overlapping by $overlap");
			}else{
				$self->{logger}->trace("The overlap between the two genes is minimal ($overlap)");
			}
			
			# Get the score for each gene
			my $queryScore=$self->getGeneScore($query);
			my $hitScore=$self->getGeneScore($hit);

			# Mark tRNAs with low scores
			$self->markLowScore( $query, $queryScore, $dataArray->[$i] );
			$self->markLowScore( $hit, $hitScore, $dataArray->[$j] );

#	 		print "Between the two there is an overlap of $overlap. The flag is set to $overlapFlag\n";
			if($overlapFlag eq 'no') { next; }
					
				
			# Figure out which gene to keep

			# same method
			if ($query->[1] eq $hit->[1]) {
				# keep the prediction with the better score
				# or if scores are the same the longest gene
				if($query->[5] eq '.' ){$query->[5]=1;}
				if($hit->[5] eq '.' ){$hit->[5]=1;}
				$self->{logger}->trace("The two predictions come from the same method");
				if ($query->[5] > $hit->[5]) {
					$removedID->{ $hitHash->{ID} } = 1;
				}
				elsif ($query->[5] == $hit->[5]) {
					# Take longer one
					if ($queryGeneLength >= $hitGeneLength) {
						$removedID->{ $hitHash->{ID} } = 1;
					}
					else {
						$removedID->{ $queryHash->{ID} } = 1;
					}
				}
				else {
					$removedID->{ $queryHash->{ID} } = 1;
				}
			}
			else {
				# Different methods
#				print "Comparing different methods\n";
				$self->{logger}->trace("The two predictions come from different method");
				my $queryType;
				my $hitType;
				my @splitA;
				$queryType = $query->[1];
				$hitType = $hit->[1];
				
				if ($queryScore >= $self->{scoreCutoffs}->{ $queryType }) {
					if ($hitScore < $self->{scoreCutoffs}->{ $hitType }) {
						$self->{logger}->trace( "Comparing $hit->[1] and $query->[1] and I drop the hit due to score ($hit->[1])");
						$removedID->{ $hitHash->{ID} } = 1;
					}
					# if both scores are above the thersholds
					else {
						if($self->{methodOrder}->{  $hit->[1] } >
						   $self->{methodOrder}->{ $query->[1]}){
						   	
						   	$self->{logger}->trace( "Comparing $hit->[1] ($hit->[3] .. $hit->[4]) and $query->[1] and I drop the query ($query->[1])");
						   	
						   	$removedID->{ $queryHash->{ID} } = 1;
						   	# since aragorn predictions are removed at the end
						   	# replace aragorn with cmsearch
						   	if($query->[1] eq $self->{ cmsearchVersion }){
						   		$hit->[8].="; Supported=".$self->{ cmsearchVersion };
						   		
						   	}
						}
						elsif($self->{methodOrder}->{  $query->[1] } >
						   $self->{methodOrder}->{ $hit->[1]}){
						   	
						   		$self->{logger}->trace( "Comparing $hit->[1] ($hit->[3] .. $hit->[4]) and $query->[1] and I drop the hit ($hit->[1])");
						   	
							$removedID->{ $hitHash->{ID} } = 1;	
							# since aragorn predictions are removed at the end
						   	# replace aragorn with cmsearch
						   	if($hit->[1] eq $self->{ cmsearchVersion }){
						   		$query->[8].="; Supported=".$self->{ cmsearchVersion };
						   	}
						}

						else {
							print "DataArray:\n";
							foreach my $d (@$dataArray) {
								print join("\t", @$d)."\n";
							}
							print "\n";

							$self->{logger}->logdie( "\n\nA ".$query->[2]." (method: ".$query->[1].
								", score: $queryScore) and a ".
								$hit->[2]."(method: ".$hit->[1].
								", score: $hitScore) are overlapping on the ".
								"sequence ".$query->[0]." and both have good scores!\n\n");
						}
					}
				}
				else {
					if ($hitScore >= $self->{scoreCutoffs}->{ $hitType }) {
						$self->{logger}->trace( "Comparing $hit->[1] and $query->[1] and I drop the query due to score ($query->[1])");
						$removedID->{ $queryHash->{ID} } = 1;
					}
					else {
						if($self->{methodOrder}->{  $hit->[1] } >
						   $self->{methodOrder}->{ $query->[1]}){
						   	$self->{logger}->trace( "Comparing $hit->[1] ($hit->[3] .. $hit->[4]) and $query->[1] and I drop the query ($query->[1])");
						   	$removedID->{ $queryHash->{ID} } = 1;
						   	# since aragorn predictions are removed at the end
						   	# replace aragorn with cmsearch
						   	if($query->[1] eq $self->{ cmsearchVersion }){
						   		$hit->[8].="; Supported=".$self->{ cmsearchVersion };
						   	}
						}
						elsif($self->{methodOrder}->{  $query->[1] } >
						   $self->{methodOrder}->{ $hit->[1]}){
						   	
						   	$self->{logger}->trace( "Comparing $hit->[1] ($hit->[3] .. $hit->[4]) and $query->[1] and I drop the hit ($hit->[1])");
						   	
						   	$removedID->{ $hitHash->{ID} } = 1;
						   	# since aragorn predictions are removed at the end
						   	# replace aragorn with cmsearch
						   	if($hit->[1] eq $self->{ cmsearchVersion }){
						   		$query->[8].="; Supported=".$self->{ cmsearchVersion };
						   	}
						}

						else {
							my $str= "DataArray:\n";
							foreach my $d (@$dataArray) {
								$str.= join("\t", @$d)."\n";
							}
							$self->{logger}->info($str);

							$self->{logger}->logdie( "\n\nA ".$query->[2]." (method: ".$query->[1].
								", score: $queryScore) and a ".
								$hit->[2]."(method: ".$hit->[1].
								", score: $hitScore) are overlapping on the ".
								"sequence ".$query->[0]." and both have bad scores!\n\n");
						}
					}
				}	
			}
		}
	}
	
	
	#remove the genes that are predicted by aragorn and have not been supported
	# by other methods
	for (my $i=0; $i<scalar(@$dataArray); $i++) {
		my $query = $dataArray->[$i];
		my $queryHash = &parseAnnotation($query->[8]);
		if ($query->[2] eq 'exon' ) { next; }
		if ($removedID->{ $queryHash->{ID} }) { next; }
		$self->{logger}->trace("Proceeding with line $i ", join("\t", @$query));
		if($query->[1] eq $self->{aragornVersion} and !defined($queryHash->{Supported})){
				$removedID->{ $queryHash->{ID} } = 1;
				$self->{logger}->trace("this entry is not supported thus it will be removed");
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
	$self->{logger}->debug( "Storing data in ", $self->{outputFile} ."\n");
	open (my $wfh, ">".$self->{outputFile}) or 
					die "Unable to open ".$self->{outputFile}.": $!\n";
	#print "Storing ".scalar(@$dataArray)." RNAs in the file...\n";
	
	my $geneCount = 0;
	foreach my $d( @$dataArray ) {
		print $wfh join("\t", @$d)."\n" ;
		if ($d->[2] ne "exon") { $geneCount++ };
	}
	close $wfh;
	
	$self->{logger}->info(  "\n$geneCount RNA genes were stored in the file.\n");
}

# get the score for this gene
# if this is a blast prediction we find the bitscore/nucleotide
# for other methods we use the score that the method returns
sub getGeneScore{
	my($self, $query)=@_;
	my $queryScore;
	my $queryGeneLength = abs($query->[4] - $query->[3]) + 1;
	if ($query->[1] eq $self->{trnaBlastVersion}) {
		$queryScore = $query->[5] / $queryGeneLength;
	}
	else {
		$queryScore = $query->[5];
	}
	if($queryScore eq "."){$queryScore=50;} # aragorn does not have a numeric value for score
	return $queryScore;
}


# mark genes with low scores by adding the line LowScore=score
# in the annotation line
sub markLowScore{
	my($self,$query,$queryScore, $darray)=@_;
	if ($query->[1] eq $self->{trnaBlastVersion} and 
		$queryScore < $self->{scoreCutoffs}->{ $self->{trnaBlastVersion} }) {
		$darray->[8] .= "; LowScore=$queryScore" if $darray->[8]!~/LowScore=$queryScore/;
	}
	elsif ($query->[1] eq $self->{tRNAscanVersion} and 
			$queryScore < $self->{scoreCutoffs}->{ $self->{tRNAscanVersion} }) {
		$darray->[8] .= "; LowScore=$queryScore" if $darray->[8]!~/LowScore=$queryScore/;
	}
}
	
1;
