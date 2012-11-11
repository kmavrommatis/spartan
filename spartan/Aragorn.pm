package Aragorn;

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/..";
use geneCallerVersions;


# Wrapper around Aragorn.
# Runs Aragorn with provided parameters
# and creates a gff3 output file.
# Written by Marcel Huntemann.

sub new{
	my $self = {};
    my ($class, @arguments) = @_;
    bless($self, $class);
    
	for (my $i=0; $i<scalar(@arguments); $i++) {
     	if( $arguments[$i] eq '-inputFile' or
     		$arguments[$i] eq '-outputFile' ) {
     				
    		my $a = substr( $arguments[$i] , 1 );
    		$self->{ $a } = $arguments[++$i];
    	}
	}
	
	$self->{bin} = $ENV{ AragornBin };
	$self->{options} = " -d -w ";
        
    # Make sure that mandatory information is given
    if( !defined($self->{ inputFile }) ) { die "Aragorn: No input file was provided!\n"; } 
    if( !defined($self->{ outputFile}) ) { $self->{outputFile} = $self->{inputFile}.".aragorn.gff"; }
    $self->{bin} = $ENV{ AragornBin };
    if( !defined($self->{bin} ) or !-x $self->{bin} ) { 
    	die "Aragorn: Cannot locate or execute Aragorn executable!\n";
    }
 	
   	$self->getVersion();
      	
    return $self;
}


sub runPrediction {
	my ($self) = @_;
	
	my $cmd = $self->{bin}.$self->{options}.$self->{inputFile};
	
	print "Aragorn: running cmd: $cmd\n";
	
	open WH, ">", $self->{outputFile} or die "Can't open ".$self->{outputFile}.": $!\n";
	
	open CMD, "$cmd |" or die "Command failed: $cmd : $!\n";
	while (my $line = <CMD>) {
		chomp($line);
		if( $line =~ /^>end\s+/ ) { last; }
		
		if( $line =~ /^>(\S+)\s*.*/ ) { # New sequence
			my $seqName = $1;
			$line = <CMD>;
			$line =~ /^(\d+)\sgenes\sfound/;
			my $numberOfFoundGenes = $1;
			
			for( my $i=1; $i<=$numberOfFoundGenes; $i++ ) {
				$line = <CMD>;
				chomp($line);
				my @sA = split /\s+/, $line;
				
				my ($rna, $permuted, $type, $start, $end, $strand, $anticodonPos, 
						$tagStart, $tagEnd, $anticodon, $intronPos, $intronLen, 
						$tagPeptide, $fragment1End, $fragment2Start);
				
				if( $sA[2] =~ /^c/ ) { $strand = "-"; }
				else { $strand = "+"; }
				
				$sA[2] =~ /(\d+),(\d+)/;
				$start = $1;
				$end   = $2;
				
				if( $sA[1] =~ /^tRNA-(\S+)/ ) { 
					$rna = "tRNA"; 
					$type = $1;
					
					$sA[3] =~ /(\d+)/;
					$anticodonPos = $1;
					if( $strand eq "+" ) {
						$anticodonPos = $start + $anticodonPos - 1; # Verified online
					}
					else {
						$anticodonPos = $end - $anticodonPos + 1; # Verified online
					}
					
					
					$sA[4] =~ /\((\D{2,3})\)/;
					$anticodon = uc($1);
					
					if( $sA[4] =~ /i\((\d+),(\d+)\)/ ) {
						$intronPos = $1;
						$intronLen = $2;
						
						if( $strand eq "+" ) {
							$fragment1End = $start + $intronPos - 2; # Verified online
							$fragment2Start = $start + $intronPos + $intronLen - 1; # Verified online
						}
						else {
							$fragment1End = $end - $intronPos - $intronLen + 1; # Verified online
							$fragment2Start = $end - $intronPos + 2; # Verified online
						}
					}
				}
				elsif( $sA[1] =~ /^tmRNA/ ) {
					$rna = "tmRNA";
					if( $sA[1] =~ /\(p\)/ ) { $permuted = 1; }
					
					$sA[3] =~ /(\d+),(\d+)/;
					if( $strand eq "+" ) {
						$tagStart = $start + $1; # Verified online
						$tagEnd   = $start + $2; # Verified online
					}
					else {
						$tagStart = $end - $1; # Verified online
						$tagEnd   = $end - $2; # Verified online
					}
					
					$sA[4] =~ s/^\s+//;
					$sA[4] =~ s/\s+$//;
					$tagPeptide = $sA[4];
				}
				else { die "Aragorn: No instructions for $sA[1]!\n"; }
				
				
				# Print to outfile in gff format
				my $geneIdentifier = $seqName."_ta".$i;
				my $noteField = "ID=".$geneIdentifier."; ";
				if( $rna eq "tRNA" ) { 
					$noteField .= "product=".$rna."_".$type."_".$anticodon."; ".
									"RNA_Class_ID=".$rna."_".$type."_".$anticodon."; ".
									"Name=".$rna."_".$type."_".$anticodon."; ".
									"Version=".$self->{version}."; tRNA_type=";
					if( $type eq "seC" ) { $noteField .= "SeC(p); "; }
					else { $noteField .= $type."; "; }
					$noteField .= "anticodon=".$anticodon."; anticodon_position=".$anticodonPos;
					$noteField .= "; intron containing" if defined($fragment1End);
				}
				else {
					$noteField .= "product=".$rna."; Name=".$rna."; ".
									"RNA_Class_ID=".$rna."; ".
									"tag_start=".$tagStart."; ".
									"tag_end=".$tagEnd."; ".
									"tag_peptide=".$tagPeptide;
					$noteField .= "; is permuted" if defined($permuted);
				}
				
				my $gffString = $seqName."\t".$self->{version}."\t".$rna."\t".
									$start."\t".$end."\t.\t". # The output has no score
									$strand."\t.\t".$noteField."\n";
				
				print WH $gffString;
				
				if( defined($fragment1End) ) {
					my $exonId = $geneIdentifier.".1";
					$noteField = "Parent=".$geneIdentifier."; ID=".$exonId."; ".
									"product=".$rna."_".$type."_".$anticodon."; ".
									"RNA_Class_ID=".$rna."_".$type."_".$anticodon."; ".
									"Name=".$rna."_".$type."_".$anticodon."; ".
									"Version=".$self->{version}."; tRNA_type=";
					if( $type eq "seC" ) { $noteField .= "SeC(p); "; }
					else { $noteField .= $type."; "; }
					$noteField .= "anticodon=".$anticodon."; anticodon_position=".$anticodonPos;
					$gffString = $seqName."\t".$self->{version}."\texon\t".
									$start."\t".$fragment1End."\t.\t".
									$strand."\t.\t".$noteField."\n";
					print WH $gffString;
					
					$exonId = $geneIdentifier.".2";
					$noteField = "Parent=".$geneIdentifier."; ID=".$exonId."; ".
									"product=".$rna."_".$type."_".$anticodon."; ".
									"RNA_Class_ID=".$rna."_".$type."_".$anticodon."; ".
									"Name=".$rna."_".$type."_".$anticodon."; ".
									"Version=".$self->{version}."; tRNA_type=";
					if( $type eq "seC" ) { $noteField .= "SeC(p); "; }
					else { $noteField .= $type."; "; }
					$noteField .= "anticodon=".$anticodon."; anticodon_position=".$anticodonPos;
					$gffString = $seqName."\t".$self->{version}."\texon\t".
									$fragment2Start."\t".$end."\t.\t".
									$strand."\t.\t".$noteField."\n";
					print WH $gffString;
				}
			}
		}
	}
	close CMD;
	
	close WH;
}


sub setSequenceTopology {
	my ($self, $topology) = @_;
	
	if( lc($topology) eq "circular" ) { $self->{options} .= "-c "; }
	elsif( lc($topology) eq "linear" ) { $self->{options} .= "-l "; }
}


sub setIntronSize {
	my ($self, $intronSize) = @_;
	
	$self->{options} .= "-i$intronSize ";
}


sub setTranslationTable {
	my ($self, $transTable) = @_;
	
	$self->{options} .= "-gc$transTable ";
}


sub getVersion{
	my ($self) = @_;
	
	if( !defined($self->{version}) ) {
		$self->{version} = geneCallerVersions::getVersion("aragorn");
	}
	
	return $self->{version};
}
1;
