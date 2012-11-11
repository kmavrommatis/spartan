package Infernal;

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/..";
use lib "$FindBin::Bin/../../helpingHands";
use geneCallerVersions;
use CommonFunc;
use parser;

use File::Basename;

use Bio::SeqIO;


# Change after changing the code
my $version = "2.2";

# find hits from the Rfam model collection
# this script reads the data from the file Rfam.desc, 
# which is no longer provided by the rfam database
# instead I download the html page with the description from the url
# http://rfam.janelia.org/browse.html
#				^
#				|
# MH:	stuck at version 10.0
# MH: Wrote script that creates desc file from downloaded (ftp) db files

sub new{
	my $self={};
    	my ($class,@arguments)=@_;
        bless($self,$class);
        
        for (my $i=0;$i<scalar(@arguments);$i++){
        	if(	$arguments[$i] eq '-inputFile' or
        		$arguments[$i] eq '-outputFile' or
        		$arguments[$i] eq '-outputFh' or
        		$arguments[$i] eq '-blastDB' or
        		$arguments[$i] eq '-path'){
        			my $a= substr( $arguments[$i] , 1 );
        			$self->{ $a  }=$arguments[++$i];
        		}
        } 
        
        # make sure that mandatory information is given
       	if(!defined($self->{ inputFile })){die "findrRNAs: No input file was provided\n";} 
       	if(defined($self->{ outputFile}) and !defined($self->{outputFh})){
       		open $self->{ outputFh },">".$self->{outputFile} or 
					die "Cannot create $self->{outputFile}\n";
       	}
       	if(!defined($self->{path})) { 
			my ($t, $fn) = CommonFunc::splitFn($self->{inputFile});
			$self->{path}= $t; 
			if(! $self->{path} or $self->{path} eq '' or $self->{path} eq '/'){
				$self->{path}='/tmp/';
			}
		}
		if( $self->{path} !~ /\/$/ ) {
			$self->{path} .= "/";
		}
		
        $self->getVersion();
        
        if( !defined($self->{blastDB}) ) {
        	$self->{blastDB} = $ENV{ rfamDir }."Rfam.fasta";
        }
       	
        
        $self->{counter} = 0;	
       
		$self->{blastMinus_b_optionValue} = 100000;
		  
        return $self;
}


sub excludeModels {
	my ($self, $models2exclude) = @_;
	
	if( $models2exclude eq "db" ) {
		
		# Create hash for allowed Rfam overlaps
		 %{$self->{excludeRNA}} = ();
		 %{$self->{ncbi_locus_type}} = ();
		open(READ, "<", $ENV{ rfam_tbl }) or die "Cannot open $ENV{ rfam_tbl }: $!";
		while(<READ>) {
			chomp $_;
			my @splitArray = split("\t", $_);
			if( $splitArray[10] ne "Rfam" ) { next; }
			if( $splitArray[12] eq "Yes" ) { $self->{excludeRNA}->{$splitArray[0]} = 1; }
			$self->{ncbi_locus_type}->{ $splitArray[0]}=$splitArray[8];
		}
		close(READ);
	}elsif( $models2exclude eq 'regular'){
		 %{$self->{excludeRNA}} = (  
        				"RF00001"=>1, "RF00002"=>1, "RF00005"=>1, "RF00177"=>1,
						"RF01315"=>1, "RF01316"=>1, "RF01317"=>1, "RF01318"=>1, "RF01319"=>1,
						"RF01320"=>1, "RF01321"=>1, "RF01322"=>1, "RF01323"=>1, "RF01324"=>1,
						"RF01325"=>1, "RF01326"=>1, "RF01327"=>1, "RF01328"=>1, "RF01329"=>1,
						"RF01330"=>1, "RF01331"=>1, "RF01332"=>1, "RF01333"=>1, "RF01334"=>1,
						"RF01335"=>1, "RF01336"=>1, "RF01337"=>1, "RF01338"=>1, "RF01339"=>1,
						"RF01340"=>1, "RF01341"=>1, "RF01342"=>1, "RF01343"=>1, "RF01344"=>1,
						"RF01345"=>1, "RF01346"=>1, "RF01347"=>1, "RF01348"=>1, "RF01349"=>1,
						"RF01350"=>1, "RF01351"=>1, "RF01352"=>1, "RF01353"=>1, "RF01354"=>1,
						"RF01355"=>1, "RF01356"=>1, "RF01357"=>1, "RF01358"=>1, "RF01359"=>1,
						"RF01360"=>1, "RF01361"=>1, "RF01362"=>1, "RF01363"=>1, "RF01364"=>1,
						"RF01365"=>1, "RF01366"=>1, "RF01367"=>1, "RF01368"=>1, "RF01369"=>1,
						"RF01370"=>1, "RF01371"=>1, "RF01372"=>1, "RF01373"=>1, "RF01374"=>1,
						"RF01375"=>1, "RF01376"=>1, "RF01377"=>1, "RF01378"=>1, "RF01379"=>1,
						"RF01959"=>1, "RF01960"=>1
					  );
	}
}


sub constrainModelsTo {
	my ($self, $models2useStr) = @_;
	
	my @models2use = split(/,\s*/, $models2useStr);
	if( scalar(@models2use) > 0 ) {
		for( @models2use ) {
			if( lc($_) eq "trnas" ) {
				$self->{models2use}{RF00005} = 1;
				$self->{models2use}{RF01852} = 1;
			}
			# Possibility to add other models as well
		}		
	}
}



sub runPrediction {
	my ($self, $rnaFn, $tmp_dir, $useGFF) = @_;
	$self->{gffOutputArray}=[];
	my %modelData;
	my %seq2model;
	$self->loadModelData(\%modelData, \%seq2model);
	my @blastHits=();
	$self->runBlast(\@blastHits, \%seq2model);
	my %sequenceHits;
	$self->findBoundaries(\%modelData,\%seq2model,\@blastHits,\%sequenceHits);
	$self->runRfam(\%sequenceHits, \%modelData);
	return "";
}



############################
# Split the input file in one sequence at a time 
# and run each model for this sequence
sub runRfam{
	my ($self, $sequenceHits, $modelData) = @_;
	my $fh=$self->{ outputFh };

	# Take one sequence at a time
	my $inStream = Bio::SeqIO->new(-file=>$self->{inputFile},-format=>'fasta');
	while (my $seq = $inStream->next_seq()) {
		# Save the sequence to a temporary file
		my $seqName = $seq->display_id();
		#print "runRfam: Processing sequence $seqName\n";
		my $models = $self->findModels($seqName, $sequenceHits);
		if (scalar(keys(%{ $models })) == 0) { next; }
		# Run cmsearch for each model
		foreach my $model2search(keys (%{ $models } )) {
			my ($sequence, $hitStart, $hitEnd) = split("\t", $models->{$model2search});
			my $offsetStart;
			my $offsetEnd;
			if ($hitStart <= 200) { $offsetStart = $hitStart-1; }
			else { $offsetStart = 200; }
			
			if ($hitEnd >= $seq->length() - 200) {
				$offsetEnd = $seq->length() - $hitEnd-1
			}
			else { $offsetEnd = 200; }
			
			my $fragmentStart = $hitStart - $offsetStart;
			my $fragmentEnd   = $hitEnd   + $offsetEnd;
			my $tempSeq = $seq->subseq($fragmentStart, $fragmentEnd);

			my $tempSeqName = $seqName;
			$tempSeqName =~ s/\//_/g;
			my $tempSeqFn = $self->{path}."$tempSeqName.$$.fna";
			my $tempSeqObj = Bio::Seq->new(-seq=>$tempSeq, -id=>$seqName);
			my $outStream = Bio::SeqIO->new(-file=>">$tempSeqFn", -format=>'fasta');
			$outStream->write_seq($tempSeqObj);
			my $outputFn = $self->{path}."$tempSeqName.$model2search.cmsearch";

#			my $cmd = $ENV{ cmsearchBin }." -g --ga --tabfile $outputFn ".
#						$ENV{ rfamDir }."models/$model2search.cm $tempSeqFn".
#						" 1> /dev/null";
						
			my $cmd = $ENV{ cmsearchBin11 }." -g --cut_ga --tblout $outputFn --cpu 0 -Z 1 ".
						$ENV{ rfamDir11 }."models/$model2search.cm $tempSeqFn".
						" 1> /dev/null";			
#			print "runRfam: executing $cmd\n";
			system($cmd);
		 	if ( $? ) { die "Command failed: $cmd : $!"; }

			$self->parseCMresults($model2search, $outputFn, $modelData,
									$fh, $fragmentStart, $fragmentEnd);
			
			unlink($outputFn);
			unlink($tempSeqFn);
			
		}
	}
	
	
	# now we have everything in the $self->{gffOutputArray}
	# resolve overlaps
	@{$self->{gffOutputArray}} = sort{ $a->[0] cmp $b->[0] || $a->[3] <=> $b->[3] } @{$self->{gffOutputArray}};
	for(my $i=0;$i<scalar( @{$self->{gffOutputArray}} )-1;$i++){
		if(!defined($self->{gffOutputArray}->[$i])){next;}
		for(my $j=$i-5;$j<$i+5;$j++){
			if( $j<0 or $j>scalar( @{$self->{gffOutputArray}} )-1  ){next;}
			if($j==$i){next;}	
			if(!defined($self->{gffOutputArray}->[$j])){next;}
			if(!defined($self->{gffOutputArray}->[$i])){next;}
#			print "Comparing:\n'".join(",", @{$self->{gffOutputArray}->[$i]}).
#					"'\n'".join(",", @{$self->{gffOutputArray}->[$j]})."'\n\n";
			if($self->{gffOutputArray}->[$i]->[0] ne $self->{gffOutputArray}->[$j]->[0]){next;}
			if( CommonFunc::isOverlapping( $self->{gffOutputArray}->[$i]->[3], $self->{gffOutputArray}->[$i]->[4],
								$self->{gffOutputArray}->[$j]->[3], $self->{gffOutputArray}->[$j]->[4]  ) > 15){
				if($self->{gffOutputArray}->[$i]->[5] < $self->{gffOutputArray}->[$j]->[5]){ undef($self->{gffOutputArray}->[$j]); }
				elsif($self->{gffOutputArray}->[$j]->[5] < $self->{gffOutputArray}->[$i]->[5]){ undef($self->{gffOutputArray}->[$i]); }
				else{
					if( abs($self->{gffOutputArray}->[$i]->[4]-$self->{gffOutputArray}->[$i]->[3]) >= 
							abs($self->{gffOutputArray}->[$j]->[4]-$self->{gffOutputArray}->[$j]->[3]) ) {
						undef($self->{gffOutputArray}->[$j]);
					}
					else { undef($self->{gffOutputArray}->[$i]); }
				}
			}
		}
	}

	# store the output	
	foreach my $line( @{$self->{gffOutputArray}}){
		if(!defined ($line)){next;}
		print $fh join("\t", @$line),"\n";
	}
	close $fh;
}



sub parseCMresults {
	my ($self, $model, $fn, $modelData, $output,
				$fragmentStart, $fragmentEnd) = @_;
				
	
	open (CM,"$fn") or die ("Cannot open file $fn\n");
#	print "parsing file $fn\n";
	while(my $line=<CM>) {
#		print $fn.":: ".$line."\n";
		chomp $line;
		
		if ($line eq "" or $line=~/^#/) { next; }
#		my ($j, $modName, $seq, $seqStart, $seqEnd, 
#			$queryStart, $queryEnd, $bit, $evalue, $gc) = split(m/\s+/, $line);
			
		my ($seq, $accession, $queryName, $modName, $modelType, $queryStart,$queryEnd, $seqStart,$seqEnd, 
		    $strand, $trunc, $pass, $gc, $bias, $bit, $evalue, $inc, $description)=split(m/\s+/, $line);
		
		
		my $modelSize = $modelData->{$model}->{seqLen};
		my $modelName = $modelData->{$model}->{modelFunction};
		my $modelDesc = $modelData->{$model}->{description};
		# Accept the RNA model
		if(abs($queryStart-$queryEnd) >= $modelSize * 0.8) {
			
			my $LseqStart = $seqStart + $fragmentStart-1;
			my $LseqEnd   = $seqEnd + $fragmentStart  -1;
#			print "parseCMRresults: Found a hit at position $seqStart - $seqEnd ".
#							"of the fragment $fragmentStart - $fragmentEnd\n";
			
			$self->{counter}++;
			my $strand;
			my $s;
			my $e;
#			my $type = "misc_RNA";
#			
#			if ($modelDesc =~/T-box/ or $modelDesc =~/ydaO\/yuaA/ or 
#				$modelDesc=~/L10_leader/ or $modelDesc=~/Hfq binding/ or 
#				$modelDesc=~/Pseudoknot/) {
#				$type = 'misc_feature';
#			}
#			elsif($modelDesc =~/PyrR/ or $modelDesc=~/FMN/ or $modelDesc =~/TPP/ or 
#									$modelDesc =~/Purine/ or $modelDesc=~/SAM/ or 
#									$modelDesc =~/Cobalamin/ or lc($modelDesc)=~/riboswitch/)	{
#				$type = 'misc_bind';
#			}
#			elsif($modelDesc =~/tmRNA/){
#				$type = 'tmRNA';
#			}
#			elsif($modelDesc =~/tRNA/){
#				$type = 'tRNA';
#			}
			my $type= $self->{ncbi_locus_type}->{ $model  };
			if ($LseqStart < $LseqEnd) { 
				$strand = "+";
				$s = $LseqStart;
				$e = $LseqEnd;
			}
			else {
				$strand = "-";
				$s = $LseqEnd;
				$e = $LseqStart;
			}		
			my $note="ID=$seq.RFAM." . $self->{counter}."; ".
								"Version=". $self->{version} ."; ".
								"RNA_Class_ID=$model; ".
								"Model=$modelName; ".
								"product=$modelDesc";
			push @{$self->{ gffOutputArray }},
								[$seq , $self->{version} ,$type ,$s ,$e, 
													$bit, $strand, 0, $note];
			
		}
		else {
#			print "parseCMRresults: hit is not long enough\n";
		}
	}
	close CM;
}



sub getVersion{
	my ($self) = @_;
	
	if(!defined($self->{version})){
		$self->{version} = geneCallerVersions::getVersion("cmsearch11");
	}
	
	return $self->{version};
}



###############################
# Find the models that have hits on 
# that sequence
sub findModels{
	my ($self, $seqName, $sequenceHits) = @_;
	
	my %models;
	foreach my $model(keys %{ $sequenceHits }) {
		if (defined($self->{excludeRNA}->{$model})) { next; }
		my $hits = $sequenceHits->{$model};
		foreach my $h(@{ $hits }) {
			if ($h->[0] ne $seqName) { next; }
			my ($seq_name,$seq_start,$seq_end) = (undef, 100000000000000000000, -1);
			if (defined($models{$model})) {
				($seq_name,$seq_start,$seq_end) = split("\t",$models{$model});
			}
			if ($seq_start > $h->[1]) { $seq_start=$h->[1]; }
			if ($seq_end   < $h->[2]) { $seq_end  =$h->[2]; }
			$models{$model} = "$h->[0]\t$seq_start\t$seq_end";
		}
	}
	return \%models;
}



#############################
# View the results of the hits
# and possilble rfam models
sub printModels{
	my ($self, $sequenceHits) = @_;
	foreach my $model(keys %{ $sequenceHits }) {
		if (defined($self->{excludeRNA}->{$model})) { next; }
		print "$model:\n";
		my $hits = $sequenceHits->{$model};
		foreach my $h(@{ $hits }) {
			print "$model\t".$h->[0]."\t".$h->[1]."\t",$h->[2];
			print "\t".$h->[3],"\t".$h->[4],"\n";
		}
	}
}



#############################
# Process the blast output
# and find the boundaries of the sequence
# that needs to be searched for each model
sub findBoundaries {
	my($self, $modelData, $seq2model, $blastHits, $sequenceHits) = @_;
	
	# for each of the hits find the Rfam model that has a hit
	for(my $i=0; $i<scalar(@{ $blastHits }); $i++) {
		my $seqName = $blastHits->[$i][0];
		
		my $modelHit = $seq2model->{ $blastHits->[$i][1] };
		my ($seqStart, $seqEnd) = ($blastHits->[$i][4], $blastHits->[$i][5]);
		my ($modelStart, $modelEnd) = ($blastHits->[$i][6], $blastHits->[$i][7]);
		
		# If this model has not been identified earlier add it in the list
		if (!defined($sequenceHits->{$modelHit})) {
			$sequenceHits->{$modelHit} = ();
		}

		push @{ $sequenceHits->{$modelHit} }, [$seqName,$seqStart,$seqEnd,$modelStart,$modelEnd];		
	}
}



##############################
# Read the file with the fasta sequences of the models
# and remember the name, model function and size for each sequence
sub loadModelData {
	my ($self, $modelData, $seq2model) = @_;
	#my $modelFunction={};
	my $modelStream = Bio::SeqIO->new(-file=>$self->{blastDB}, -format=>'fasta');
	while(my $seq = $modelStream->next_seq()) {
		my $name = $seq->display_id();
#		my $model = $seq->desc(); # OLD FORMAT!
#		$model =~/(\S+)\;(\S+)\;/; # OLD FORMAT!
		$name =~/(\S+)\;(\S+)\;\S+/;
		my ($rfamName, $function) = ($1, $2);
		if( defined($self->{models2use}) ) {
			if( defined($self->{models2use}{$rfamName}) ) {
				$seq2model->{$name} = $rfamName;
			}
		}
		else {
			$seq2model->{$name} = $rfamName;
		}
		
	}

	my $rfamDesc = $ENV{ rfamDir }."Rfam.desc";
	open (DESC, $rfamDesc) or die "Cannot open $rfamDesc\n";

	while( <DESC> ) {
		chomp($_);
		if( $_ eq "" ) { next; }
		
		my ($rfamAcc, $rfamName, $rfamDesc, $avgSeedLen) = split "\t", $_;
		
		if( defined($self->{models2use}) ) {
			if( defined($self->{models2use}{$rfamAcc}) ) {
				$modelData->{ $rfamAcc } = {	modelFunction=>$rfamName, 
						 		 				seqLen=>$avgSeedLen, 
						  						description=>$rfamDesc };
			}
		}
		else {
			$modelData->{ $rfamAcc } = {	modelFunction=>$rfamName, 
						 		 			seqLen=>$avgSeedLen, 
						  					description=>$rfamDesc };
		}
	}
	close DESC;
}



################################
# Compare the input file to the database of 
# Rfam using blast
sub  runBlast {
	my ($self, $blastHits, $seq2model) = @_;
	my $blastOut = $self->{path}."runInfernal.$$.blout";
						
	my $blastcmd = "$ENV{ blastallBin } -p blastn -i ".$self->{inputFile}." -d ".$self->{blastDB}.
	
						" -e 10 -W7 -F F -K 5 -m8 -b ". $self->{blastMinus_b_optionValue};

#	print "runInfernal: Running blast for fast search \n $blastcmd\n";	
	
	## parse the blast output
#	open BLAST, $blastOut or die "runInfernal: Cannot open the blast output file $blastOut\n";
#	print "runInfernal: parsing the blast output\n";
	my $numberOfHits = 0;
	my $rerun = "No";
	open(BLASTCMD, "$blastcmd|") or die "Command failed: $blastcmd : $!";
	while (my $line = <BLASTCMD>) {
		if ($numberOfHits ==  $self->{blastMinus_b_optionValue}) {
			$rerun = "Yes";
			@{ $blastHits } = ();
			last;
		}
		print $line;
		chomp $line;
		my ($query, $hit, $pid, $alLen, $i, $j, 
				$queryStart, $queryEnd, $hitStart, 
						$hitEnd, $evalue, $bitScore) = split("\t", $line);
		$hit=~/^(RF\d+);/; my $rfid=$1;
		if($self->{excludeRNA}->{$rfid} ){ print "This line will be excluded \n";next; }
		if( defined($self->{models2use}) ) {
			if( defined($seq2model->{$hit}) and defined($self->{models2use}{$seq2model->{$hit}}) ) {
				push @{ $blastHits }, [$query,$hit,$pid,$alLen,$queryStart,$queryEnd,$hitStart,$hitEnd];
				$numberOfHits++;
			}
		}
		else {
			push @{ $blastHits }, [$query,$hit,$pid,$alLen,$queryStart,$queryEnd,$hitStart,$hitEnd];
			$numberOfHits++;
		}
	}
	close(BLASTCMD);
	
	if ($rerun eq "Yes") {
		$self->{blastMinus_b_optionValue} *= 10;
#		print "runInfernal: Too many hits! Rerunning blast with -b ".$self->{blastMinus_b_optionValue}."\n";
		$blastcmd = "$ENV{ blastallBin } -p blastn -i ". $self->{inputFile}." -d ".$self->{blastDB}.
						" -e 10 -W7 -F F -K 5 -m8 -b  ".$self->{blastMinus_b_optionValue};
		open(BLASTCMD, "$blastcmd|") or die "Command failed: $blastcmd : $!";
		while (my $line = <BLASTCMD>) {
			chomp $line;
			my ($query, $hit, $pid, $alLen, $i, $j, 
					$queryStart, $queryEnd, $hitStart, 
							$hitEnd, $evalue, $bitScore) = split("\t", $line);
			if( defined($self->{models2use}) ) {
				if( defined($self->{models2use}{$seq2model->{$hit}}) ) {
					push @{ $blastHits }, [$query,$hit,$pid,$alLen,$queryStart,$queryEnd,$hitStart,$hitEnd];
				}
			}
			else {
				push @{ $blastHits }, [$query,$hit,$pid,$alLen,$queryStart,$queryEnd,$hitStart,$hitEnd];
			}
		}
		close(BLASTCMD);
	}
	
	unlink("$self->{inputFile}.*");
	unlink("formatdb.log");
#	print "runBlast: After parsing ".scalar(@{ $blastHits })." lines with hits were detected\n";
#	close $blastOut;
	#unlink("$blastOut") if (-e $blastOut);
}
1;
