package rRNA;

use strict;
use warnings;
use Log::Log4perl;
use FindBin;
use lib "$FindBin::Bin";
use geneCallerVersions;
use CommonFunc;
use hmmsearchParser;
use cmsearchParser;
use IPC::Open2;

use Bio::SeqIO;

my $logger =  Log::Log4perl->get_logger('rRNAm');
my $isMetagenome;

my $intronMaxSize=	1500;
my($gapOnModelMax, $gapOnModelMin);
my $minExonDistance=30;
my $partialFlag='off';
sub findRRNA{
	my ($contigFn, $domains, $type, $dataArray, $partial) = @_;
	
	
	my $counter = 0;
	my $count=0;my $count2=0;
	my $removedID={};
	my $taxassignment={};
	if(defined($partial) and $partial ==1 ){$partialFlag = 'on';}
	foreach my $domain (@$domains) {
		$logger->debug("Looking for rRNA with hmmer file:$contigFn, Domain:$domain, RNA type:$type");
		$count = findRNA_hmmer($contigFn, $domain, $type, $dataArray,$taxassignment);
		$counter += $count;
		$logger->debug("Looking for rRNA with cmsearch file:$contigFn, Domain:@$domains, RNA type:$type");
		$count2 =findRNA_infernal($contigFn, $domain, $type, $dataArray,$taxassignment);
		$counter += $count2;	
	}
	
	
	if($count>0){
		$logger->info( "Hmmer found $count rRNAs ($type)!");
	}
	if($count2>0){
		$logger->info( "Cmsearch found $count2 rRNAs ($type)!");
	}
	
}


sub findRNA_hmmer{
	
	$ENV{ rrnaHmmsDir }=$FindBin::Bin ."/hmms/";

	my ($contigsFn, $domain, $type, $dataArray, $taxassignment) = @_;
	if($type eq 'tsu'){return 0;}
	$domain=lc($domain);
	my $programVersion = geneCallerVersions::getVersion("hmmer");
	$logger->debug( "Predicting rRNA genes with HMMs and $programVersion\n");
	my $evalue_cutoff=0.01;
	if($type eq 'tsu'){$evalue_cutoff=0.1;}
	my $filePrefix="bac";
	if($domain eq 'b'){$filePrefix="bac"}
	elsif($domain eq 'a'){$filePrefix="arc"}
	elsif($domain eq 'e'){$filePrefix="euk"}
#print "Domain is $domain so prefix becomes $filePrefix\n";
	my $hmmFn=$type."_".$filePrefix.".hmm";

	my $rnaCounter = 0;
	my $sequenceName="";
	my $parser='off';
	my $description;
	my $modelSize=0;
	my $modelName;
	my $strand;
	my $sequenceSize = &getSequenceSize($contigsFn);

	my $resultsFull={}; # { alifrom } => \@rec
	my $resultsPartial={};
	my $specificFlag;
	if($type eq 'tsu'){$specificFlag='--max --incE 0.1';}
	else{$specificFlag='--nobias --incE 0.01 --domE 0.1 ';}
	my $cmdline = $ENV{ hmmsearchBin }." -E $evalue_cutoff --noali $specificFlag -Z 1 --cpu 0 ".
					$ENV{ rrnaHmmsDir }."$hmmFn -";
	$logger->info( "Executing command $cmdline\n");
	
	
	runHMMPrediction( $cmdline , $contigsFn, $sequenceSize ,$resultsFull, $resultsPartial);
	gatherPhylogeny( $resultsPartial, $taxassignment);
	gatherPhylogeny( $resultsFull, $taxassignment);
	#resultsFull have all the records that have full name
	# resultsPartial have all the records that have partial hits
	# 	we need to check if there are any exons in the partial list
	# if a record appears to have introns it is merged and then added
	# in the resultsFull list.
	processPartial($resultsPartial,$resultsFull,$type);
	
	#process the full models
	my $coordinates=[];
	my $gene_id=1;
	foreach my $k(keys (%$resultsFull) ) {
		my $ar=$resultsFull->{$k};
		foreach my $r( @$ar) {
			my $description= $r->[2];
			my $geneIdentifier=$r->[0]."_".uc($domain)."_rh". substr($type,0,1) . $gene_id;
			my $groupData="ID=$geneIdentifier; Version=$programVersion; RNA_Class_ID=$description; ".
							"Type=$description; product=$description; Name=$description";
			$gene_id++;
			push @$dataArray, [$r->[0],  # name of the sequence 
								$programVersion,	# name of the method 
								'rRNA' ,		# type of molecule 
								$r->[19], $r->[20], # coordinates of molecule on the contig
								$r->[9],	# score
								$r->[5],	# strand
								".",
								$groupData
							   ];
			# treat exons			
			if($r->[24] and scalar( @{$r->[24]} )>1  ){
				for(my $i=scalar( @{$r->[24]} )-1 ;$i>=0;$i-- ){
					my $e = $r->[24]->[$i];
					my $exonId="$geneIdentifier." . ($i+1);
					
					my $groupData="ID=$exonId;Parent=$geneIdentifier; Version=$programVersion; ".
							"Type=$description; product=$description; Name=$description";
						push @$dataArray, [$r->[0],  # name of the sequence 
								$programVersion,	# name of the method 
								'exon' ,		# type of molecule 
								$e->[0], $e->[1], # coordinates of exon on the contig
								$r->[9],	# score
								$r->[5],	# strand
								".",
								$groupData
							   ];
				}
				
			}
			$rnaCounter++;
		}
	}

	return $rnaCounter;
}

# in the case of tsu the infernal model cannot predict
# the phylogenetic signal of the hit.
# We are storing all the hits that hmmsearch finds, 
# and if we have a hit that is kept from infernal
# we apply this phylogenetic information to it
sub gatherPhylogeny{
	my($resultsPartial, $taxassignment)=@_;
	
	foreach my $seqName( keys(%$resultsPartial) ){
			
			foreach my $r( @{ $resultsPartial->{$seqName}   }   ){
				
				# instead of score (ie bitscore) use the size of the alignment on the model
				my ($envStart,$envEnd,$score, $type)=($r->[19], $r->[20],$r->[14]-$r->[13], $r->[2] );
				if(!defined( $taxassignment->{ $seqName })){ 
					$taxassignment->{ $seqName }=[];
				}
#				print "******** loagind taxassignment with $seqName, $envStart,$envEnd,$score, $type \n";
				push @{ $taxassignment->{ $seqName } }, [$envStart,$envEnd,$score, $type];
			}
			
	}
}
sub checkPhylogeny{
	my( $taxassignment, $targetName, $SeqFrom, $SeqTo)=@_;
	my @temp;
	
#	print "****** checking for $targetName in the region $SeqFrom to $SeqTo\n";
	
	if(!defined( $taxassignment->{ $targetName }  )){ return undef;}
	foreach my $r( @{ $taxassignment->{ $targetName }  }  ){
		if( CommonFunc::isOverlapping( $r->[0],$r->[1] , $SeqFrom ,$SeqTo  ) > 50 ){
			push @temp, [ $r->[ 2 ], $r->[3] ];
		}
	}
	# the @temp array holds the score and annotation of the hits in that location
	# we want to get the best hit, i.e. best score
	@temp = sort{ $b->[0] <=> $a->[0]} @temp;
	return $temp[0][1];
}



sub findRNA_infernal{
	my ($contigsFn, $domain , $type, $dataArray, $taxassignment) = @_;
# 	print "findRNA_infernal: The domai is $domain\n";
	if($type ne 'tsu'){return 0;}
	my $programVersion = geneCallerVersions::getVersion("cmsearch11");
	$logger->debug( "Predicting rRNA genes with covariance models and $programVersion\n");
	my $evalue_cutoff=0.1;

	my $filePrefix="bac";
	if(lc($domain) eq 'b'){$filePrefix="bac"}
	elsif(lc($domain) eq 'a'){$filePrefix="arc"}
	elsif(lc($domain) eq 'e'){$filePrefix="euk"}
#print "Domain is $domain so prefix becomes $filePrefix\n";
	my $covFn=$ENV{rrnaHmmsDir}."/".$type."_".$filePrefix.".cm";

	#my $covFn=$ENV{rrnaHmmsDir}."/5s.cm";

	my $rnaCounter = 0;
	my $sequenceName="";
	my $parser='off';
	my $description;
	my $modelSize=0;
	my $modelName;
	my $strand;
	my $sequenceSize = &getSequenceSize($contigsFn);

	my $resultsFull={}; # { alifrom } => \@rec
	my $resultsPartial={};
	my $specificFlag;
	if($type eq 'tsu'){$specificFlag=' --cut_ga ';}
	else{$specificFlag='--nobias';}
	my $cmdline =  $ENV{ cmsearchBin11 }.
				" --notextw $specificFlag  --cpu 0  ".
				"$covFn $contigsFn";
	$logger->info( "Executing command $cmdline");
	runCMPrediction( $cmdline ,$contigsFn,  $sequenceSize ,$resultsFull, $resultsPartial);
	processPartial($resultsPartial,$resultsFull,$type);
	my $gene_id=1;
	foreach my $k(keys (%$resultsFull) ) {
		my $ar=$resultsFull->{$k};
		foreach my $r( @$ar) {
			my $description= $r->[2];
			
			# check if this hit exists as a partial hit from hmmsearch 
			# and if so transfer its annotation to the cmsearch hit.
			# this will transfer the information about the taxonomic assignment of the 
			# hit
# 			my $hmmDescription=checkPhylogeny( $taxassignment, $r->[0], $r->[19], $r->[20]);
# 			if(defined($hmmDescription)){ $description= $hmmDescription;}
		
			
			my $geneIdentifier=$r->[0]."_rc". substr($type,0,1) . $domain. $gene_id;
			my $groupData="ID=$geneIdentifier; Version=$programVersion; RNA_Class_ID=$description; ".
							"Type=$description; product=$description; Name=$description";
			$gene_id++;
			push @$dataArray, [$r->[0],  # name of the sequence 
								$programVersion,	# name of the method 
								'rRNA' ,		# type of molecule 
								$r->[19], $r->[20], # coordinates of molecule on the contig
								$r->[9],	# score
								$r->[5],	# strand
								".",
								$groupData
							   ];
			# treat exons			
			if($r->[24] and scalar( @{$r->[24]} )>1  ){
				for(my $i=scalar( @{$r->[24]} )-1 ;$i>=0;$i-- ){
					my $e = $r->[24]->[$i];
					my $exonId="$geneIdentifier." . ($i+1);
					
					my $groupData="ID=$exonId;Parent=$geneIdentifier; Version=$programVersion; ".
							"Type=$description; product=$description; Name=$description";
						push @$dataArray, [$r->[0],  # name of the sequence 
								$programVersion,	# name of the method 
								'exon' ,		# type of molecule 
								$e->[0], $e->[1], # coordinates of exon on the contig
								$r->[9],	# score
								$r->[5],	# strand
								".",
								$groupData
							   ];
				}
				
			}
			$rnaCounter++;
		}
	}
#	open (CMSEARCH, "$cmdline|") or $logger->logdie( "Command failed: $cmdline $!");
#	my $cm=cmsearchParser->new( '-cmOutStream'=>\*CMSEARCH );
#	my $gene_id=1;
#	while(my $hit=$cm->nextHit()){
#		
#		
#		
#		
#		if($hit->{Qual} eq '?'){next;}
#		my $description= $cm->{queryDescription};
#		
#		# check if this hit exists as a partial hit from hmmsearch 
#		# and if so transfer its annotation to the cmsearch hit.
#		# this will transfer the information about the taxonomic assignment of the 
#		# hit
#		my $hmmDescription=checkPhylogeny( $taxassignment, $hit->{targetName}, $hit->{SeqFrom}, $hit->{SeqTo});
#		if(defined($hmmDescription)){ $description= $hmmDescription;}
#		
#		my $geneIdentifier=$hit->{targetName}."_cm" . $gene_id;
#		my $groupData="ID=$geneIdentifier; Version=$programVersion; ".
#						"RNA_Class_ID=$description; Type=$description;product=$description; Name=$description";
#		$gene_id++;
#		push @$dataArray, [$hit->{targetName},  # name of the sequence 
#							$programVersion,	# name of the method 
#							'rRNA' ,		# type of molecule 
#							$hit->{SeqFrom}, $hit->{SeqTo}, # coordinates of molecule on the contig
#							$hit->{Score},	# score
#							$hit->{SeqStrand},	# strand
#							".",
#							$groupData
#						   ];
#		$rnaCounter++;
#	}
#	close(CMSEARCH);

	return $rnaCounter;
}



sub getSequenceSize{
	my ($contigFn) = @_;
	my $returnHash = {};
	my $inStream = Bio::SeqIO->new(-file=>$contigFn, -format=>'fasta');
	while(my $seqObj = $inStream->next_seq()) {
		my $length = $seqObj->length();
		$returnHash->{$seqObj->display_id()} = $length;
	}
	return $returnHash;
}


sub processPartial{
	my ($resultsPartial, $resultsFull,$type)=@_;
	if($type =~/ssu/){$gapOnModelMax=30;$gapOnModelMin=-300;}
	if($type =~/lsu/){$gapOnModelMax=50;$gapOnModelMin=-500;}
	if($type =~/tsu/){$gapOnModelMax=15;$gapOnModelMin=0;}
	$logger->trace("------------ process partial ------------------");
	#process the partial hits
	foreach my $k(keys (%$resultsPartial) ){
		
		$logger->trace( "Processing hits from sequence $k\n");
		my $ar=$resultsPartial->{$k}; 
		
		my $coordinates=[];
		
		@$ar=sort{ $a->[16] <=> $b->[16]} @$ar; # sort based on the alifrom coordinate
#		for(my $i=0;$i<scalar(@$ar);$i++){
#			
#		}
		$logger->trace("Checking hits for merging");
		#find successive hits to the same model, that have distance shorter than a determined value
		#	[0] $sequenceName,$sequenceSize->{$sequenceName},$description,$modelName,$modelSize,$strand,
		#	[6] $sp,$no,$test,$score,$bias,$c_eval,$i_eval,[13] $hmmfrom,[14] $hmmto,$dom1,
		#       [16] $alifrom,$alito,$dom2,$envfrom,$envto,$dom3,$acc, [coordinates]
		for(my $i=1;$i<scalar(@$ar);$i++){
			my $gapOnModel;
			my $intronSize;
			$logger->trace("Comparing " , join("\t", @{ $ar->[$i-1]}[0..22] ) );
			$logger->trace("       to " , join("\t", @{ $ar->[$i  ]}[0..22] ) );
#				
#			if( $ar->[$i]->[13]  > $ar->[$i-1]->[14]){
#				$gapOnModel=  $ar->[$i]->[13]- $ar->[$i-1]->[14] ;
#			}else{
#				$gapOnModel=  $ar->[$i-1]->[13]- $ar->[$i]->[14] ;
#			}
			$gapOnModel=findGap( [ $ar->[$i]->[13], $ar->[$i]->[14] ] , [$ar->[$i-1]->[13], $ar->[$i-1]->[14] ] );
#			if( $ar->[$i]->[5] eq '+' ){
#				$intronSize=$ar->[$i]->[16] - $ar->[$i-1]->[17] ;
#			}else{
#				$intronSize=$ar->[$i]->[16] - $ar->[$i-1]->[17] ;
#			}
			$intronSize=findGap(  [$ar->[$i]->[16] , $ar->[$i]->[17]], [$ar->[$i-1]->[16] , $ar->[$i-1]->[17]] );
			
			
			$logger->trace(  "Intron size " ,$ar->[$i]->[16]," - ", $ar->[$i-1]->[17] ,"=", $intronSize," maximum allowd size is $intronMaxSize");
			$logger->trace(  "Gap on model $ar->[$i-1]->[13]/$ar->[$i]->[13] to $ar->[$i-1]->[14]/$ar->[$i]->[14] =", $gapOnModel, " max/min for gap is $gapOnModelMax/$gapOnModelMin" );
			
			if (  ($ar->[$i]->[3] eq $ar->[$i-1]->[3]) and   #the hits are to the same model
				  (($intronSize < $intronMaxSize) and   # the size of the intron should be reasonable 
				  ( $intronSize >0)) and
				  (($gapOnModel <$gapOnModelMax) and  #the hits to the model should be reasonable 
			   	  ( $gapOnModel >$gapOnModelMin) )){
#			   	$logger->trace( "findRNA_hmmer: checking $i to ". ($i - 1) );
						
				my $r = &merge($ar->[$i], $ar->[$i-1]);
				my $passRes = 0;
				if($passRes = &pass($r)==1){
					push @{ $resultsFull->{$k} }, $r;
					$logger->trace( "The result from this line is $passRes");
				}
			 }
		}
	}
	$logger->trace("------------------------------------------------");
}


# this sub merges a second record into the first.
# it assumes the records are provided in order of aliFrom and simply assigns
# the ali/hmm/end-to coordinates to that of the 2nd record.
# the maximum e-score is used.
# when the two components to be merged are too close
# then one new record is created without keeping the individual ones
sub merge {
	my ($rec0,$rec1)=@_;

	$logger->trace(  "------------------ Merging:---------------------");
	$logger->trace(	"rec0\t". join(',',@$rec0[0..22]));
	my $msg="   envcoords: ";
	foreach my $ar( @{$rec0->[24]} ){$msg.=   "$ar->[0], $ar->[1] :";}
	$msg.=   "\thmmcoords: ";
	foreach my $ar( @{$rec0->[23]} ){$msg.=    "$ar->[0], $ar->[1] :";}
	$msg.=  "  with:";
	$logger->trace($msg);
	$logger->trace("rec1\t" . join(',',@$rec1[0..22]));
	
	$msg=   "   envcoords: ";
	foreach my $ar( @{$rec1->[24]} ){$msg.=    "$ar->[0], $ar->[1] :";}
	$msg.=   "\thmmcoords: ";
	foreach my $ar( @{$rec1->[23]} ){$msg.=    "$ar->[0], $ar->[1] :";}
	$logger->trace($msg);
	
	# keep track of the coordinates of each exon
	if(!defined($rec0->[25])  ){
		$rec0->[25]=[];
		$logger->trace( "=============  adding the initial coordinates for the exons ", $rec0->[19]," _ ", $rec0->[20] , "\n");
		push @{$rec0->[25]}, [$rec0->[19],$rec0->[20] ];
	}
	$logger->trace( "=============  adding the adi/nal coordinates for the exons ", $rec1->[19]," _ ", $rec1->[20] , "\n");
	
	
	my $distanceBetweenExons=findGap([ $rec0->[19], $rec0->[20] ], [ $rec1->[19], $rec1->[20] ]);
	
	# hmmfrom, hmmto = hmmfrom1, hmmto2
	if($rec0->[5] eq '+'){
		$rec0->[13]=$rec1->[13];
	}else{
		$rec0->[14]=$rec1->[14];
	}
	$logger->trace( "rec0 hmmfrom becomes ", $rec0->[13],"\n");
	# alifrom, alito = alifrom1, alito2
	$rec0->[16]=$rec1->[16];
	$logger->trace( "rec0 alifrom becomes ", $rec0->[16],"\n");
	# envfrom, envto = envfrom1, endto2
	$rec0->[19]=$rec1->[19];
	$logger->trace( "rec0 envrom becomes ", $rec0->[19],"\n");
	# score = score1 + score2
	$rec0->[9] += $rec1->[9];
	# c_eval = min(c_eva1,c_eval2)
	if ($rec0->[11] =~ /^\d+\.?\d*e\-(\d+)$/) {
		my $e1=$1;
		if ($rec1->[11] =~ /^\d+\.?\d*e\-(\d+)$/ and $1<$e1) {
			$rec0->[11]=$rec1->[11];
		}
	}
	# i_eval = min(i_eva1,i_eval2)
	if ($rec0->[12] =~ /^\d+\.?\d*e\-(\d+)$/) {
		my $e1=$1;
		if ($rec1->[12] =~ /^\d+\.?\d*e\-(\d+)$/ and $1<$e1) {
			$rec0->[12]=$rec1->[12];
		}
	}
	# add coordinates
	$logger->trace( "rec0 is now ", $rec0->[13],"-",$rec0->[14]," ", $rec0->[16],"-",$rec0->[17]," ", $rec0->[19],"-",$rec0->[20]);

	push @{$rec0->[25]}, [ $rec1->[19], $rec1->[20] ];
	
#	if($distanceBetweenExons< $minExonDistance){
#		$logger->trace("** The distance between the exons is $distanceBetweenExons**");
#	}
	#hmm coordinates
	foreach my $ar (@{$rec1->[23]}) {
		push @{$rec0->[23]}, $ar;
	}
	#env coordinates
	foreach my $ar (@{$rec1->[24]}) {
		
		if(findGap( [  $rec0->[24]->[0]->[0] , $rec0->[24]->[0]->[1] ], [ $ar->[0], $ar->[1] ]) < $minExonDistance  ){
			$logger->trace("** The distance between the exons is $distanceBetweenExons**");
			( $rec0->[24]->[0]->[0], $rec0->[24]->[0]->[1] )=
				mergeCoordinates(  [  $rec0->[24]->[0]->[0] , $rec0->[24]->[0]->[1] ], [ $ar->[0], $ar->[1] ] );
		}else{
			push @{$rec0->[24]}, $ar;
		}
	}

	
	
	
	return $rec0;
}

# parse the hmm output
# and return a hash of results that will 
# contain all the hits;


sub runHMMPrediction{
	my ($cmdline, $contigsFn, $sequenceSize, $resultsFull,$resultsPartial)=@_;
	my $parser='off';
	my $description;
	my $modelName;
	my $modelSize;
	my $strand;
	my $results={};
	my $sequenceName;

	require IPC::Open2;
	my $inStream=Bio::SeqIO->new('-file'=>$contigsFn, '-format' => 'Fasta');
	my $orientation=0; # indicates the orientation of the sequence as it will be searched by hmmsearch
			# if orientation is 0 the script reads a sequence from the file 
			# 1 indicates processing the forward strand
			# -1 indicates processing the revcom sequence
	my $seq;
	for(;;){
		if($orientation ==0){
			$seq=$inStream->next_seq();
			$orientation =1;
		}
		if(!defined($seq)){last;}
		#print "passing sequence ", $seq->display_id(),"\n";

		my $pid= open2(  my $Reader, my $Writer, "$cmdline " );
		if (!defined($pid)){die "Cannot execute $cmdline\n";}
		my $sequenceString="";
		if($orientation ==1){
			$sequenceString=CommonFunc::wrapSeq($seq->seq(), 60);
		}else{
			$sequenceString=CommonFunc::wrapSeq( $seq->revcom->seq(), 60 );
		}
		print $Writer ">".$seq->display_id()."\n". $sequenceString ."\n";
		close $Writer;
		for(;;){
			my $cm;
			if($cmdline=~/hmmsearch/){  $cm=hmmsearchParser->new( "-hmmOutStream"=>$Reader) ; }
			if(!defined($cm)){last;}
			while (my $hit= $cm->nextHit() ){
	#			$hit->printHit();
				my $rec=[];
				my $hmmcoordinates=[];	
				@$hmmcoordinates=([ $hit->{MdlFrom}, $hit->{MdlTo}]);
				my $envcoordinates=[];
				my ($envfrom,$envto, $seqfrom, $seqto);
				if($orientation ==1){
					($envfrom,$envto, $seqfrom, $seqto)=
						($hit->{EnvFrom},$hit->{EnvTo},$hit->{SeqFrom},$hit->{SeqTo});
					
				}else{	
					($envfrom,$envto, $seqfrom, $seqto)=
						($sequenceSize->{  $hit->{targetName} }-$hit->{EnvTo}+1,
						$sequenceSize->{  $hit->{targetName} }-$hit->{EnvFrom}+1,
						$sequenceSize->{  $hit->{targetName} }-$hit->{SeqTo}+1,
						$sequenceSize->{  $hit->{targetName} }-$hit->{SeqFrom}+1);
				}
				@$envcoordinates=([ $envfrom, $envto]); #start and end on the envelope
				my $strand;
				if($orientation ==1){$strand ='+';}
				else{ $strand = '-';}
				@$rec=( $hit->{targetName},
						$sequenceSize->{  $hit->{targetName} },
						$cm->{queryDescription},
				$cm->{queryName},
				$cm->{querySize},
				$strand,
				" ",
				$hit->{Rank},
				$hit->{Qual},
				$hit->{Score},
				$hit->{Bias},
				$hit->{cEvalue},
				$hit->{iEvalue},
				$hit->{MdlFrom},
				$hit->{MdlTo},
				$hit->{MdlBoundaries},
				$seqfrom,
				$seqto,
				$hit->{SeqBoundaries},
				$envfrom,
				$envto,
				$hit->{EnvBoundaries},
				$hit->{Acc}, 
				$hmmcoordinates, 
				$envcoordinates);
		# save record to evaluate later (indexed by contig->alignment start position)
				my $passRes=0;
				if($passRes=pass($rec)==1){
					push @{ $resultsFull->{$hit->{targetName}} }, $rec;
					$logger->trace( "loadFile:  Loading to full list.\n") ;
				}
				else{
					push @{ $resultsPartial->{$hit->{targetName}}}, $rec;
					$logger->trace( "loadFile:  Loading to partial list.\n") ;
				}
			}
		}
		close $Reader  or die "Something bad happened with the $cmdline after it started\n";;
		waitpid($pid, 0);
		if($orientation ==1){$orientation =-1;}
		else{$orientation =0;}
	} # end of processing $contigFn
	print "End of processing all input sequences \n";    
}




###############
sub runCMPrediction{
	my ($cmdline, $contigsFn, $sequenceSize, $resultsFull,$resultsPartial)=@_;
	my $parser='off';
	my $description;
	my $modelName;
	my $modelSize;
	my $strand;
	my $results={};
	my $sequenceName;

	
	open (my $Reader, "$cmdline|") or die "Command failed: $cmdline $!";

	for(;;){
		my $cm;
		if($cmdline=~/cmsearch/){  $cm=cmsearchParser->new( "-cmOutStream"=>$Reader) ; }
		if(!defined($cm)){last;}
		while (my $hit= $cm->nextHit() ){
#			$hit->printHit();
			my $rec=[];
			my $hmmcoordinates=[];	
			@$hmmcoordinates=([ $hit->{MdlFrom}, $hit->{MdlTo}]);
			my $envcoordinates=[];
			@$envcoordinates=([ $hit->{EnvFrom}, $hit->{EnvTo}]);
			my $strand=$hit->{SeqStrand} ;
			
			@$rec=( $hit->{targetName},
					$sequenceSize->{  $hit->{targetName} },
					$cm->{queryDescription},
			$cm->{queryName},
			$cm->{querySize},
			$strand,
			" ",
			$hit->{Rank},
			$hit->{Qual},
			$hit->{Score},
			$hit->{Bias},
			$hit->{cEvalue},
			$hit->{iEvalue},
			$hit->{MdlFrom},
			$hit->{MdlTo},
			$hit->{MdlBoundaries},
			$hit->{SeqFrom},
			$hit->{SeqTo},
			$hit->{SeqBoundaries},
			$hit->{EnvFrom},
			$hit->{EnvTo},
			$hit->{EnvBoundaries},
			$hit->{Acc}, 
			$hmmcoordinates, 
			$envcoordinates);
            # save record to evaluate later (indexed by contig->alignment start position)
            
			my $passRes=0;
			if($passRes=pass($rec)==1){
				push @{ $resultsFull->{$hit->{targetName}} }, $rec;
				$logger->trace( "runCMPrediction:  Loading to full list.\n") ;
			}
			else{
				push @{ $resultsPartial->{$hit->{targetName}}}, $rec;
				$logger->trace( "runCMPrediction:  Loading to partial list.\n") ;
			}
			
		}
	}

	close $Reader  ;#or die "Something bad happened with the $cmdline after it started\n";;
}



# find the gap between two pairs of coordinates
sub findGap{
	my($firstPair, $secondPair)=@_;
#	$logger->trace("findGap: comparing ", $firstPair->[0],"-",$firstPair->[1],"  to  ", 
#										  $secondPair->[0],"-",$secondPair->[1]);
#	
	if($firstPair->[0] > $firstPair->[1]){
		# swap them
		my $t= $firstPair->[0];
		$firstPair->[0]= $firstPair->[1];
		$firstPair->[1]=$t;
	}
	if($secondPair->[0] > $secondPair->[1]){
		# swap them
		my $t= $secondPair->[0];
		$secondPair->[0]= $secondPair->[1];
		$secondPair->[1]=$t;
	}
	
	if($firstPair->[0] > $secondPair->[0]){
		#swap them
		my $t= $firstPair;
		$firstPair= $secondPair;
		$secondPair= $t;	
	}
	$logger->trace("findGap: comparing ", $firstPair->[0],"-",$firstPair->[1],"  to  ", 
										  $secondPair->[0],"-",$secondPair->[1]);
	return ($secondPair->[0] - $firstPair->[1]);
	
}


# give the merged coordinates for a set of coordinates
sub mergeCoordinates{
	my($firstPair, $secondPair)=@_;
	
	if($firstPair->[0] > $firstPair->[1]){
		# swap them
		my $t= $firstPair->[0];
		$firstPair->[0]= $firstPair->[1];
		$firstPair->[1]=$t;
	}
	if($secondPair->[0] > $secondPair->[1]){
		# swap them
		my $t= $secondPair->[0];
		$secondPair->[0]= $secondPair->[1];
		$secondPair->[1]=$t;
	}
	if($firstPair->[0] > $secondPair->[0]){
		#swap them
		my $t= $firstPair;
		$firstPair= $secondPair;
		$secondPair= $t;	
	}
	$logger->trace("mergeCoordinates: merging ", $firstPair->[0],"-",$firstPair->[1],"  to  ", 
										  $secondPair->[0],"-",$secondPair->[1]);
	return ( $firstPair->[0], $secondPair->[1]);
}



# check that a hit, either a simple hit or a composite i.e. hit with introns
# is good enough to be considered valid
sub pass {
   		my ($rec)=@_;
   		
        my( $sequenceName,$sequenceSize,$description,
            $modelName,$modelSize,$strand,$sp,$no,$test,$score,$bias,$c_eval,$i_eval,
            $hmmfrom,$hmmto,$dom1,$alifrom,$alito,$dom2,$envfrom,$envto,$dom3,$acc,$hmmcoordinates,$envcoordinates)=@$rec;
        my $line=join(' ', $sp,$no,$test,$score,$bias,$c_eval,$i_eval,$hmmfrom,$hmmto,$dom1,$alifrom,$alito,$dom2,$envfrom,
            $envto,$dom3,$acc);
            
           
        my $modelFlex=10; # nucleotides from the ends of the models that we allow missing in order to accept the gene.
        my $hmmPerc=0.80; # percentage of model that is aligned of hte sequence to be considered full
        my $hitPerc=0.80; # percentage of hit that is aligned to be considered full; 
        my $modelFlex2=10; # nucleotides from the ends of the sequence that we allow missing
        my $fullGene='no';
		my $fullModel='no';
        if($modelName =~/ssu/){
        	$modelFlex=30; $modelFlex2=40;
        	$hitPerc=0.8;
        }
        if($modelName =~/lsu/){
        	$modelFlex=50; $modelFlex2=60;
#        	$hitPerc=0.8;
        }
        if($modelName =~/tsu/){
        	$modelFlex=20; $modelFlex2=20;
			$hitPerc=0.8;
        }
        my $sizeOnSequence= abs($envto-$envfrom)+1;
        my $sizeOnModel = abs($hmmto-$hmmfrom) +1;
        $logger->trace("------------------- pass -----------------------\n");
        $logger->trace("Model $modelName size: $modelSize\n");
        $logger->trace("$sequenceName ($sequenceSize)\t$line\n");
        
        if(scalar(@$envcoordinates)>1){
        	$logger->trace( "We are checking a split location \n");
        	$sizeOnSequence =0;
        	$envfrom=$envcoordinates->[0][0];
        	$envto=$envcoordinates->[-1][1];
        	foreach my $coordPair( @$envcoordinates ){
        		$logger->trace("$coordPair->[0] - $coordPair->[1]    (envelope coordinates)");
        		if($coordPair->[0] < $envfrom){ $envfrom=$coordPair->[0];}
        		if($coordPair->[1] > $envto){ $envto=$coordPair->[1];}
        		$sizeOnSequence+= $coordPair->[1] - $coordPair->[0]+1;
        	}
        	$logger->trace( "The sum size of the envelope is $sizeOnSequence\n");
        }
        if(scalar(@$hmmcoordinates)>1){
        	$logger->trace( "We are checking a split location \n");
        	$hmmfrom=$hmmcoordinates->[0][0];
        	$hmmto=$hmmcoordinates->[-1][1];
        	$sizeOnModel =0;
        	foreach my $coordPair( @$hmmcoordinates ){
        		$logger->trace( "$coordPair->[0] - $coordPair->[1]    (coordinates on model)");
        		if($coordPair->[0] < $hmmfrom){ $hmmfrom=$coordPair->[0];}
        		if($coordPair->[1] > $hmmto){ $hmmto=$coordPair->[1];}
        		$sizeOnModel+= $coordPair->[1] - $coordPair->[0]+1;
        	}
        	$logger->trace( "The sum size of the model is $sizeOnModel\n");
        }
        my $returnValue=0;
        # Tests to decide if the line is worth processing
        
        #check if the ends of the hits on the model are good to retain
        my ($keepModelTo, $keepModelFrom)=('n','n');
        if($hmmfrom < $modelFlex or $dom1 =~/^\[/ ){ $keepModelFrom = 'y';}
        if($hmmto> $modelSize - $modelFlex or $dom1 =~/\]$/ ){ $keepModelTo='y';}
#        $logger->trace("For model we keep start $keepModelFrom, we keep end $keepModelTo\n");
        if($keepModelFrom eq 'y'){ $logger->trace("Model 5' is OK. Threshold for rejecting is $modelFlex.");}
        else{ $logger->trace("Model 5' is TRUNCATED. Threshold is $modelFlex.");}
        if($keepModelTo eq 'y'){ $logger->trace("Model 3' is OK. Threshold for rejecting is $modelFlex.");}
        else{ $logger->trace("Model 3' is TRUNCATED. Threshold is $modelFlex.");}
        
        
        #check if the hit on the sequence is close to the boundaries of the sequence (i.e. partial hit)
   
        my ($keepHitTo, $keepHitFrom)=($keepModelTo, $keepModelFrom);
        if($keepModelTo eq 'n' and ($envfrom < $modelFlex2 ) ){ $keepHitFrom = 'y';}
        if($keepModelFrom eq 'n' and ($envto> $sequenceSize - $modelFlex2 )  ){ $keepHitTo='y';}
        # regardless of the above if hmmsearch has identified the hits as full size (highly unlikely to activate this code but just in case)
        if($dom3 =~/^\[/ ){ $keepHitFrom = 'y';}
        if($dom3 =~/\]$/) {$keepHitTo='y';}
#        $logger->trace("For hit we keep start $keepHitFrom, we keep end $keepHitTo\n");
        if($keepHitFrom eq 'y'){ $logger->trace("Hit 5' is OK. Threshold is $modelFlex2. ");}
        else{ $logger->trace("Hit 5' is TRUNCATED. Threshold is $modelFlex2.(used for partial predictions)");}
        if($keepHitTo eq 'y'){ $logger->trace("Hit 3' is OK. Threshold is $modelFlex2. ");}
        else{ $logger->trace("Hit 3' is TRUNCATED. Threshold is $modelFlex2.(used for partial predictions)");}
        
        # check if the size of the hit (env from to env to ) is long enough
	# Kostas Aug 25 2012
	# This filter was removed since we identified that there were some rRNAs
	# that had a small exon (<50 nt) close to the 3' of the gene. This filter	
	# was allowing the remaining gene (i.e. all the other exons ) to be considered
	# long enough.
	# Since all the predictions are expected to be very close to the ends of the model
	# a percentage filter is obsolete.
#         my $sizeOnSequence=abs( $envto - $envfrom ) +1;


	# check the size of the model and predicted genes
	# this is to ensure that large chunks of the model are not missing inside the model
	# when the boundaries of the model align very well

	# the gene size is appropriate
    if( $sizeOnSequence > $hitPerc * $modelSize){
       	$fullGene='yes';
    }else{ $fullGene ='no';}
    $logger->trace("Since the size of the hit is $sizeOnSequence and the limit is ", $hitPerc * $modelSize," the flag for full gene is set to $fullGene\n");
	# the model hit is appropriate        
	if( $sizeOnModel > $hmmPerc * $modelSize){
		$fullModel='yes';
	}else{
		$fullModel='no';
	}
	$logger->trace("Since the size of the model is $sizeOnModel and the limit is ",  $hmmPerc * $modelSize," the flag for full model is set to $fullModel\n");

	# check the boundaries
	

	# the quality of the prediction is of poor quality
#        if ($test eq '?') {
#        	$logger->trace("found ?. This means it is not a significant result\n");
#            $returnValue= 0; # not significant result
#        }
#        elsif($keepModelFrom eq 'y' and $keepModelTo eq 'y') {
#        	$logger->trace("we found a hit to the model from start to finish\n");
#            $returnValue= 1; # we have a full hit to the model
#		}
#        elsif($keepHitFrom eq 'y' and $keepHitTo eq 'y'){
##        	if($modelName !~/tsu/ or ($fullGene eq 'yes' and $modelName=~/tsu/) ){
##			if($fullGene eq 'yes'){
#			$logger->trace("the hit is covering the full size of the sequence\n");
#        		$returnValue=1;
##        	}
#         }
#        elsif( $keepModelTo eq 'y' and $keepHitFrom eq 'y'){
#        	$logger->trace( "we have exceeded the start of the sequence and the end of the model is correct\n");
#        	$returnValue=1;
#        } elsif( $keepModelFrom eq 'y' and $keepHitTo eq 'y'){
#        	$logger->trace("we have exceeded the end of the sequence and the start of the model is correct\n");
#        	$returnValue=1;
#        } 


	# final test based on sequence size
	if($partialFlag eq 'on'){
			# We don't care about the size as long as the boundaries of the query (on the genome) have been marked as 'y' for both sides
			$logger->trace("This run looks for partial and full size genes on the genome");
			$logger->trace("For partial genes we do allow genes or models that are not full length");
			if($keepHitFrom eq 'n' or $keepHitTo eq 'n'){
				$returnValue=0;
				$logger->trace("But in this case the hit on the sequence is not good enough");
			}else{
				if($keepModelFrom eq 'n'){ $rec->[19] = "<". $rec->[19];}
				if($keepModelTo   eq 'n'){ $rec->[20] = ">". $rec->[20];}
				$returnValue=1;
			}

	}else{
		$logger->trace("This run looks only for full size genes on the genome");
		# we need to verify that both the model and the sequence are complete and that the size of the gene is correct
		if($fullGene eq 'no' or $fullModel eq 'no' or 
			$keepModelFrom eq 'n' or $keepModelTo eq 'n' or
			$keepHitFrom eq 'n' or $keepHitFrom eq 'n'){
			$returnValue=0;
			$logger->trace("For full only genes we don't allow genes or models that are not full length");
		}else{
			$returnValue=1;
		}
	}
	

	# Kostas Aug 25 2012
	# This filter was removed since we identified that there were some rRNAs
	# that had a small exon (<50 nt) close to the 3' of the gene. This filter	
	# was allowing the remaining gene (i.e. all the other exons ) to be considered
	# long enough.
	# Since all the predictions are expected to be very close to the ends of the model
	# a percentage filter is obsolete.
# 	elsif( abs($sizeOnModel) > $hmmPerc*$modelSize and abs($sizeOnSequence) >$seqPerc*$modelSize) {
#         	$logger->trace( "we still have the full hit to the model $sizeOnModel>",$hmmPerc*$modelSize," $sizeOnSequence>", $seqPerc*$modelSize);
#             $returnValue= 1; # we still have the full model
#         } 
# 	elsif(  abs($sizeOnSequence) >$seqPerc*$modelSize) {
#         	$logger->trace( "we still have a hit with envelope comparable to the model size ", $seqPerc*$modelSize);
#             $returnValue= 1; # we still have the full model
#         } 

    $logger->trace( "The size of the envelope is $sizeOnSequence. The model is $modelSize and the coverage of the model $sizeOnModel");
	if($returnValue==1){
        	$logger->trace("Return will be 'full gene'\n");
	}else{
		$logger->trace("Return will be 'partial gene'\n");
	}
    # the quality of the prediction is of poor quality
    if ($test eq '?') {
      	$logger->trace("found ?. This means it is not a significant result and this will superseed all other results.\n");
        $returnValue= 0; # not significant result
    }
        
		$logger->trace("------------------------------------------------");
		return $returnValue;
}	
1;
