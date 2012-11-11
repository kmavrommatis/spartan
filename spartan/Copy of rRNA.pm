package rRNA;

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/..";
use lib "$FindBin::Bin/../../helpingHands";
use geneCallerVersions;
use CommonFunc;

use IPC::Open2;

use Bio::SeqIO;


sub findRRNA{
	my ($contigFn, $domains, $type, $dataArray) = @_;
	
	my $counter = 0;

	my $removedID={};
	foreach my $domain (@$domains) {
		my $count = &findRNA_hmmer($contigFn, $domain, $type, $dataArray);
		$counter += $count;
	}
	
	print "\nHmmer found $counter rRNAs ($type)!\n\n";
}


sub findRNA_hmmer{
	my ($contigsFn, $domain, $type, $dataArray, $hmmsearchOpt) = @_;
	$domain=lc($domain);
	my $programVersion = geneCallerVersions::getVersion("hmmer");
	print "Predicting rRNA genes with HMMs and $programVersion\n";
	my $evalue_cutoff=0.01;
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
	if($type eq 'tsu'){$specificFlag='max';}
	else{$specificFlag='nobias';}
	my $cmdline = $ENV{ hmmsearchBin }." -E $evalue_cutoff --noali --$specificFlag -Z 1 --cpu 0 ".
					$ENV{ rrnaHmmsDir }."$hmmFn $contigsFn";
	print "Executing command $cmdline\n";
	&loadFile( $cmdline , $sequenceSize ,$resultsFull, $resultsPartial);
	#resultsFull have all the records that have full name
	# resultsPartial have all the records that have partial hits
	# 	we need to check if there are any exons in the partial list

	&processPartial($resultsPartial,$resultsFull,$type);
	
	#process the full models
	my $coordinates=[];
	my $gene_id=1;
	foreach my $k(keys (%$resultsFull) ) {
		my $ar=$resultsFull->{$k};
		foreach my $r( @$ar) {
			my $description= $r->[2];
			my $geneIdentifier=$r->[0]."_".uc($domain)."_rh". substr($type,0,1) . $gene_id;
			my $groupData="ID=$geneIdentifier; Version=$programVersion; ".
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
	my $intronSize=	1500;
	my $gapOnModelMax=10;
	my $gapOnModelMin=0;
	if($type =~/ssu/){$gapOnModelMax=30;$gapOnModelMin=-300;}
	if($type =~/lsu/){$gapOnModelMax=50;$gapOnModelMin=-500;}
	
	
	#process the partial hits
	foreach my $k(keys (%$resultsPartial) ){
		print "Processing hits from sequence $k\n" if $ENV{ verbose } eq "Yes";
		my $ar=$resultsPartial->{$k}; 
		
		my $coordinates=[];
		@$ar=sort{ $a->[16] <=> $b->[16]} @$ar;
		#find successive hits to the same model, that have distance shorter than a determined value
		#			$sequenceName,$sequenceSize->{$sequenceName},$description,$modelName,$modelSize,$strand,
		#							   $sp,$no,$test,$score,$bias,$c_eval,$i_eval,$hmmfrom,$hmmto,$dom1,
		#                			   $alifrom,$alito,$dom2,$envfrom,$envto,$dom3,$acc, [coordinates]
		for(my $i=1;$i<scalar(@$ar);$i++){
			if ($ENV{ verbose } eq "Yes"){
				print ">>> Intron size " ,$ar->[$i]->[16]," - ", $ar->[$i-1]->[17] ,"=", $ar->[$i]->[16] - $ar->[$i-1]->[17] ,"\n";
				print ">>> Gap on model " ,$ar->[$i]->[13]," - ", $ar->[$i-1]->[14] ," =", $ar->[$i]->[13]- $ar->[$i-1]->[14] , "\n";
			}
			if (  ($ar->[$i]->[3] eq $ar->[$i-1]->[3]) and   #the hits are to the same model
				  (($ar->[$i]->[16] - $ar->[$i-1]->[17] < $intronSize) and   # the size of the intron should be reasonable 
				  ($ar->[$i]->[16] - $ar->[$i-1]->[17] >0)) and
				  (($ar->[$i]->[13]- $ar->[$i-1]->[14] <$gapOnModelMax) and  #the hits to the model should be reasonable 
			   	  ($ar->[$i]->[13]- $ar->[$i-1]->[14] >$gapOnModelMin) )){
				if ($ENV{ verbose } eq "Yes"){
			   		print "findRNA_hmmer: checking $i to ". ($i - 1) ." \n" ;
						
				}
				my $r = &merge($ar->[$i], $ar->[$i-1]);
				my $passRes = 0;
				if($passRes = &pass($r)==1){
					push @{ $resultsFull->{$k} }, $r;
					print "The result from this line is $passRes\n" if $ENV{ verbose } eq "Yes";
				}
			 }
		}
	}
}


# this private sub merges a second record into the first.
# it assumes the records are provided in order of aliFrom and simply assigns
# the ali/hmm/end-to coordinates to that of the 2nd record.
# the maximum e-score is used.
sub merge {
	my ($rec0,$rec1)=@_;

	if( $ENV{ verbose } eq "Yes" ) {
		print  "Merging:\n";
		print "rec0\t",join(',',@$rec0[0..22]);
		print  "\nenvcoords: ";
		foreach my $ar( @{$rec0->[24]} ){print  "$ar->[0], $ar->[1] :";}
		print  "\thmmcoords: ";
		foreach my $ar( @{$rec0->[23]} ){print  "$ar->[0], $ar->[1] :";}
		print "  with\nrec1\t",join(',',@$rec1[0..22]);
		print  "\nenvcoords: ";
		foreach my $ar( @{$rec1->[24]} ){print  "$ar->[0], $ar->[1] :";}
		print  "\thmmcoords: ";
		foreach my $ar( @{$rec1->[23]} ){print  "$ar->[0], $ar->[1] :";}
		print  "\n";
	}
	
	# keep track of the coordinates of each exon
	if(!defined($rec0->[25])  ){
		$rec0->[25]=[];
		print "=============  adding the initial coordinates for the exons ", $rec0->[19]," _ ", $rec0->[20] , "\n";
		push @{$rec0->[25]}, [$rec0->[19],$rec0->[20] ];
	}
		print "=============  adding the adi/nal coordinates for the exons ", $rec1->[19]," _ ", $rec1->[20] , "\n";
	push @{$rec0->[25]}, [ $rec1->[19], $rec1->[20] ];
	# hmmfrom, hmmto = hmmfrom1, hmmto2
	$rec0->[13]=$rec1->[13];
	print "rec0 hmmfrom becomes ", $rec0->[13],"\n";
	# alifrom, alito = alifrom1, alito2
	$rec0->[16]=$rec1->[16];
	print "rec0 alifrom becomes ", $rec0->[16],"\n";
	# envfrom, envto = envfrom1, endto2
	$rec0->[19]=$rec1->[19];
	print "rec0 envrom becomes ", $rec0->[18],"\n";
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
	print "rec0 is now ", $rec0->[13],"-",$rec0->[14]," ", $rec0->[16],"-",$rec0->[17]," ", $rec0->[19],"-",$rec0->[20],"\n";
	#hmm coordinates
	foreach my $ar (@{$rec1->[23]}) {
		push @{$rec0->[23]}, $ar;
	}
	#env coordinates
	foreach my $ar (@{$rec1->[24]}) {
		push @{$rec0->[24]}, $ar;
	}
	
	
	return $rec0;
}

# parse the hmm output
# and return a hash of results that will 
# contain all the hits;


sub loadFile{
	my ($cmdline, $sequenceSize, $resultsFull,$resultsPartial)=@_;
	my $parser='off';
	my $description;
	my $modelName;
	my $modelSize;
	my $strand;
	my $results={};
	my $sequenceName;
	open (HMMER, "$cmdline|") or die "Command failed: $cmdline $!";
    while (my $line=<HMMER>) {
		chomp $line;
        if (!$line or $line=~/^\s*#/ or $line=~/^\s+$/ or $line=~/-----/) {
            next;
        } 
        elsif ($line=~/Internal pipeline/) { # end of hmmer output , show some statistics
            $parser='off';
        }
        elsif ($line=~/Description:\s*(.+)/) {
            $description=$1; 
		} 
		elsif($line=~/Query:\s+(\S+)\s+\[M=(\d+)\]/){
            $modelName=$1;
            $modelSize=$2;
		    $strand= $modelName =~/rc$/ ? '-':'+';
        } 
        elsif ($line=~/^>> (\S*)/) {
            $sequenceName = $1;
            $parser='on';
        } 
        elsif ($parser eq 'on') {
#            print "Processing line:$sequenceName=> $line\n" ;

            my @row=split(/\s+/,$line);
            if (@row != 17) {
                warn("Incomplete record: $line\n");
                next;
            }
            my($sp,$no,$test,$score,$bias,$c_eval,$i_eval,$hmmfrom,$hmmto,$dom1,
                $alifrom,$alito,$dom2,$envfrom,$envto,$dom3,$acc)=@row;

#            my $score1000=int(($alito-$alifrom+1)/$modelSize*1000+0.5);
#            $score1000=1000 if $score1000 > 1000;
#            print "score1000= $score1000\n";
            my $hmmcoordinates=[];
            @$hmmcoordinates=([ $hmmfrom, $hmmto]);
            my $envcoordinates=[];
            @$envcoordinates=([ $envfrom, $envto]);
            my $rec=[];
            @$rec=( $sequenceName,$sequenceSize->{$sequenceName},$description,
            			$modelName,$modelSize,$strand,$sp,$no,$test,$score,$bias,
            			$c_eval,$i_eval,$hmmfrom,$hmmto,$dom1,$alifrom,$alito,$dom2,
            			$envfrom,$envto,$dom3,$acc, $hmmcoordinates, $envcoordinates);
                # modelName is 3
                # E-score is 9
                # ali is 16-17
                # hmm is 13-14
                # env is 19-20
            # save record to evaluate later (indexed by contig->alignment start position)
            my $passRes=0;
            if($passRes=pass($rec)==1){
            	push @{ $resultsFull->{$sequenceName} }, $rec;
            	print "loadFile:  Loading to full list.\n" if $ENV{ verbose } eq "Yes";
            }
            else{
            	push @{ $resultsPartial->{$sequenceName}}, $rec;
            	print "loadFile:  Loading to partial list.\n" if $ENV{ verbose } eq "Yes";
            }
        }
    }
    close HMMER;
}



sub pass {
   		my ($rec)=@_;
   		
   		my $verbose;
   		if($ENV{ verbose } eq "Yes"){$verbose=1;}
   		else{$verbose=0;}
   		
        my( $sequenceName,$sequenceSize,$description,
            $modelName,$modelSize,$strand,$sp,$no,$test,$score,$bias,$c_eval,$i_eval,
            $hmmfrom,$hmmto,$dom1,$alifrom,$alito,$dom2,$envfrom,$envto,$dom3,$acc,$hmmcoordinates,$envcoordinates)=@$rec;
        my $line=join(' ', $sp,$no,$test,$score,$bias,$c_eval,$i_eval,$hmmfrom,$hmmto,$dom1,$alifrom,$alito,$dom2,$envfrom,
            $envto,$dom3,$acc);
            
           
        my $modelFlex=10; # nucleotides from the ends of the models that we allow missing in order to accept the gene.
        my $hmmPerc=0.80; # percentage of model that is aligned of hte sequence
        my $seqPerc=0.80; # percentage of sequence that is aligned on the model
        my $modelFlex2=10; # nucleotides from the ends of the sequence that we allow missing
        if($modelName =~/ssu/){
        	$modelFlex=20; $modelFlex2=20;
        }
        if($modelName =~/lsu/){
        	$modelFlex=50; $modelFlex2=20;
        }
        
        my $sizeOnSequence= abs($envto-$envfrom)+1;
        my $sizeOnModel = abs($hmmto-$hmmfrom) +1;
        if($verbose>0){print "------------------- pass -----------------------\n";}
        if($verbose>0){print "Model $modelName size: $modelSize\n";}
        if($verbose>0){print "$sequenceName ($sequenceSize)\t$line\n";}
        
        if(scalar(@$envcoordinates)>1){
        	print "We are checking a split location \n" if $verbose>0;
        	$sizeOnSequence =0;
        	$envfrom=$envcoordinates->[0][0];
        	$envto=$envcoordinates->[-1][1];
        	foreach my $coordPair( @$envcoordinates ){
        		print "$coordPair->[0] - $coordPair->[1]\n";
        		if($coordPair->[0] < $envfrom){ $envfrom=$coordPair->[0];}
        		if($coordPair->[1] > $envto){ $envto=$coordPair->[1];}
        		$sizeOnSequence+= $coordPair->[1] - $coordPair->[0]+1;
        	}
        	print "The size of the envelope is $sizeOnSequence\n";
        }
        if(scalar(@$hmmcoordinates)>1){
        	print "We are checking a split location \n" if $verbose>0;
        	$hmmfrom=$hmmcoordinates->[0][0];
        	$hmmto=$hmmcoordinates->[-1][1];
        	$sizeOnModel =0;
        	foreach my $coordPair( @$hmmcoordinates ){
        		print "$coordPair->[0] - $coordPair->[1]\n";
        		if($coordPair->[0] < $hmmfrom){ $hmmfrom=$coordPair->[0];}
        		if($coordPair->[1] > $hmmto){ $hmmto=$coordPair->[1];}
        		$sizeOnModel+= $coordPair->[1] - $coordPair->[0]+1;
        	}
        	print "The size of the model is $sizeOnModel\n"if $verbose>0;
        }
        my $returnValue=0;
        # Tests to decide if the line is worth processing
        
        #check if the ends of the hits on the model are good to retain
        my ($keepModelTo, $keepModelFrom)=('n','n');
        if($hmmfrom < $modelFlex or $dom1 =~/^\[/ ){ $keepModelFrom = 'y';}
        if($hmmto> $modelSize - $modelFlex or $dom1 =~/\]$/ ){ $keepModelTo='y';}
        if($verbose >0){print "for model do we keep start = $keepModelFrom, do we keep end = $keepModelTo\n";}
        
        my ($keepHitTo, $keepHitFrom)=('n','n');
        if($envfrom < $modelFlex2 or $dom3 =~/^\[/ ){ $keepHitFrom = 'y';}
        if($envto> $sequenceSize - $modelFlex2 or $dom3 =~/\]$/ ){ $keepHitTo='y';}
        if($verbose >0){print "for hit do we keep start = $keepHitFrom, do we keep end= $keepHitTo\n";}
        
        if ($test eq '?') {
        	if($verbose>0){print "found ?. This means it is not a significant result\n";}
            $returnValue= 0; # not significant result
        } elsif($keepModelFrom eq 'y' and $keepModelTo eq 'y') {
        	if($verbose>0){print "we found a full hit to the model\n";}
            $returnValue= 1; # we have a full hit to the model
        } elsif($keepHitFrom eq 'y' and $keepHitTo eq 'y'){
        	if($verbose>0){print "the hit is covering the full size of the sequence\n";}
        	$returnValue=1;
        } elsif( $keepModelTo eq 'y' and $keepHitFrom eq 'y'){
        	if($verbose>0){ print "we have exceeded the start of the sequence and the end of the model is correct\n";}
        	$returnValue=1;
        } elsif( $keepModelFrom eq 'y' and $keepHitTo eq 'y'){
        	if($verbose>0){ print "we have exceeded the end of the sequence and the start of the model is correct\n";}
        	$returnValue=1;
#        } elsif( abs($sizeOnModel) > $hmmPerc*$modelSize and abs($sizeOnSequence) >$seqPerc*$modelSize) {
#        	if($verbose>0){print "we still have the full hit to the model $sizeOnModel>",$hmmPerc*$modelSize," $sizeOnSequence>", $seqPerc*$modelSize ,"\n"; }
#            $returnValue= 1; # we still have the full model
        } 

        if($verbose>0){
        	print "The size of the envelope is $sizeOnSequence. The model is $modelSize and the coverage of the model $sizeOnModel\n";
        	print "Return will be $returnValue\n";
        }
        
		if($verbose>0){print "------------------------------------------------\n";}
		return $returnValue;
}	
1;
