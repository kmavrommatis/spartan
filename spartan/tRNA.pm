package tRNA;

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin";
use geneCallerVersions;
use CommonFunc;
use Log::Log4perl;
use IPC::Open2;
use Log::Log4perl;
use Bio::SeqIO;


my $logger =  Log::Log4perl->get_logger('tRNAm');
sub findTRNA{
	my ($contigsFn, $domains, $dataArray, $trnaScanFlag, 
				$trnaBlastFlag, $aragornFlag,$infernalFlag, $transTable, $circular) = @_;
	
	if(!defined($trnaBlastFlag)) { $trnaBlastFlag='off'; }
	if(!defined($trnaScanFlag)) { $trnaScanFlag='on'; }
	if(!defined($aragornFlag)) { $aragornFlag='off'; }
	if(!defined($infernalFlag)) { $infernalFlag='on'; }
	$logger->info( "Processing with tRNAscan $trnaScanFlag.\n");
	$logger->info(  "Processing with BLAST $trnaBlastFlag.\n");
	$logger->info(  "Processing with Aragorn $aragornFlag.\n");
	$logger->info(  "Processing with Infernal $infernalFlag.\n");
	my $counter = 0;
	foreach my $domain(@$domains){
		$counter += &findTRNA_trnascan($contigsFn, $domain, 
										$dataArray) if $trnaScanFlag ne 'off';
	}
	$logger->info(  "\tRNAscan found $counter tRNAs!");
	
	$counter = &findTRNA_blast($contigsFn, 
									$dataArray) if $trnaBlastFlag ne 'off';
	$logger->info(  "\ttrnaBLAST found $counter tRNAs!");
	
	$counter = &findTRNA_aragorn($contigsFn, $dataArray, $transTable, 
										$circular) if $aragornFlag ne 'off';
	$logger->info(  "\tAragorn found $counter tRNAs!");
	
	$counter = &findTRNA_infernal($contigsFn, $dataArray) if $infernalFlag ne 'off';
	$logger->info(  "\tInfernal found $counter tRNAs!");
}


# use blast to find the tRNAs
sub findTRNA_blast {
	my ($contigFn, $dataArray) = @_;
	
	my $programVersion = geneCallerVersions::getVersion("trnablast");
	$logger->info(  "Predicting tRNA genes with BLAST ($programVersion)\n");
	my $counter = 0;
	my $gene_id = 0;
	#open contig file and retrieve sequences
	# mask everything but the ending 150
	my $inStream=Bio::SeqIO->new(-file=>$contigFn,-format=>'fasta');
	while(my $seqObj=$inStream->next_seq()){
		my $length=$seqObj->length();
		
# 		print "Sequence: ", $seqObj->display_id()," length ", $seqObj->length(),"\n";
		# we need to search 1-150
		# and seqlen-150 to seqlen
		my $segments=[];
		if($seqObj->length() > 300){
			push @$segments ," -L 1,150";
			push @$segments," -L ".($seqObj->length()-150) .",". $seqObj->length();
		}else{push @$segments  , "-L 1,". $seqObj->length();}

		foreach my $s(@$segments){
			my $cmd = $ENV{ blastallBin };
			$cmd .= " -m8 -a1 -p blastn -d $ENV{ trnaDb } ";
			$cmd .= $s;
			$logger->debug(  "Running the blast command $cmd\n" );
			my $pid=open2(my $blastReaderFh, my $blastWriterFh, "$cmd "); 
# 			print "Created child $pid\n";
			if (!defined($pid)){$logger->logdie( "Cannot execute $cmd\n");}
# 			print "Running Blast\n";
			print $blastWriterFh ">". $seqObj->display_id()."\n".
							CommonFunc::wrapSeq( $seqObj->seq())."\n";
			close $blastWriterFh  or $logger->logdie( "Command $cmd failed (write to pipe)\n");;
# 			print "Command send to blast ($pid)\n";
			# get only the first line from the blast results
			my $readLine='no'; # changes to yes after the first line has been read, and avoids reloading hits
			my $slack=10;
			while(my $line=<$blastReaderFh>)
# 			if(defined($line) and $line ne "")
			{
				if ($readLine eq 'yes'){next;}
				chomp $line;
				my ($contig,$trna,$alnpid,$alnLen, $j1,$j2,$contig_start,
								$contig_end, $trna_start,$trna_end,$eval,$bit) = split(/\t\s*/, $line);
				
				$readLine='yes';
				my($toid,$amino,$anticodon,$size)=split("_",$trna);
				my $strand='+';
				if($trna_start>$trna_end) {
					$strand='-';
					my $t=$trna_start;
					$trna_start=$trna_end;
					$trna_end=$t;
				}
				# we need some significant hit
				if($alnpid <85){next;}
				if($alnLen <40){next;}
				# the hit needs to be at the end of the sequence
				if($contig_start >$slack and $contig_end < $seqObj->length()-$slack) {
					$logger->debug( "Reject 1. Hit is in the middle of the sequence\n");
					next;
				}
				# the hit should cover a fragment of the trna 
				# trna:    =========>
				# aln:        |||||||
				# contig:     =====================
				if    ($contig_start< $slack and $trna_end > $size-$slack and $alnLen < $length){}
				# trna:                     =========>
				# aln:                      |||||||
				# contig:     =====================
				elsif ($contig_end > $length-$slack and $trna_start < $slack and $alnLen <$length){}
				# trna:             =========>
				# aln:                |||||
				# contig:             =====
				elsif ($contig_start < $slack and $contig_end > $length-$slack ){}
				else{next;}
				my $geneIdentifier = $contig."_tb".$gene_id;
				my $groupData = "ID=$geneIdentifier; product=tRNA_${amino}_${anticodon}; ".
								"RNA_Class_ID=tRNA_${amino}_${anticodon}; ".
								"Name=tRNA_${amino}_${anticodon}; Version=$programVersion; ".
								"tRNA_type=$amino; anticodon=$anticodon";
				push @$dataArray, [$contig, 
									$programVersion,
									"tRNA",
									$contig_start, 
									$contig_end, 
									$bit, 
									$strand,
									".",
									$groupData
								   ];
				$counter++;
				$gene_id++;
	# 				print "Line accepted\n";
			}
			close $blastReaderFh or $logger->logdie( "Something bad happened with the $cmd after it started\n");
			
			waitpid($pid,0);
# 			print "Child $pid finished\n";
		}

	}
	
	return $counter;
}



#use tRNAScan to find tRNA

sub findTRNA_trnascan{
	my ($contigFn, $domain , $dataArray )=@_;
	my $programVersion = geneCallerVersions::getVersion("trnascan");
	$logger->info( "Predicting tRNA genes with $programVersion\n");
	$domain=lc($domain);
	my $gene_id=0;
	my $options = "-q";           # output without headers
	if    ( $domain eq "bacteria"   or $domain eq 'b' ) {$options .= " -B";} 
	elsif ( $domain eq "archaea"    or $domain eq 'a' ) {$options .= " -A";} 
	elsif ( $domain eq "eukaryota"  or $domain eq 'e' ) {$options .= " -E";}	
	elsif ( $domain eq "metagenome" or $domain eq 'm' ) {$options .= " -BA";}
	elsif ( $domain eq "unknown"    or $domain eq 'u' ) {$options .= " -G";}

#	my $logFn = "$temp_dir/trna_search." .$$. "_cove.log";unlink($logFn) if -e $logFn;
#	my $statFn = "$temp_dir/trna_search." .$$. "_cove.stat";unlink($statFn) if -e $statFn;
#	my $outFn=$temp_dir."/trnaScan_out.$$.txt";

	my $cmdline  = $ENV{trnascanBin}." $options";
	
#	$cmdline .= " -l $logFn -m $statFn";
	$cmdline .= " $contigFn";

 	$logger->debug( "tRNAscan command line will be: $cmdline\n");
#        system( "$cmdline 1>$outFn 2>/dev/null");   # note - no error catching here!!
    my @output=`$cmdline`;
	if ( $? ) { $logger->logdie( "Command failed: $cmdline $!"); }
	my $parse='off';
	
	my $counter=0;
	
	foreach my $line (@output)
	{
		chomp $line;
#		print $line."\n";
		if($line =~/^---/){$parse='on'; next;}
		next if $parse eq 'off';
		
		my ($contig, $aa, $start, $end, $amino, $anticodon, $intronStart,$intronEnd, $score)=
			split ("\t", $line);
		$contig =~ s/ //g;
		$start =~ s/^\s+|\s+$//g; # Remove leading and trailing whitespaces
		$end =~ s/^\s+|\s+$//g; # Remove leading and trailing whitespaces
		
		my ($fstart,$fend,$fintronStart,$fintronEnd,$fstrand)=($start,$end,$intronStart,$intronEnd,"+");
		if($end<$start) { 
			($fstart,$fend,$fintronStart,$fintronEnd,$fstrand)=($end,$start,$intronEnd,$intronStart,"-");
		}
		my $coordinates=[];
		if($intronStart >0){ # we have a broken tRNA gene
			push @$coordinates, [$fstart,  $fintronStart, $fstrand];
			push @$coordinates, [$fintronEnd, $fend, $fstrand];

		}
		else{
			push @$coordinates, [ $fstart, $fend, $fstrand];
		}
		my $geneIdentifier=$contig."_".uc($domain)."_ts".$gene_id;
		my $groupData = "ID=$geneIdentifier; product=tRNA_${amino}_${anticodon}; ".
						"RNA_Class_ID=tRNA_${amino}_${anticodon}; ".
						"Name=tRNA_${amino}_${anticodon}; Version=$programVersion; ".
						"tRNA_type=$amino; anticodon=$anticodon";
		$groupData.= ";intron containing" if scalar(@$coordinates)>1;
		$logger->debug("Adding gene from tRNAsnan: $geneIdentifier, tRNA_${amino}_${anticodon} at $fstart $fend");
		push @$dataArray, [$contig, $programVersion,"tRNA",
			 $fstart, $fend, $score, $fstrand,".",$groupData];
		$counter++;
		
        if (@$coordinates>1) {
            my $exonCounter=1;
            foreach my $c( @$coordinates){
                my $exonId="$geneIdentifier.$exonCounter";
                my $groupData = "Parent=$geneIdentifier; ID=$exonId; ".
                				"product=tRNA_${amino}_${anticodon}; ".
                				"RNA_Class_ID=tRNA_${amino}_${anticodon}; ".
                				"Name=tRNA_${amino}_${anticodon}; Version=$programVersion; ".
                				"tRNA_type=$amino; anticodon=$anticodon";
                push @$dataArray,[$contig, $programVersion,"exon",$c->[0], 
                					$c->[1], $score, $c->[2],".",$groupData];
                $exonCounter++;
            }
        }
		$gene_id++;	
	}
	return $counter;
}


#use tRNAScan to find tRNA

sub findTRNA_aragorn{
	my ($contigFn, $dataArray, $transTable, $circular) = @_;
	
	my $programVersion = geneCallerVersions::getVersion("aragorn");
	$logger->info( "Predicting tRNA genes with $programVersion\n");
	
	require Aragorn;
	
	my $intronSize = 1000;
	
	my $outFile = $contigFn.".aragorn.out";
	
	my $aragornObj = Aragorn->new(	"-inputFile"  => $contigFn,
									"-outputFile" => $outFile
								);
	
	if( defined($circular) ) {
		$aragornObj->setSequenceTopology( "circular" );
	}
	else {
		$aragornObj->setSequenceTopology( "linear" );
	}
	$aragornObj->setIntronSize( $intronSize );
	$aragornObj->setTranslationTable( $transTable ) if defined($transTable);
	$aragornObj->runPrediction();

	# Put results in data array
	open FH, "<", $outFile or $logger->logdie( "Cannot open aragorn output: $!\n");
	my $counter = 0;
	while(<FH>) {
		chomp $_;
		$logger->debug("Adding gene from Aragorn: $_");
		push @$dataArray, [split("\t", $_)];
		$counter++;	
	}
	close(FH);
	unlink($outFile);
	
	return $counter;
}





sub findTRNA_infernal{
	my ($contigsFn, $dataArray) = @_;
	require cmsearchParser;
	my $programVersion = geneCallerVersions::getVersion("cmsearch11");
	$logger->info("Predicting tRNA genes with covariance models and $programVersion\n");
	my $evalue_cutoff=0.01;

	my $covFn=$ENV{TRNAmodel};

	my $rnaCounter = 0;
	my $sequenceName="";
	my $parser='off';
	my $description;
	my $modelSize=0;
    my $modelName;
    my $strand;
	my $sequenceSize = &getSequenceSize($contigsFn);

    my $resultsFull={}; # { alifrom } => \@rec
	my $specificFlag;
	my $cmdline =  $ENV{ cmsearchBin11 }.
				" --cut_ga --notextw -Z 1 --cpu 0 ".
				"$covFn $contigsFn";
	$logger->debug( "Executing command $cmdline\n");
	open (HMMER, "$cmdline|") or die "Command failed: $cmdline $!";

	my $cm=cmsearchParser->new( '-cmOutStream'=>\*HMMER );
	my $gene_id=1;
	while(my $hit=$cm->nextHit()){
		
#		$hit->printHit();
		if($hit->{Qual} eq '?'){next;}
		my ($aa, $res)=findAntiCodon( $hit );
		
		my $description= $cm->{queryDescription} . '_' . $res . '_'. $aa;
		my $geneIdentifier=$hit->{targetName}."_cm" . $gene_id;
		my $groupData="ID=$geneIdentifier; Version=$programVersion; ".
						"RNA_Class_ID=$description; product=$description; Name=$description";
		$gene_id++;
		$logger->debug("Adding gene from Infernal: $geneIdentifier, $hit->{SeqFrom}, $hit->{SeqTo}");
		push @$dataArray, [$hit->{targetName},  # name of the sequence 
							$programVersion,	# name of the method 
							'tRNA' ,		# type of molecule 
							$hit->{SeqFrom}, $hit->{SeqTo}, # coordinates of molecule on the contig
							$hit->{Score},	# score
							$hit->{SeqStrand},	# strand
							".",
							$groupData
						   ];
		$rnaCounter++;
	}
	close(HMMER);

	
	return $rnaCounter;
}

# parse the hmm output
# and return a hash of results that will 
# contain all the hits;



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
sub findAntiCodon{
	my($hit)=@_;
	my $trnaSeq= $hit->{ ModelString };
	my $genomeSeq=$hit->{ TargetString };
	my %residues=(
		"GAA" => "Phe" ,"AAA" => "Phe" ,"TAA" => "Leu" ,"CAA" => "Leu" ,
		"GAG" => "Leu" ,"TAG" => "Leu" ,"AAG" => "Leu" ,"CAG" => "Leu" ,
		"GGA" => "Ser" ,"TGA" => "Ser" ,"AGA" => "Ser" ,"CGA" => "Ser" ,
		"GCT" => "Ser" ,"ACT" => "Ser" ,"GTA" => "Tyr" ,"ATA" => "Tyr" ,
		"GCA" => "Cys" ,"ACA" => "Cys" ,"CCA" => "Trp" ,"GGG" => "Pro" ,
		"TGG" => "Pro" ,"AGG" => "Pro" ,"CGG" => "Pro" ,"GTG" => "His" ,
		"ATG" => "His" ,"TTG" => "Gln" ,"CTG" => "Gln" ,"GCG" => "Arg" ,
		"TCG" => "Arg" ,"ACT" => "Arg" ,"CCG" => "Arg" ,"TCT" => "Arg" ,
		"CCT" => "Arg" ,"GAT" => "Ile" ,"AAT" => "Ile" ,"TAT" => "Ile" ,
		"CAT" => "Met" ,"GGT" => "Thr" ,"TGT" => "Thr" ,"AGT" => "Thr" ,
		"CGT" => "Thr" ,"GTT" => "Asn" ,"ATT" => "Asn" ,"TTT" => "Lys" ,
		"CTT" => "Lys" ,"GAC" => "Val" ,"TAC" => "Val" ,"AAC" => "Val" ,
		"CAC" => "Val" ,"CGG" => "Ala" ,"TGC" => "Ala" ,"AGC" => "Ala" ,
		"CGC" => "Ala" ,"GTC" => "Asp" ,"ATC" => "Asp" ,"TTC" => "Glu" ,
		"CTC" => "Glu" ,"GCC" => "Gly" ,"TCC" => "Gly" ,"ACC" => "Gly" ,
		"CCC" => "Gly" ,"CTA" => "Pyl" ,"TCA" => "Sec(p)", "Undet"=>"???"
	);
	
	my $i=index( uc($trnaSeq),"GACUUAAAAUC" );
	my $antiCodon;
	if($i<0){$antiCodon="Undet";}
	else{$antiCodon= substr( $genomeSeq, $i+4, 3);}
	$antiCodon=~s/U/T/g;
	my $res;
	if($residues{$antiCodon}){$res=$residues{ $antiCodon };}
	else{$res="Pseudo";}
	
#	print "The trna is $trnaSeq,\nthe genome is $genomeSeq \nthe anticodon is $antiCodon (Position $i)\n";
	return ($antiCodon,$res);
}
1;
