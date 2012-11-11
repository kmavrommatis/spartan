package cmsearchParser;

use strict;
use warnings;



sub new{
	my ($class, @arguments) = @_;
	my $self={};
    
   
   	for (my $i=0;$i<scalar(@arguments);$i++){
		if(	$arguments[$i] eq '-cmOutFile' or
		    $arguments[$i] eq '-cmOutStream'
		){
			my $a= substr( $arguments[$i] , 1 );
			$self->{ $a  }=$arguments[++$i];
		}
	} 
   
   	
    if(!defined($self->{cmOutStream})){
    	die "Please provide the output file from cmdsearch";
    }
	if(eof( $self->{cmOutStream}  )){ return undef;} # return undef if we are at teh end of the file
    bless($self, $class);
    $self->_init();
    return $self;
	
	
}


sub DESTROY{
	my ($self)=@_;
}



sub _init{
	my($self)=@_;
	my $rfh=$self->{cmOutStream};
	while(my $line=<$rfh>){
		chomp $line;
		if( $line=~/^#/){ $self->parseComment( $line ) };
		if($line=~/^Query:\s+(\S+)\s+\[CLEN=(\d+)\]/){
			$self->{queryName }=$1;
			$self->{querySize }=$2;
		}
		if($line=~/Accession:\s+(\S+)/){
			$self->{queryAccession}=$1;
		}
		if($line=~/Description:\s+(.+)/){
			$self->{queryDescription}=$1;
		}
		if($line=~/Hit alignments:/){
			last;
		}
	}
	
	
}


{
my $hit;
sub nextHit{
	my($self)=@_;

	# if there are introns left return those 
	if(defined($hit) and  $self->nextIntron( $hit ) ){

# 		print "Returning the next intron\n";
		return $hit;
	}

	my $rfh=$self->{cmOutStream};
	my $currentName;
	while(my $line=<$rfh>){
		if($line=~/Internal CM pipeline/){return undef;}
		if(eof($rfh) ){return undef;}	
		if( $line =~/^>> (\S+)/ ){
			$hit=cmHit->new( $rfh , $line);
				
			# check if the targetstring indicates an intron
			
			# Kostas ( Sep 06. After discussion with Jim decided not to look for  introns)
#			$self->identifyIntrons( $hit );
			#if there are introns then return the first intron
			$self->nextIntron( $hit );
			
			return $hit;
			
		}
	}

}
}




# if a hit has several introns this subroutine retunrs one by one the introns
# as separate hits
sub identifyIntrons{
	my($self, $hit)=@_;
	$self->{ExonCoordinates}=[]; # store the exon specific data (start,end and sequence of exon)
	$self->{ExonIndex}=0;
	my $s = $hit->{TargetString};	
	my $m = $hit->{ModelString};
	my($exonStart ,$exonEnd )=($hit->{EnvFrom}, $hit->{EnvTo});
	my($modelStart,$modelEnd)=($hit->{MdlFrom}, $hit->{MdlTo});
	
	if($hit->{SeqStrand} eq '-'){ ($exonEnd,$exonStart)=($hit->{EnvFrom}, $hit->{EnvTo}); } 
	# get the sequence of the exons on the genome sequence
	my @exons= split( /\*\[\s*\d+\]\*/, $hit->{ TargetString } );
	my @intronSize= $hit->{TargetString}=~/\*\[\s*(\d+)\]\*/g ; 
	my @exonsModel= split( /\*\[\s*\d+\]\*/, $hit->{ ModelString } );
	my @intronSizeModel= $hit->{ModelString}=~/\*\[\s*(\d+)\]\*/g ; 
	
# 	print join("\n", @exons),"\n";
# 	print join(" = ", @intronSize) ,"\n";
# 
# 	print join("\n", @exonsModel),"\n";
# 	print join(" = ", @intronSizeModel) ,"\n";
	

	

	# check for false Positive introns
	for(my $i=0; $i<scalar(@exons)-1; $i++){
		# if the $intronSize[$i] and $intronSizeModel[$i] have the same length we are not considering this an intron
		if($intronSize[$i]<  $intronSizeModel[$i] or
			abs( $intronSize[$i] - $intronSizeModel[$i]) < 10 ){
			$exons[$i+1] = $exons[$i] . 'N'x $intronSize[$i] . $exons[$i+1];
			$exonsModel[$i+1] = $exonsModel[$i] . 'N'x $intronSizeModel[$i] . $exonsModel[$i+1];
			splice ( @exons, $i,1);
			splice ( @exonsModel, $i,1);
		}
	}

	# we need to store the exons in the array 
	# push @{$self->{ExonCoordinates}},[ $exonStart,$exonEnd, $s , $modelStart,$modelEnd, $modelBoundaries ];
	for(my $i=0; $i<scalar(@exons); $i++){
		$exons[$i]=~s/[-~\.\*]//g;
		$exonsModel[$i]=~s/[-~\.\*]//g;
		if( $hit->{SeqStrand} eq '+'){
			$exonEnd = $exonStart + length( $exons[$i] ) -1 ;
			if($i==scalar(@exons)-1){ $exonEnd = $hit->{ SeqTo};} 
		}else{
			$exonEnd = $exonStart - length( $exons[$i] ) -1 ;
			if($i==scalar(@exons)-1){ $exonEnd = $hit->{ SeqFrom};} 
		}
		$modelEnd= $modelStart + length( $exonsModel[$i]  ) -1;
		if($i==scalar(@exons)-1){ $modelEnd = $hit->{ MdlTo};} 
		my $modelBoundaries=$hit->{MdlBoundaries};
		if($modelStart > $hit->{MdlFrom}) { substr( $modelBoundaries , 0 ,1 )=".";} 
		if($modelEnd < $hit->{MdlTo}) { substr( $modelBoundaries , 1 ,1 )=".";}
# 		print " We identified an exon at Position $exonStart - $exonEnd \n( $exons[$i] )\n($exonsModel[$i])\n on model $modelStart - $modelEnd\n";
		push @{$self->{ExonCoordinates}},[ $exonStart,$exonEnd, $exons[$i] , $modelStart,$modelEnd, $modelBoundaries ];
		if($hit->{SeqStrand} eq '+' ){
			$exonStart = $exonEnd + $intronSize[$i] if ($i<scalar(@exons)-1);
		}else{
			$exonStart = $exonEnd - $intronSize[$i] if ($i<scalar(@exons)-1);
		}
		$modelStart=$modelEnd + $intronSizeModel[$i] if ($i<scalar(@exons)-1);
	}
}

sub nextIntron{
	my($self, $hit)=@_;
	# don't do anything if there are no introns
	if( ! $self->{ ExonCoordinates } or scalar( @{ $self->{ ExonCoordinates } } )==0 ){ return undef; }

	# flush data  if there are introns but we have returned all of them
	if( scalar( @{ $self->{ ExonCoordinates } })>0  and
	    $self->{ExonIndex} == scalar( @{ $self->{ ExonCoordinates } } ) ){
		$self->{ExonCoordinates}=undef;
		$self->{ExonIndex}=0;
		return undef;
	}

	# else update the data with the exon specific information
	my $totalExons= scalar( @{ $self->{ ExonCoordinates } });
	
	if( $self->{ ExonIndex } < $totalExons ){
		# update the coordinates with the exon 
		
		my($exonStart,$exonEnd, $exonSequence, $modelStart,$modelEnd,$modelBoundaries)= @{$self->{ExonCoordinates}->[ $self->{ExonIndex} ++ ]};
		if($modelStart < $modelEnd){ ($hit->{ MdlFrom},$hit->{ MdlTo  }) = ($modelStart,$modelEnd) ;}
		       else{ ($hit->{MdlTo} , $hit->{MdlFrom} ) = ($modelStart,$modelEnd);}
		if($exonStart<$exonEnd){ ($hit->{ SeqFrom},$hit->{ SeqTo}) = ($exonStart,$exonEnd)  ;}
			else   { ($hit->{ SeqTo},$hit->{ SeqFrom}) = ($exonStart,$exonEnd)  ;}
		$hit->{EnvFrom}=$hit->{SeqFrom};# for consistency with hmm models parser
		$hit->{EnvTo }= $hit->{SeqTo};# for consistency with hmm models parser
		$hit->{MdlBoundaries}=$modelBoundaries;
		return $self->{ExonIndex}; 
	}else{
		return undef;
	}


}



# parse the comment lines and get metadata about the program and the run
sub parseComment{
	my($self, $line)=@_;
	
	if($line =~/# (INFERNAL .+)/){
		$self->{programVersion}=$1;
	}
	if($line =~/# query CM file:\s+(\S+)/){
		$self->{queryCMFile}=$1;
	}
	if($line =~/# target sequence database:\s+(\S+)/){
		$self->{targetSequenceFile}=$1;
	}		
}


1;

package cmHit;

use strict;
use warnings;

sub new{
	
	my ($class, @arguments) = @_;
	my $self={};
    bless($self, $class);
   
   	
    $self->{rfh}=$arguments[0];
	
    
    $self->nextHit($arguments[1]);
#    print "Inside cmHit\n";
    return $self;
}


sub targetName{
	my($self)=@_;
	return $self->{targetName};
}

sub printHit{
	my($self)=@_;
	print STDOUT 
	$self->{ targetName } . "\n" .
	$self->{ StructureString } . "\n" .
	$self->{ ModelString } . "\n" .
	$self->{ TargetString } . "\n" .
		$self->{ Rank } . "\t" .
		$self->{ Qual } . "\t" .
		$self->{ Eval } . "\t" .
		$self->{ Score} . "\t" .
		$self->{ Bias } . "\t" .
		$self->{ Mdl  } . "\t" .
		$self->{ MdlFrom} . "\t" .
		$self->{ MdlTo  } . "\t" .
		$self->{ MdlBoundaries  } . "\t" .
		$self->{ SeqFrom}  . "\t" .
		$self->{ SeqTo} . "\t" .
		$self->{ SeqStrand} . "\t" .
		$self->{ SeqBoundaries} . "\t" .
		$self->{ Acc} . "\t" .
		$self->{ Trunc} . "\t" .
		$self->{ GC} . "\n"; 
}


sub nextHit{
	my($self,$line)=@_;

	my $rfh=$self->{rfh};
	


#	print "Processing line $line";
	if($line=~/^>> (\S+)/){ 
		$self->{ targetName }=$1;
	}
	$line=<$rfh>;$line=<$rfh>; # throw these lines
	$line=<$rfh>;
	
	# parse the stats of the alignment
	if($line=~/\s+\((\d+)\)\s+([!?])\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
#		print "Parsing line $line\n";
		$self->{ Rank }=$1;
		$self->{ Qual }=$2;
		$self->{ Eval }=$3;
		$self->{iEvalue}=$self->{cEvalue}=$self->{Eval}; # for consistency with hmm models parser
		$self->{ Score}=$4;
		$self->{ Bias }=$5;
		$self->{ Mdl  }=$6;
		if($7 < $8){ ($self->{ MdlFrom},$self->{ MdlTo  }) = ($7,$8) ;}
		       else{ ($self->{MdlTo} , $self->{MdlFrom} ) = ($7,$8);}
		$self->{ MdlBoundaries  }=$9;
		if($10<$11){ ($self->{ SeqFrom},$self->{ SeqTo}) = ($10,$11)  ;}
			else   { ($self->{ SeqTo},$self->{ SeqFrom}) = ($10,$11)  ;}
		
		$self->{ SeqStrand}=$12 ; 
		$self->{ SeqBoundaries}=$13;
		$self->{EnvFrom}=$self->{SeqFrom};# for consistency with hmm models parser
		$self->{EnvTo }= $self->{SeqTo};# for consistency with hmm models parser
		$self->{EnvBoundaries}=$self->{SeqBoundaries};# for consistency with hmm models parser
		$self->{ Acc}=$14;
		$self->{ Trunc}=$15;
		$self->{ GC}=$16;
	}
	$line=<$rfh>;$line=<$rfh>;
	$line=<$rfh>;
#	print "Processing line $line";
	if($line=~/\s+(\S+)\s+(\S+)/){
		$self->{ StructureString }=$1;
#		print "Structure String = ". $self->{StructureString}."\n";
	}
	$line=<$rfh>;
#	print "Processing line $line";
	if($line=~/\s+\S+\s+\d+\s+(.+)\s+\d+/){
		my $t=$1;
		$self->{ ModelString }=$1;
#		print "Structure String = ". $self->{ModelString}."\n";
	}
	$line=<$rfh>;
	$line=<$rfh>;
#	print "Processing line $line";
	if($line=~/\s+\S+\s+\d+\s+(.+)\s+\d+/){
		$self->{ TargetString }=$1;
#		print "Structure String = ". $self->{TargetString}."\n";
	}
	$line=<$rfh>;$line=<$rfh>;

}







1;

