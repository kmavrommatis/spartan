package hmmsearchParser;

use strict;
use warnings;



sub new{
	my ($class, @arguments) = @_;
	my $self={};
    
   
   	for (my $i=0;$i<scalar(@arguments);$i++){
		if(	
		    $arguments[$i] eq '-hmmOutStream'
		){
			my $a= substr( $arguments[$i] , 1 );
			$self->{ $a  }=$arguments[++$i];
		}
	} 
    if(!defined($self->{hmmOutStream})){
    	die "Please provide the output file from hmmsearch";
    }
    if(eof( $self->{hmmOutStream}  )){ return undef;} # return undef if we are at teh end of the file
    bless($self, $class);
    $self->_init();
    return $self;
}


sub DESTROY{
	my ($self)=@_;
}



sub _init{
	my($self)=@_;
	my $rfh=$self->{hmmOutStream};
	my $line=<$rfh>;
	if( $line !~ /# hmmsearch/){ # we are not at the beginning of a new model hit. We need to reach that point
#		print "We are not at the beginning of the file, we need to go further down\n";
		while($line=<$rfh>){
			chomp $line;
			if($line eq "//"){ last ; }
		}
	}
	while(my $line=<$rfh>){
		chomp $line;
#		print "_init: working with line $line\n";
		
		
		
		if( $line=~/^#/){ $self->parseComment( $line ) };
		if($line=~/^Query:\s+(\S+)\s+\[M=(\d+)\]/){
			$self->{queryName }=$1;
			$self->{querySize }=$2;
		}

		if($line=~/Description:\s+(.+)/){
			$self->{queryDescription}=$1;
		}
		if($line=~/^Domain annotation /){
			last;
		}
	}
}



{
	my $currentName="";
sub nextHit{
	my($self)=@_;
    my $rfh=$self->{hmmOutStream};
    
    while(my $line=<$rfh>){
    	if($line=~/Internal pipeline statistics/){return undef;}
    	if($line=~/No targets detected that satisfy reporting thresholds/){return undef;}
#    	print "Processing line $line";
    	if(eof($rfh) ){return undef;}	
    	if( $line =~/^>> (\S+)/ ){
    		$currentName=$1;
    		$line=<$rfh>; 
    		if ($line =~/score  bias/){ 
#    			print "Skipping crap\n";
    			$line=<$rfh>;
    			
    			$line=<$rfh>;
    		}else{
    			$currentName=undef;
    		}
#    		print "Passing line $line\n";
    		
    	}
		if( defined($currentName) and 
    			$line ne "\n" ){
#    		print "Processing line for $currentName\n";
    		my $hit=hmmHit->new( $rfh , $line, $currentName);
    		if(!$hit->{targetName}){
#    			print "Hit does not have targetName\n";
    			$currentName=undef;
    			return undef;}
			return $hit;
    	}
    }
   

}
}
# parse the comment lines and get metadata about the program and the run
sub parseComment{
	my($self, $line)=@_;
	
	if($line =~/# (HMMER .+)/){
		$self->{programVersion}=$1;
	}
	if($line =~/# query HMM file:\s+(\S+)/){
		$self->{queryHMMFile}=$1;
	}
	if($line =~/# target sequence database:\s+(\S+)/){
		$self->{targetSequenceFile}=$1;
	}		
}


1;

package hmmHit;

use strict;
use warnings;

sub new{
	
	my ($class, @arguments) = @_;
	my $self={};
    bless($self, $class);
   	
    $self->{rfh}=$arguments[0];
    $self->{targetName}=$arguments[2];
    $self->nextHit($arguments[1]);
    
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
	$self->{ Rank } . "\t" .
	$self->{ Qual } . "\t" .
	$self->{ cEvalue } . "\t" .
	$self->{ iEvalue } . "\t".
	$self->{ Score} . "\t" .
	$self->{ Bias } . "\t" .
	$self->{ Mdl  } . "\t" .
	$self->{ MdlFrom} . "\t" .
	$self->{ MdlTo  } . "\t" .
	$self->{ MdlBoundaries  } . "\t" .
	$self->{ SeqFrom}  . "\t" .
	$self->{ SeqTo} . "\t" .
	$self->{ SeqBoundaries} . "\t" .
	$self->{ EnvFrom}  . "\t" .
	$self->{ EnvTo} . "\t" .
	$self->{ EnvBoundaries} . "\t" .
	
	$self->{ Acc} . "\n" ;
}


sub nextHit{
	my($self,$line)=@_;
	my $rfh=$self->{rfh};

	
#	print "nextHit: processing line $line\n";
	chomp $line;
	if($line=~/\s+(\d+)\s+([!?])\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
#		print "nextHit: parsing line\n";
		$self->{ Rank }=$1;
		$self->{ Qual }=$2;
		$self->{ Score }=$3;
		$self->{ Bias}=$4;
		$self->{ cEvalue}=$5;
		$self->{ iEvalue}=$6;
		$self->{ Eval }= $self->{iEvalue}; # for consistency with cm parser
		$self->{Mdl}="hmm";
		if($7 < $8){ ($self->{ MdlFrom},$self->{ MdlTo  }) = ($7,$8) ;}
		       else{ ($self->{MdlTo} , $self->{MdlFrom} ) = ($7,$8);}
		$self->{ MdlBoundaries  }=$9;
		if($10<$11){ ($self->{ SeqFrom},$self->{ SeqTo}) = ($10,$11)  ;}
			else   { ($self->{ SeqTo},$self->{ SeqFrom}) = ($10,$11)  ;}
		$self->{SeqStrand}='?'; # for consistency with cm parser
		$self->{ SeqBoundaries  }=$12;
		
		if($13<$14){ ($self->{ EnvFrom},$self->{ EnvTo}) = ($13,$14)  ;}
			else   { ($self->{ EnvTo},$self->{ EnvFrom}) = ($13,$14)  ;}
		$self->{ EnvBoundaries  }=$15;
		$self->{ Acc}=$16;
		if($self->{SeqBoundaries} =~/[\[\]]/){$self->{Trunc}='yes';} # for consistency with cm parser
		else{$self->{Trunc}='no';} # for consistency with cm parser
		$self->{GC}="50%"; # for consistency with cm parser
	}
	if($line eq ""){ $self->{targetName}=undef;return;}
	if($line =~/No targets detected/){ $self->{targetName}=undef;return;}
	
}
1;

