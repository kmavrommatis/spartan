#!/jgi/tools/bin/perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::RealBin";
use hmmsearchParser;


my $inFn='/house/homedirs/k/kmavromm/house/Tests/rna/16.hmmout';

open (my $fh, $inFn);
#my $cm=hmmsearchParser->new( "-hmmOutStream"=>$fh);

#		$self->{ Rank }=$1;
#		$self->{ Qual }=$2;
#		$self->{ Eval }=$3;
#		$self->{ Score}=$4;
#		$self->{ Bias }=$5;
#		$self->{ Mdl  }=$6;
#		$self->{ MdlFrom}=$7;
#		$self->{ MdlTo  }=$8;
#		$self->{ MdlBoundaries  }=$9;
#		$self->{ SeqFrom} =$10;
#		$self->{ SeqTo}=$11;
#		$self->{ SeqStrand}=$12;
#		$self->{ SeqBoundaries}=$13;
#		$self->{ Acc}=$14;
#		$self->{ Trunc}=$15;
#		$self->{ GC}=$16;

while (my $cm=hmmsearchParser->new( "-hmmOutStream"=>$fh)){

	while (my $hit= $cm->nextHit() ){
		
		
		
		print"
$hit->{targetName}, $cm->{queryDescription},
$cm->{queryName}, $cm->{querySize},
$hit->{Rank}, $hit->{Qual}, $hit->{Score}, $hit->{Bias},
$hit->{cEvalue}, $hit->{iEvalue}, $hit->{MdlFrom}, $hit->{MdlTo},
$hit->{MdlBoundaries}, $hit->{SeqFrom}, $hit->{SeqTo}, $hit->{SeqBoundaries},
$hit->{EnvFrom}, $hit->{EnvTo}, $hit->{EnvBoundaries}, $hit->{Acc}\n";	
	}
}
close($fh);



