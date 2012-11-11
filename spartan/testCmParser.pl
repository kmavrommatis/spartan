#!/jgi/tools/bin/perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::RealBin";
use cmsearchParser;


my $inFn='/house/homedirs/k/kmavromm/house/Tests/rna/rna.cmout';

open (my $fh, $inFn);
my $cm=cmsearchParser->new( "-cmOutStream"=>$fh);

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

while (my $hit= $cm->nextHit() ){
	$hit->printHit();
	print $hit->{Rank}, "\t", $hit->{Qual},"\t", $hit->{ Mdl} ."\n";
	print $cm->{querySize},"\n";
#	my ($aa,$res)= findAntiCodon( $hit);
#	print "$hit->{SeqFrom}, $hit->{SeqTo} is $aa / $res\n";
}

close($fh);




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
		"CCC" => "Gly" ,"TCA" => "Sec(p)" 
	);
	
	my $i=index( uc($trnaSeq),"GACUUAAAAUC" );
	if($i<0){return ;}
	my $antiCodon= substr( $genomeSeq, $i+4, 3);
	$antiCodon=~s/U/T/g;
	my $res=$residues{ $antiCodon };
#	print "The trna is $trnaSeq,\nthe genome is $genomeSeq \nthe anticodon is $antiCodon (Position $i)\n";
	return ($antiCodon,$res);
}