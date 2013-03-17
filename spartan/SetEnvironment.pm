package SetEnvironment;
use FindBin;
use File::Basename;
sub setEnv {
	
		$ENV{ PIPELINE_VERSION } = "4.0.0";
		my $hmmsearchBin=`which hmmsearch`; chomp $hmmsearchBin;
		$ENV{ hmmsearchBin } = $hmmsearchBin;
		if ( ! -e $ENV{ hmmsearchBin } ) { die "The file ".$ENV{ hmmsearchBin }." doesn't exist!" };
		my $cmsearchBin=`which cmsearch`; chomp $cmsearchBin;
		$ENV{ cmsearchBin11 } = $cmsearchBin;
		if ( ! -e $ENV{ cmsearchBin11 } ) { die "The file ".$ENV{ cmsearchBin11 }." doesn't exist!" };

		my $trnascanBin = `which tRNAscan-SE`; chomp $trnascanBin;
		my ($bin, $dir) =fileparse($trnascanBin);
		$ENV{ TRNASCAN_DIR}=$dir;
		if ( ! -d $ENV{ TRNASCAN_DIR } ) { die "The directory ".$ENV{ TRNASCAN_DIR }." doesn't exist!" };
		$ENV{trnascanBin}=$trnascanBin;
		my $aragornBin=`which aragorn`;chomp $aragornBin;
		$ENV{ AragornBin } = $aragornBin;
		if ( ! -e $ENV{ AragornBin } ) { die "The file ".$ENV{ AragornBin }." doesn't exist!" };
		
		my $blastallBin= `which blastall`;chomp $blastallBin;
		$ENV{blastallBin}= $blastallBin;
		if ( ! -e $ENV{ blastallBin } ) { die "The file ".$ENV{ blastallBin }." doesn't exist!" };
	
		my $grepBin= `which grep`; chomp $grepBin;
		$ENV{ grepBin }=$grepBin;
		
		
		$ENV{ rrnaHmmsDir } = "$FindBin::RealBin/hmms/";
		if ( ! -e $ENV{ rrnaHmmsDir } ) { die "The file ".$ENV{ rrnaHmmsDir }." doesn't exist!" };
#		$ENV{ trnaDb } = $ENV{ PIPELINE_PATH }."/data/tRNA_db/tRNA.NR.fna";
#		if ( ! -e $ENV{ trnaDb } ) { die "The file ".$ENV{ trnaDb }." doesn't exist!" };
		$ENV{ TRNAmodel } = "$FindBin::RealBin/hmms//trna.cm";
		if ( ! -e $ENV{ TRNAmodel } ) { die "The file ".$ENV{ TRNAmodel }." doesn't exist!" };
		
		
}
1;
