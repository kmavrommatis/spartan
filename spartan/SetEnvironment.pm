package SetEnvironment;

sub setEnv {
	
		$ENV{ PIPELINE_VERSION } = "4.0.0";

		$ENV{ hmmsearchBin } = "/usr/local/bin/hmmsearch";
		if ( ! -e $ENV{ hmmsearchBin } ) { die "The file ".$ENV{ hmmsearchBin }." doesn't exist!" };
		$ENV{ rrnaHmmsDir } = "/Users/kmavrommatis/Documents/Projects/spartan/hmms/";
		if ( ! -e $ENV{ rrnaHmmsDir } ) { die "The file ".$ENV{ rrnaHmmsDir }." doesn't exist!" };
#		$ENV{ trnaDb } = $ENV{ PIPELINE_PATH }."/data/tRNA_db/tRNA.NR.fna";
#		if ( ! -e $ENV{ trnaDb } ) { die "The file ".$ENV{ trnaDb }." doesn't exist!" };
		$ENV{ TRNASCAN_DIR } ="/usr/local/trnascan/";
		if ( ! -e $ENV{ TRNASCAN_DIR } ) { die "The file ".$ENV{ TRNASCAN_DIR }." doesn't exist!" };
		$ENV{ trnascanBin } =$ENV{TRNASCAN_DIR}."/tRNAscan-SE";
		$ENV{ cmsearchBin11 } = "/usr/local/bin/cmsearch";
		if ( ! -e $ENV{ cmsearchBin11 } ) { die "The file ".$ENV{ cmsearchBin11 }." doesn't exist!" };
		$ENV{ TRNAmodel } = "/Users/kmavrommatis/Documents/Projects/spartan/hmms//trna.cm";
		if ( ! -e $ENV{ TRNAmodel } ) { die "The file ".$ENV{ TRNAmodel }." doesn't exist!" };
		$ENV{ AragornBin } ="/usr/local/bin/aragorn";
		if ( ! -e $ENV{ AragornBin } ) { die "The file ".$ENV{ AragornBin }." doesn't exist!" };
		$ENV{blastallBin}="/sw/bin/blastall";
		$ENV{ grepBin }="/usr/bin/grep";
}
1;
