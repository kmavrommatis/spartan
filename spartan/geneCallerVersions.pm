package geneCallerVersions;

use strict;
use warnings;


sub getVersion {
	my ($geneCaller) = @_;
	
	my $version = "";
	
	if ($geneCaller eq 'hmmer') {
		my $cmd = $ENV{ hmmsearchBin }." -h";
		open(my $pipe, "$cmd |") or die "Cannot run $cmd\n";
		while(my $line = <$pipe>) {
			chomp $line;
			if($line !~ /^#/) {next;}
			if($line =~ /(HMMER.+);/) {
				$version=$1;
			}
		}
		close $pipe;
	}
	
	elsif ($geneCaller eq 'aragorn') {
		my $cmd = $ENV{ AragornBin }." -h";
		open(my $pipe, "$cmd |")or die "Cannot run $cmd\n";
		while(my $line = <$pipe>) {
			chomp $line;
			if (lc($line) =~ /^\s*$geneCaller/) {
				my @sA = split /\s+/, $line;
				
				$version = "Aragorn ".$sA[1];
				
				last;
			}
		}
		close $pipe;
	}
	
	elsif ($geneCaller eq 'trnascan') {
		my $cmd = $ENV{trnascanBin}." 2>&1";
		open(my $pipe, "$cmd |")or die "Cannot run $cmd\n";
		while(my $line = <$pipe>) {
			chomp $line;
			if (lc($line) =~ /^\s*$geneCaller/){
				$version = $line;

			}
		}
		close $pipe;
	}
	
	elsif ($geneCaller eq 'trnablast') {
		$geneCaller = 'blastall';
		my $cmd = $ENV{ blastallBin }." 2>&1";
		open(my $pipe, "$cmd |") or die "Cannot run $cmd\n";
		while(my $line = <$pipe>) {
			chomp $line;
			if (lc($line) =~ /^\s*$geneCaller/){
				$version = $line;
				
			}
		}
		close $pipe;
		$version =~ s/\s*arguments://;
	}
	
	elsif ($geneCaller eq 'rfam') {
		my $cmd = $ENV{ cmsearchBin }." -h | ".$ENV{ grepBin }." INFERNAL";
		$version = `$cmd`;
		chomp($version);
		$version = substr($version, 2);
	}
	elsif ($geneCaller eq 'cmsearch11') {
		my $cmd = $ENV{ cmsearchBin11 }." -h | ".$ENV{ grepBin }." INFERNAL";
		$version = `$cmd`;
		chomp($version);
		$version = substr($version, 2);
	}
	elsif ($geneCaller eq 'crt') {
		my $cmd = "$ENV{ javaBin } -Xmx512m -jar $ENV{ crtBin } -version";
		my $version = `$cmd`;
		chomp($version);
	}
	
	elsif ($geneCaller eq 'pilercr') {
		my $cmd = "$ENV{ pilercrBin } 2>&1";
		my @outputArray = `$cmd`;
		foreach my $l (@outputArray) {
			chomp $l;
			if ($l =~ /^pilercr/) { 
				$version = $l;
				last; 
			}
		}
	}
	
	elsif ($geneCaller eq 'multiblastx') {
		$version = "MultiBlastX 1.5";
	}
	
	elsif ($geneCaller eq 'prodigal') {
		my $cmd = "$ENV{ ProdigalBin } -v 2>&1";
		open( my $pipe, "$cmd |") or die "Could not get the version for $geneCaller. Command: $cmd. $!\n";;
		while(my $l = <$pipe>) {
			chomp $l;
			if (lc($l) =~ /$geneCaller/) { $version = $l; }
		}
		close $pipe or die "Something happened when trying to close $cmd\n";
	}
	
	elsif ($geneCaller eq 'metagene') {
		$version = "Metagene Annotator 1.0";
	}
	
	elsif ($geneCaller eq "genemark") {
		my $cmd = "$ENV{ GeneMarkBin }";
		my @outputArray = `$cmd`;
		foreach my $l (@outputArray) {
			chomp $l;
			if (lc($l) =~ /$geneCaller/ or lc($l) =~ /version/) { 
				$version = $l;
				last; 
			}
		}
	}
	
	elsif ($geneCaller eq "fraggenescan") {
		$version = "FragGeneScan 1.16";
	}
	
	
	return $version;
}
1;
