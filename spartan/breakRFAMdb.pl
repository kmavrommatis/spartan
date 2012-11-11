#!/jgi/tools/bin/perl
use strict;
use warnings;


# use this script to break down a big Rfam database 
# to files, each file contains the model for one type
# of Rfam entry

# get a list of the RFams that we don't want to compare 
# in a genome.
my %excludeRNA= (  
        				"RF00001"=>1, "RF00002"=>1, "RF00005"=>1, "RF00177"=>1,
						"RF01315"=>1, "RF01316"=>1, "RF01317"=>1, "RF01318"=>1, "RF01319"=>1,
						"RF01320"=>1, "RF01321"=>1, "RF01322"=>1, "RF01323"=>1, "RF01324"=>1,
						"RF01325"=>1, "RF01326"=>1, "RF01327"=>1, "RF01328"=>1, "RF01329"=>1,
						"RF01330"=>1, "RF01331"=>1, "RF01332"=>1, "RF01333"=>1, "RF01334"=>1,
						"RF01335"=>1, "RF01336"=>1, "RF01337"=>1, "RF01338"=>1, "RF01339"=>1,
						"RF01340"=>1, "RF01341"=>1, "RF01342"=>1, "RF01343"=>1, "RF01344"=>1,
						"RF01345"=>1, "RF01346"=>1, "RF01347"=>1, "RF01348"=>1, "RF01349"=>1,
						"RF01350"=>1, "RF01351"=>1, "RF01352"=>1, "RF01353"=>1, "RF01354"=>1,
						"RF01355"=>1, "RF01356"=>1, "RF01357"=>1, "RF01358"=>1, "RF01359"=>1,
						"RF01360"=>1, "RF01361"=>1, "RF01362"=>1, "RF01363"=>1, "RF01364"=>1,
						"RF01365"=>1, "RF01366"=>1, "RF01367"=>1, "RF01368"=>1, "RF01369"=>1,
						"RF01370"=>1, "RF01371"=>1, "RF01372"=>1, "RF01373"=>1, "RF01374"=>1,
						"RF01375"=>1, "RF01376"=>1, "RF01377"=>1, "RF01378"=>1, "RF01379"=>1,
						"RF01959"=>1, "RF01960"=>1
					  );

my $wfh=\*STDOUT;
my $modName="";my $nc; my $ga;
my @a;
while(my $l=<>){
	chomp $l;
	
	if($l=~/^INFERNAL/){
		if(scalar(@a)>0){
			if(!$excludeRNA{$modName}){
#			open(my $wfh, ">".$modName.".cm") or die ("Cannot create $modName.cm\n");
				foreach my $i(@a){
					print $wfh $i."\n";
					if($i =~/^GA  / and !defined($nc)){
						$i=~s/GA  /NC  /;
						print $wfh $i."\n";
					}
				}
			}
#			close $wfh;
			
			$modName=undef;
			$ga=undef; $nc=undef;
			@a=();
		}
	}
	push @a, $l;
	if($l=~/ACC\s+(.+)/){
		$modName=$1;
	}
	if($l=~/^GA\s+(\S+)/){ $ga=$1;}
	if($l=~/^NC\s+(\S+)/){ $nc=$1;}

	
}

#open(my $wfh, ">".$modName.".cm") or die ("Cannot create $modName.cm\n");
if(!$excludeRNA{$modName}){
	foreach my $i(@a){
		print $wfh $i."\n";
		if($i =~/^GA  / and !defined($nc)){
			$i=~s/GA  /NC  /;
			print $wfh $i."\n";
                }

	}
}
#close $wfh;

close $wfh;
