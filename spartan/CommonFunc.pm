package CommonFunc;
use Exporter();
@ISA=qw(Exporter);
@EXPORT=qw(
isOverlapping 
swap 
runCmd);
use strict;
use FindBin qw( $RealBin );
use lib $RealBin;
use warnings;
use Bio::Perl;

sub runCmd{
	my ($cmd,$verbose)=@_;
	if(!defined($verbose)){$verbose=0;}
	if($verbose>0){ warn "Executing command: $cmd\n";}
	system($cmd);
	if ($?) {die "command: $cmd failed\n"};
}

#wrap a string to multiple lines of length wrapLen
sub wrapSeq{
   my( $seq, $wrapLen ) = @_;

   if( !defined($wrapLen) or $wrapLen eq "" ) {
      $wrapLen = 50;
   }
   my $i;
   my $s2;
   my $len = length( $seq );
	if ($len == 0){return $seq;}
   for( $i = 0; $i < $len; $i += $wrapLen ) {
       my $s = substr( $seq, $i, $wrapLen );
       $s2 .= $s . "\n";
   }

   chomp $s2;
   return $s2;
}



#swap the values of two variables
sub swap{
	my ($var1,$var2)=@_;
	my $temp=$var2;
	$var2=$var1;
	$var1=$temp;
	return($var1,$var2);
}





# decide if two sets of coordinates overlap
sub isOverlapping{
	my ($s1,$e1,$s2,$e2)=@_;
	my $verbose=1;
	my $overlap=0;
	if ($e1< $s1){ ($s1,$e1)= swap($s1,$e1);}
	if ($e2< $s2){ ($s2,$e2)= swap($s2,$e2);}



	if($s2<= $s1 and $e2>= $e1){ $overlap= $e1-$s1 +1;}
	if($s2>= $s1 and $e2<= $e1){ $overlap= $e2-$s2 +1;}
	if($s1<= $s2 and $e1>= $s2 and $e1<= $e2){$overlap= $e1-$s2+1;}
	if($s1>= $s2 and $s1<= $e2 and $e1>= $e2){$overlap= $e2-$s1+1;}
	
 	print "isOverlapping:  $s1 - $e1 with $s2 - $e2 overlap $overlap\n" if $verbose==1;
	return $overlap;
}



# splits a filename. Returns the path and filename
sub splitFn{
	
	my ($fn)=@_;
	my ($path,$filename,$fnBase,$extension)=splitFileName($fn);
	return($path,$filename,$fnBase,$extension);
}


1;
