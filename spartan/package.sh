#!/bin/bash



package_dir=$1;
mkdir -p $package_dir


# copy the scripts/modules 
for file in findrRNA.pl findrRNAs.pm cmsearchParser.pm hmmsearchParser.pm findrRNA.config rRNA.pm ; do
cp $file $package_dir/$file
done



cd ../
for file in geneCallerVersions.pm; do
	cp $file $package_dir/$file;
done

cd ../
for file in SetParametersPipeline.pm ; do
	cp $file $package_dir/$file;
done

cd helpingHands
for file in CommonFunc.pm ; do
	cp $file $package_dir/$file;
done

cd ../
mkdir -p $package_dir/hmms
cp geneCalling/rnaAnnotation/hmms/*.hmm $package_dir/hmms
cp geneCalling/rnaAnnotation/hmms/*.cm $package_dir/hmms


echo 'To install the program copy this directory in the location you want to use it e.g. $LOCATION
and modify the line \
$ENV{ rrnaHmmsDir } = $ENV{ PIPELINE_PATH }."/data/rrna_hmms/"
in the file SetParametersPipeline.pm to
$ENV{ rrnaHmmsDir } = "$LOCATION/hmms"
' > $package_dir/README
