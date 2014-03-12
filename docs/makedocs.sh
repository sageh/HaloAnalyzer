#!/bin/sh

for i in ../AMIGA/*.pm 
do
echo "Processing $i..."
pod2html --infile $i --outfile AMIGA/`basename $i .pm`.html --htmldir AMIGA --htmlroot ../
done

for i in ../*.pod 
do
echo "Processing $i..."
pod2html --infile $i --outfile `basename $i .pod`.html --htmldir ./ --htmlroot ./
done

