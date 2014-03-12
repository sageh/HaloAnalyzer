#!/bin/sh

# If the module documentation directories do not exist, then create them
if [ ! -e AMIGA ]; then
	echo "Creating subdirectory AMIGA...";
	mkdir AMIGA;
fi

if [ ! -e Tools ]; then
	echo "Creating subdirectory Tools...";
	mkdir Tools;
fi

if [ ! -e visualize ]; then
	echo "Creating subdirectory visualize...";
	mkdir visualize;
fi

# Go through the module documentation and the C-code POD documentation and
# create a complete set of html documentation.
for i in ../AMIGA/*.pm 
do
echo "Processing $i..."
pod2html --infile $i --outfile AMIGA/`basename $i .pm`.html --htmldir AMIGA --htmlroot ../
done

for i in ../Tools/*.pm 
do
echo "Processing $i..."
pod2html --infile $i --outfile Tools/`basename $i .pm`.html --htmldir Tools --htmlroot ../
done

for i in ../visualize/*.pod 
do
echo "Processing $i..."
pod2html --infile $i --outfile visualize/`basename $i .pod`.html --htmldir visualize --htmlroot ../
done

for i in ../*.pod 
do
echo "Processing $i..."
pod2html --infile $i --outfile `basename $i .pod`.html --htmldir ./ --htmlroot ./
done

