#! /bin/sh

# dependencies
# must have samtools and blast command line tools.
# homebrew is the easiest way to install 
# (once you have homebrew installed), use these commands to install:
# brew install samtools
# brew install homebrew/science/blast

# To run, drag the directories/files on top of the shell script in the following order, or write them out and copy paste them in.

#genome-dir=$1
#probe-file=$2
#results-dir=$3
#taxon-names-file=$4

echo ""
if [ -d $1 ]
then echo "Genome directory: $1"
else echo "$1 is not a directory...exiting"; exit; fi
echo ""
if [ -f $2 ]
then echo "Probe fasta file: $2"
else echo "$1 has a problem...exiting"; exit; fi
echo ""
if [ -d $3 ]
then echo "blast results directory already exists...exiting"; exit
else mkdir $3
echo "Blast results directory: $3"; fi
echo ""
if [ -f $4 ]
then echo "Search taxon names file: $4"
else echo "$4 has a problem...exiting"; exit; fi
echo ""

# Index the probe file, and get names of probeseqs
samtools faidx $2
cut -f1 $2.fai > query-names.txt


# Make new directories to store results
for taxon in $(cat "$4"); do # Loopthrough taxon name file 
	mkdir $3/$taxon
	echo "Making output directory for $taxon"
done
echo ""

# Run the blast searches
for q in $(cat query-names.txt); do
	echo "Running searches on $q"
	samtools faidx $2 $q > query-current.txt #extract current query
	for taxon in $(cat "$4"); do # Loop through taxo name file 
		tblastn -query query-current.txt -db $1/$taxon.blastdb -out $3/$taxon/$q.txt -evalue 1e-10 -outfmt 6
		if [ -s $3/$taxon/$q.txt ]; then
			echo "$taxon $q Blasted at E-10"
		else
			tblastn -query query-current.txt -db $1/$taxon.blastdb -out $3/$taxon/$q.txt -evalue 1e-1 -outfmt 6
			if [ -s $3/$taxon/$q.txt ]; then
				echo "$taxon $q Blasted at E-1"
			fi
		fi
		if ! [ -s $3/$taxon/$q.txt ]; then
			echo "$taxon $q could not find a blast hit"
		fi
	done
	rm query-current.txt
done
rm query-names.txt
rm $2.fai



