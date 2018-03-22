# Ortholog Assembly and Concatenation
# by Benjamin Karin


# This program takes the results of 
# tblastn-probes-to-genomes.sh or blast-probes-to-genomes.sh
# and combines blast hits to the same contig into one contig.
#
# It also extracts the sequences from the genomes, aligns them,
# and trims them to a reference sequence.


############################################################################################

# USER INPUT - Please adjust this section

############################################################################################

# Dependencies
library(ape)
# Must have samtools istalled on computer. 
# I use homebrew to install everything on mac.

# Directory holding probefile and TaxonNames file (one directory up) and TaxonNames.txt file
# Program will also will save a spreadsheet of the results here
wd <- "~/Desktop/bird_phylogenomics_data"
setwd(wd)

# Set probefile name, must be in working directory
probefile <- "Gallus-Exons.fas"
#probefile <- "relec-additions-chicken.fas"
#probefile <- "AHE-galGal-seqs.fas"

# Read in TaxonNames file, must be in working directory
Taxa <- read.table("TaxonNames.txt",header=F)

# Set Genome Directory and file extension
genomes.path <- paste0(wd,"/assembly")
genome.ext <- ".fa"

# Set BlastResults Directory
blastresults.dir <- paste0(wd,"/BlastResults-RELEC-Gallus/")
#blastresults.dir <- paste0(wd,"/BlastResults-AHE-Gallus/")

# Output spreadsheet name for condensed blast hits
outtable <- "Blast.Hits.Results-RELEC-additions-10mar2017.csv"

# Flanking sequence to extract (Default is 500 bp).
# This will be trimmed later to reference (probefile seq)
flank <- 600




############################################################################################

# RUN THE REST OF THE PROGRAM FROM HERE:

############################################################################################
## CONDENSE SEEGMENTED BLAST HITS

# Get names of queries by samtools indexing probe file to create name file (by creating query_names.txt)
system(paste0("samtools faidx ",probefile)) # make index file
probe_fai <- read.table(paste0(probefile,".fai")) # read in index file
query_names <- probe_fai[1] # read in names from index file
no.queries <- nrow(query_names) # Find number of probes from probefile

# Setup output dataframe
results <- matrix(nrow=no.queries,ncol=nrow(Taxa))
results <- as.data.frame(results)
colnames(results) <- as.vector(Taxa$V1)
rownames(results) <- as.vector(query_names$V1)

for(j in 1:nrow(Taxa)) { # Loop through the number of taxa
  current.taxon <- as.vector(Taxa[j,1])
  current.taxon.dir <- paste(blastresults.dir,current.taxon,sep="")
  if (dir.exists(paste0(current.taxon.dir,"/condensed"))) { } 
  else { dir.create(paste0(current.taxon.dir,"/condensed")) } # make a directory for the condensed blast files if it doesn't already exist
  for (i in 1:no.queries) { # Loop through the number of genes
    query_current <- query_names[i,1]
    infile <- paste(current.taxon.dir,"/",query_current,".txt",sep="")
    if (file.info(infile)$size == 0) {
      print(paste(current.taxon,query_current, "had no BLAST hit!"))
      results[i,j] <- "0"
    } else {
      blast.result.i <- read.table(file=infile, header=F)
      if (nrow(blast.result.i) == 1 ) {
        write.table(blast.result.i, file=paste(current.taxon.dir,"/condensed/",query_current,".txt",sep=""), row.names=F,col.names=F, quote=F, sep="\t")
        print(paste(current.taxon,query_current,"had just 1 hit. Moved it as is."))
        results[i,j] <- " "
      } else if (nrow(blast.result.i) > 1 ) {
        if (length(levels(blast.result.i[,2]))==1) {
          pos.vals <- c(blast.result.i[,9],blast.result.i[,10])
          pos1 <- min(pos.vals)
          pos2 <- max(pos.vals)
          blast.result.i.condense <- blast.result.i[1,]
          blast.result.i.condense[,9] <- pos1
          blast.result.i.condense[,10] <- pos2
          write.table(blast.result.i.condense, file=paste(current.taxon.dir,"/condensed/",query_current,".txt",sep=""), row.names=F,col.names=F, quote=F, sep="\t")
          print(paste(current.taxon,query_current,"was condensed"))
          results[i,j] <- "  "
        } else if (nrow(blast.result.i[blast.result.i[,2]==blast.result.i[1,2],]) > 1) { #Take the first hit and condense the positions in that contig it there were more than one hit on a contig.
          level1 <- levels(blast.result.i[,2])[1]
          blast.result.i.level1 <- blast.result.i[blast.result.i[,2]==blast.result.i[1,2],]
          pos.vals <- c(blast.result.i.level1[,9],blast.result.i.level1[,10])
          pos1 <- min(pos.vals)
          pos2 <- max(pos.vals)
          blast.result.i.condense <- blast.result.i.level1[1,]
          blast.result.i.condense[,9] <- pos1
          blast.result.i.condense[,10] <- pos2
          write.table(blast.result.i.condense, file=paste(current.taxon.dir,"/condensed/",query_current,".txt",sep=""), row.names=F,col.names=F, quote=F, sep="\t")
          results[i,j] <- paste(nrow(blast.result.i),"hits,",nrow(blast.result.i.level1),"condensed",sep=" ")
        } else { # Just take the first one
          write.table(blast.result.i[1,], file=paste(current.taxon.dir,"/condensed/",query_current,".txt",sep=""), row.names=F,col.names=F, quote=F, sep="\t")
          print(paste(current.taxon,query_current,"has hits on multiple contigs. Took first hit"))
          results[i,j] <- paste(nrow(blast.result.i), "hits, used 1st")
        } 
      } else print("what happened here?")
    }
  }
}

write.csv(results, file=outtable)


############################################################################################
## EXTRACT AND ADD FLANKING SEQUENCE

# make the fasta extraction directory
out.d <- paste0(blastresults.dir,"_fasta_genome_extract")
dir.create(out.d)
#file.remove(paste(out.d,"/*",sep="")) #remove any old output file

# make single fasta output files for each locus with the original probe as reference sequence at the top.
reference.seqs <- read.FASTA(probefile) # For Reference Alignments (reference as first sequence)
for (i in 1:no.queries) {
  query_current <- as.vector(query_names[i,1])
  #system(paste("touch ",out.d,"/",query_names[i,1],".fas",sep="")) # FOR BLANK FILES (no reference)
  write.dna(reference.seqs[query_current], file = paste(out.d,"/",query_names[i,1],".fas",sep=""),format="fasta",colsep="")
  system(paste("sed -i '' 's/>",query_current,".*/>",query_current,"__reference/g' ",out.d,"/",query_current,".fas",sep=""))      
}

# Extract the sequences
for(j in 1:nrow(Taxa)) {
  current.taxon = as.vector(Taxa[j,1])
  current.taxon.dir <- paste0(blastresults.dir,current.taxon)
  current.faidx <- read.table(paste0(genomes.path,"/",current.taxon,genome.ext,".fai"),header=F)
  for (i in 1:no.queries) { #the number of genes
    query_current <- as.vector(query_names[i,1])
    infile <- paste0(current.taxon.dir,"/condensed/",query_current,".txt")
    if (file.exists(infile)) {
      blast.result.i <- read.table(file=infile, header=F)
    if (file.info(infile)$size == 0) {
      print(paste(current.taxon,query_current, "had no BLAST hit!"))
    } else {
      #Extract positions from condensed files
      pos1 <- blast.result.i[1,9]
      pos2 <- blast.result.i[1,10]
      contig.name <- as.vector(blast.result.i[1,2]) # NEEDS TO COME FROM THE BLAST FILE
      faidx.row <- current.faidx[current.faidx[,1]==contig.name,]
      contig.length <- faidx.row[1,2]
      
      #if position 1 is greater than pos 2, set 1 as max and 2 as min, and vice versa.
      if (pos1 > pos2) {  
        max.pos <- pos1
        min.pos <- pos2
        TO.REVERSE <- TRUE
      } else if (pos1 < pos2) {
        max.pos <- pos2
        min.pos <- pos1
        TO.REVERSE <- FALSE
      } else {
        print("step 1: cannot read positions")
      }
      #max.extract <- max.pos  # FOR NO FLANK
      #min.extract <- min.pos
      #If min position is less than flank, set min=1, if greater than flank, subtract flanking region
      if (min.pos <= flank) {
        print(paste(current.taxon,query_current,"min position close to contig end"))
        min.extract <- 1
      } else if (min.pos > flank) {
        min.extract <- min.pos - flank
      } else {
        print("step 2: cannot read positions")
      }
      ##If length - max position is less than flank, then set max=length, otherwise add flanking region
      leftover.contig <- contig.length - max.pos
      if (leftover.contig <= flank) {
        print(paste(current.taxon,query_current,"max position close to contig end"))
        max.extract <- contig.length
      } else if (leftover.contig > flank) {
        max.extract <- max.pos + flank
        print(paste(current.taxon,query_current,"successfully extracted"))
      } else {
        print("step 3: cannot read positions")
      }
  
      # Extract the contig using the new positions using samtools faidx into a temporary fasta file
      system(paste0("samtools faidx ",genomes.path,"/",current.taxon,genome.ext," ",contig.name,":",min.extract,"-",max.extract," > ",out.d,"/current-aln.fasta")) 
      if (TO.REVERSE==T) {
        seq <- read.dna(paste0(out.d,"/current-aln.fasta"),format="fasta",as.character=F)
        write.dna(complement(seq),file=paste0(out.d,"/current-aln.fasta"),format="fasta",colsep="")
      }
      system(paste0("sed -e 's/>.*/>",query_current,"__",current.taxon,"/g' ",out.d,"/current-aln.fasta > ",out.d,"/current-aln_edit.fasta"))
      system(paste0("cat ",out.d,"/current-aln_edit.fasta >> ",out.d,"/",query_current,".fas"))
      system(paste0("rm ",out.d,"/current*"))
    }
    } else {
      print(paste(current.taxon,query_current,"BLAST file does not exist"))
    }
  }
}



############################################################################################
## Align the single gene fastas

out.d2 <- paste0(blastresults.dir,"_alignments") # set a new output directory
dir.create(out.d2) # make the fasta extraction directory
out.d <- paste0(blastresults.dir,"_fasta_genome_extract") #should already be set from previous step

for (i in 1:no.queries) { #the number of genes
  locus <- query_names[i,1]
  system(paste0("mafft --auto --thread 4 --adjustdirectionaccurately ",out.d,"/",locus,".fas > ",out.d2,"/",locus,".fas"))
}





  
  

