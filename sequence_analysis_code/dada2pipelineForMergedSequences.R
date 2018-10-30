library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library("optparse")

option_list = list(
  make_option(c("-i", "--inputPath"), type="character", default=NULL, 
              help="Path to directory with unmodified FASTQ files", metavar="character"),
  make_option(c("-s", "--sampleFile"), type="character", default=NULL, 
              help="Path to sample file (list of sample ids)", metavar="character"),
  make_option(c("-w", "--workingDir"), type="character", default=NULL, 
              help="Path to output files", metavar="character"),
  make_option(c("-f", "--fTrunc"), type="numeric", default=400, 
              help="truncation length for forward read", metavar="numeric"),
  make_option(c("-e", "--maxEE"), type="numeric", default=2, 
              help="expected number of errors per read", metavar="numeric"),
  make_option(c("-q", "--truncQ"), type="numeric", default=2, 
              help="remove reads that have quality scores below this theshold", metavar="numeric"),
   make_option(c("-l", "--learn"), type="character", default=NULL, 
              help="Path to file containing list of sample ids used for learning error rates)", metavar="character"),
  make_option(c("-t", "--taxonomy"), type="character", default=NULL, 
              help="Path to working directories (with filterd subdirectory)", metavar="character"),
  make_option(c("-m", "--mergeOverlap"), type="numeric", default=100, 
              help="min. overlap for merging sequences", metavar="numeric")
); 


## Make sure options are properly set
opt_parser = OptionParser(option_list=option_list);
params = parse_args(opt_parser);

if (is.null(params$inputPath)){
  print_help(opt_parser)
  stop("Provide a path to the directory containing demultiplexed FASTQ files", call.=FALSE)
}

if (is.null(params$sampleFile)){
  print_help(opt_parser)
  stop("Provide a sample file", call.=FALSE)
}

if (is.null(params$workingDir)){
  print_help(opt_parser)
  stop("Provide an working directory path", call.=FALSE)
}

if (is.null(params$learn)){
  print_help(opt_parser)
  stop("provide list of sample IDS for learning error rates", call.=FALSE)
}

if (is.null(params$taxonomy)){
  params$taxonomy <- file.path(params$workingDir,'rdp_train_set_14.fa.gz');
  print("did not provide a taxnomy file, we will assume rdp_train_set_14.fa.gz is in the working directory")
}


# define path variable, this is constant
# parse input variables that are constant during script
path <- params$inputPath;

# get the sample files (all samples)
samples <- read.table(file.path(params$sampleFile),header = FALSE,sep=',');
samples <- samples$V1;
fns <- as.character(paste0(samples,".fastq"));
fns <- file.path(path, fns);

learn.samples <- read.table(file.path(params$learn),header = FALSE,sep=',');
learn.samples <- learn.samples$V1;
learn.fns <- as.character(paste0(learn.samples,".fastq"));
learn.fns <- file.path(path, fns);


# Make directory and filenames for the filtered fastqs
filt_path <- paste0(params$workingDir,'/filtered')
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtfns <- file.path(filt_path, paste0(samples, "_filt.fastq.gz"))
# Filter
for(i in seq_along(fns)) {
  fastqFilter(fns[i],filtfns[i],
                    trimLeft=10, truncLen=params$fTrunc, 
                    maxN=0, maxEE=params$maxEE, truncQ=params$truncQ, 
                    compress=TRUE, verbose=TRUE)
}


## LEARN ERROR MODEL ##
learn.filtfns <- file.path(filt_path, paste0(learn.samples, "_filt.fastq.gz"))

# dereplicate the samples
learn.dereps <- derepFastq(learn.filtfns, verbose=TRUE)

#Perform joint sample inference and error rate estimation (takes a few minutes):
learn.dadas <- dada(learn.dereps , err=NULL, selfConsist = TRUE,MAX_CONSIST=20)
saveRDS(learn.dadas,file.path(params$workingDir,"learn_dadas.rds"));


# take out the error matrix and write it to an RDS file
err <- learn.dadas[[1]]$err_out;
saveRDS(err, file.path(params$workingDir,"error_rates.rds"));


## INFERENCE ##
# get sample ID's


dadas <- vector("list", length(samples))
dereps <- vector("list", length(samples))

#names(dadaFs) <- samples
#names(dadaRs) <- samples
#names(derepFs) <- samples
#names(derepRs) <- samples

for(s in 1:length(samples)){
    cat("Processing:", s, "\n")
    dereps[[s]] <- derepFastq(filtfns[[s]]);
    dadas[[s]] <- dada(dereps[[s]], err=err);
}

saveRDS(dadas,file.path(params$workingDir,"dadas.rds"));
saveRDS(dereps,file.path(params$workingDir,"dereps.rds"));



names(dadas) <- as.character(samples)
names(dereps) <- as.character(samples);


seqtab <- makeSequenceTable(dadas);
seqtab2 <- removeBimeraDenovo(seqtab,tableMethod="consensus",verbose=TRUE);
saveRDS(seqtab2,'/Users/joshuagoldford/Dropbox/phd/data/Sequencing_projects/Marine_June2017/dada2_merged/seqtab.rds')

csvfile <- '/Users/joshuagoldford/Dropbox/phd/data/Sequencing_projects/Marine_June2017/dada2_merged/seqtab.csv';
write.table(t(seqtab2),file=csvfile,append=FALSE,quote=FALSE,sep=",",eol = "\n", na = "NA", dec = ".",col.names=NA,row.names=TRUE);


database = list()
database$base <- "~/Dropbox/phd/data/Sequencing_projects/dada2_scripts/assignTaxonomy";
database$train <- file.path(database$base,"silva_nr_v123_train_set.fa.gz");
database$species <- file.path(database$base,"silva_species_assignment_v123.fa.gz");

seqs <- colnames(seqtab2);
tt <- assignTaxonomy(seqs, database$train)
#tt.plus <- addSpecies(tt, database$species,allowMultiple=FALSE, verbose=TRUE)
saveRDS(tt,'/Users/joshuagoldford/Dropbox/phd/data/Sequencing_projects/Marine_June2017/dada2_merged/taxonomy_trimmed.rds');
csvfile <- '/Users/joshuagoldford/Dropbox/phd/data/Sequencing_projects/Marine_June2017/dada2_merged/taxonomy_trimmed.csv';
write.table(tt,file=csvfile,append=FALSE,quote=FALSE,sep=",",eol = "\n", na = "NA", dec = ".",col.names=NA,row.names=TRUE);
#write.table(tt.plus,file="/Users/joshuagoldford/Dropbox/phd/research/coalescence/taxonomy_species_silva_61817.csv",append=FALSE,quote=FALSE,sep=",",eol = "\n", na = "NA", dec = ".",col.names=NA,row.names=TRUE);




