#In this home work, we are going to run RAxML, an external program.
#You will need to understand some of the arguments to RAxML and what they mean.
#Things for you to change or fix are set off with ______UNDERSCORES______
#Remember they might not have quotes in the correct answer (the seed, for example)

InferMorphologyTree_exercise <- function(in.place=FALSE, input.path=NULL, input.file = "binary.phy", output.path=NULL, output.name = "morpho1", random.seed=76185, model="ASC_BINGAMMA", other='--asc-corr=lewis') {
	if(!in.place) {
		if(is.null(input.path)) {
			fpath <- system.file("extdata", input.file, package="PhyloMethLikelihoodTrees")
			system(paste("cp ", fpath, " ", output.path,"/", input.file, sep=""))
		} else {
			system.call <- paste("cp ", input.path, "/", input.file, " ", output.path,"/", input.file, sep="")
			cat(system.call, file="~/Desktop/test_call.txt")
			system(system.call)	
		}
		cur.wd <- getwd()
		setwd(output.path)
	} 
	raxml.call <- paste("raxmlHPC -m ", model, " -p ", random.seed, " -s ", input.file, " ", other, " -n ", output.name, sep="")
	cat(raxml.call, file="~/Desktop/test_raxml.call")
	cat(getwd(), file="~/Desktop/test_getwd_testthat.call")
	status <- system(raxml.call)
	save(list=ls(), file='~/Desktop/test_dump.RData')
	parsimony.tree <- read.tree(paste("RAxML_parsimonyTree.",output.name, sep=""))
	ml.tree <- read.tree(paste("RAxML_bestTree.",output.name, sep=""))
	if(!in.place) {
		setwd(cur.wd)
	} else {
		if(!is.null(output.path)) {
			try(system(paste("mv RAxML* ", output.path, sep="")))
		} else {
			system("rm RAxML*")
		}
	}
	return(list(parsimony.tree=parsimony.tree, ml.tree=ml.tree))
}

InferDNATreeWithBootstrappingAndPartitions_exercise <- function (in.place=FALSE, input.path=NULL, input.file = "dna.phy", input.partition = "dna12_3.partition.txt", output.path=NULL, output.name = "dna1", random.seed=12345, boot.seed=12345, model="GTRGAMMA", boot=100) {
	if(!in.place) {
		if(is.null(input.path)) {
			fpath <- system.file("extdata", input.file, package="PhyloMethLikelihoodTrees")
			system(paste("cp ", fpath, " ", output.path,"/", input.file, sep=""))
			fpath <- system.file("extdata", input.partition, package="PhyloMethLikelihoodTrees")
			system(paste("cp ", fpath, " ", output.path,"/",input.partition, sep=""))
		} else {
			system(paste("cp ", input.path, "/", input.file, " ", output.path,"/", input.file, sep=""))
			system(paste("cp ", input.path, "/", input.partition, " ", output.path,"/",input.partition, sep=""))
		}
		cur.wd <- getwd()
		setwd(output.path)
	} 

	raxml.call <- paste("raxmlHPC -k -f a", " -m ", model, " -p ", random.seed, " -x ", boot.seed, " -q ", input.partition, " -# ", boot, " -s ", input.file, " -n ", output.name, sep="")
	status <- system(raxml.call)
	
	ml.tree <- read.tree(paste("RAxML_bestTree.",output.name, sep=""))
	ml.with.bs.tree <- read.tree(paste("RAxML_bipartitions.",output.name, sep=""))
	bs.trees <- read.tree(paste("RAxML_bootstrap.",output.name, sep=""))
	if(!in.place) {
		setwd(cur.wd)
	} else {
		if(!is.null(output.path)) {
			try(system(paste("mv RAxML* ", output.path, sep="")))
		} else {
			system("rm RAxML*")
		}
	}
	return(list( ml.tree=ml.tree, ml.with.bs.tree=ml.with.bs.tree, bs.trees))
}