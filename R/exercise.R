#installed.packages("phangorn")
#help(package = "phangorn")
require(ape)
require(phangorn)
require (ips)
##################################### Morphology :/ (this face means "I Don't if I did this exercise right")


InferMorphologyTree_exercise <- function(in.place=FALSE, input.path=NULL, input.file = "binary.phy", output.path=NULL, 
                                         output.name = "morpho1", random.seed=12345, #found this on Exelixis Lab site (or I can use "startingTree.txt" BUT we don't have one.)
                                         model="ASC_BINGAMMA", # Looked for this at raxml -h after I saw the "solution".
                                         other='--asc-corr=lewis') {
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
####### DNA :)

InferMorphologyTree_solution <- function(input.path = "inst/extdata/", input.file = "binary.phy", output.path="~/Desktop/", output.name = "morpho1", random.seed=12345, model="ASC_BINGAMMA", other='--asc-corr=lewis') {
  system(paste("cp ", input.path, input.file, " ", output.path,input.file, sep=""))
  InferMorphologyTree_solution <- function(input.path=NULL, input.file = "binary.phy", output.path="~/Desktop", output.name = "morpho1", random.seed=12345, model="ASC_BINGAMMA", other='--asc-corr=lewis') 
    if(is.null(input.path)) {
      fpath <- system.file("extdata", input.file, package="PhyloMethLikelihoodTrees")
      system(paste("cp ", fpath, " ", output.path,"/", input.file, sep=""))
    } else {
      system.call <- paste("cp ", input.path, "/", input.file, " ", output.path,"/", input.file, sep="")
      cat(system.call, file="test_call.txt")
      system(system.call)
    }
  
  cur.wd <- getwd()
  setwd(output.path)
  raxml.call <- paste("raxmlHPC -m ", model, " -p ", random.seed, " -s ", input.file, " ", other, " -n ", output.name, sep="")
  status <- system(raxml.call)
  system(paste("cp -r ", tempdir(), " ~/Desktop/test_tempdir", sep=""))
  save(list=ls(), file='~/Desktop/test_dump.RData')
  parsimony.tree <- read.tree(paste("RAxML_parsimonyTree.",output.name, sep=""))
  ml.tree <- read.tree(paste("RAxML_bestTree.",output.name, sep=""))
  setwd(cur.wd)
  return(list(parsimony.tree=parsimony.tree, ml.tree=ml.tree))
}





InferDNATreeWithBootstrappingAndPartitions_exercise <- function (in.place=FALSE, input.path=NULL, input.file = "dna.phy",
                                                                 input.partition = "dna12_3.partition.txt", output.path=NULL, 
                                                                 output.name = "dna1", random.seed=12345, boot.seed=12345,
                                                                 model="GTRGAMMA", boot=100) { #found at Exelixis Lab site
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
  
  raxml.call <- paste("raxmlHPC -f a  -k -m GTRGAMMA"," -m ", model, " -p ", random.seed, " -x ", boot.seed, " -q ",
                      input.partition, " -# ", boot, " -s ", input.file, " -n ", output.name, sep="")
  ### in terminal raxml -h (help manual) I found: 
  # "-f a": rapid Bootstrap analysis and search for best-scoring ML tree in one program run 
  # "-k": Specifies that bootstrapped trees should be printed with branch lengths.
  # "-m GTRGAMMA":  GTR + Optimization of substitution rates + Optimization of site-specific
  #evolutionary rates which are categorized into numberOfCategories distinct 
  #rate categories for greater computational efficiency.  Final tree might be evaluated
  #under GTRGAMMA, depending on the tree search option.
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

