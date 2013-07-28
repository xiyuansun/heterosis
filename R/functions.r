heterosis_mcmc = function(){
  dyn.load(paste(.libPaths(), "/heterosis/libs/heterosis.so"))
}

gpu_heterosis_mcmc = function(){
  dyn.load(paste(.libPaths(), "/heterosis/libs/heterosis.so"))
}

checkVersion = function(){
  v = R.Version()
  major = as.numeric(v$major)
  minor = as.numeric(v$minor)

  if(major <= 2 && minor < 15)
    stop(paste("R version 2.15.0 or later is required. This is ", 
               v$version.string, ".", sep=""))
}

oneParm = function(kind, parmName, parmNum, outFile, nfiles){
  library(coda, quietly = T)

  l = mcmc.list()
  for(i in 1:nfiles){
    cmd = paste("cut -d ' ' -f", parmNum, " ", kind, "chain", i, ".txt", sep="")
    con = pipe(cmd)
    l[[i]] = mcmc(as.numeric(scan(con, what = character(), quiet = T)[-1]))
    close(con)
  }

  gelman = paste(round(unlist(gelman.diag(l, autoburnin = F, multivariate = F)),3), collapse = " ")
  write(paste(parmName, gelman, collapse = " "), outFile, append = T, sep = " ")
}

readHeader = function(kind){
  con = file(paste(kind, "chain1.txt", sep = ""), "r") 
  h = readLines(con, 1)
  close(con)
  strsplit(h, split = " ")[[1]]
}

oneKind = function(kind, outFile, mainDir){

  h = readHeader(kind)

  nfiles = 0
  files = list.files()

  for(file in files)
     nfiles = nfiles + length(grep(kind, file))

  if(!nfiles)
    return()

  for(i in 1:length(h))
    oneParm(kind, h[i], i, outFile, nfiles)
}

gelmanFactors = function(mainDir = "out", parms = F){

  checkVersion()
  cwd = getwd()

  if(!file.exists(mainDir))
    stop(paste("Could not open directory, ", mainDir, ".", sep =""))

  setwd(mainDir)

  kinds = c("hyper-")

  if(parms)
    kinds = c(kinds, "parms-")

  found = 0;
  for(kind in kinds)
    if(file.exists(paste(kind, "chain1.txt", sep = "")))
      found = found + 1

  if(!found)
    stop("No parameters found.")
 
  outFile = "gelman-factors.txt"    
  write("parameter gelman-point-est 95%-upper-bd", outFile, sep = " ")

  for(kind in kinds)
    if(file.exists(paste(kind, "chain1.txt", sep = ""))){
      if(kind == "hyper-"){
        print(paste("Computing Gelman factors on hyperparameters."))
      } else {
        print(paste("Computing Gelman factors on non-hyper parameters."))
      }

      oneKind(kind, outFile, mainDir)
    }

  print(paste("Please find Gelman factors in gelman-factors.txt within ", 
              mainDir, ".", sep = ""))
  setwd(cwd)
}

#options <- commandArgs(trailingOnly = TRUE)
#gelmanFactors(options[1])