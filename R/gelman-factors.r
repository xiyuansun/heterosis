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

gelmanFactors = function(mainDir){
  cwd = getwd()

  if(!file.exists(mainDir)){
    print(paste("ERROR: could not open directory, ", mainDir, ".", sep =""))
    return();
  }

  setwd(mainDir)

  if(!file.exists("diagnostics"))
    dir.create("diagnostics")

  kinds = c("hyper-", "parms-")

  found = 0;
  for(kind in kinds)
    if(file.exists(paste(kind, "chain1.txt", sep = "")))
      found = found + 1

  if(!found){
    print("ERROR: no parameters found.")
    return()
  }
 
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

options <- commandArgs(trailingOnly = TRUE)
gelmanFactors(options[1])