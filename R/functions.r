heterosis_mcmc = function(data = "data.txt", group = "group.txt", out = "out",
                          chains = 2, iter = 10, burnin = iter/2,
                          phiPrior = 0, alphaPrior = 0, deltaPrior = 0,
                          seed = FALSE, joint = FALSE, verbose = FALSE, dic = FALSE,
                          hyper = FALSE, parms = FALSE, probs = FALSE, rates = FALSE,
                          time = FALSE, debug = FALSE, sigmaC0 = 10, d0 = 1000, aTau = 100,
                          aAlpha = 1, aDelta = 1, bTau = 100, bAlpha = 1, bDelta = 1,
                          gammaPhi = 2, gammaAlpha = 2, gammaDelta = 2, sigmaPhi0 = 2,
                          sigmaAlpha0 = 2, sigmaDelta0 = 2, sigmaC = NULL, d = NULL,
                          tau = NULL, thetaPhi = NULL, thetaAlpha = NULL, thetaDelta = NULL,
                          sigmaPhi = NULL, sigmaAlpha = NULL, sigmaDelta = NULL,
                          piAlpha = NULL, piDelta = NULL){
 
  argv = c("./heterosis-mcmc")

  if (is.character(data)){
    
    if(file.exists(data)){
      argv = c(argv, "--data", data)
    } else {
      stop(paste("Data file", data, "does not exist."))
    }

  } else if(is.matrix(data) || is.array(data) || is.data.frame(data)){
    write.table(data, "data.txt", row.names = F, col.names = F)
    argv = c(argv, "--data", "data.txt")
  } else {
    stop("bad data argument.")
  }

  if (is.character(group)){
    
    if(file.exists(group)){
      argv = c(argv, "--group", group)
    } else {
      stop(paste("Group file", group, "does not exist."))
    }

  } else if(is.numeric(group)){
    write(group, "group.txt")
    argv = c(argv, "--group", "group.txt")
  } else {
    stop("bad group argument.")
  }

  if(joint)   argv = c(argv, "--joint")
  if(verbose) argv = c(argv, "--verbose")
  if(dic)     argv = c(argv, "--dic")
  if(hyper)   argv = c(argv, "--hyper")
  if(parms)   argv = c(argv, "--parms")
  if(probs)   argv = c(argv, "--probs")
  if(rates)   argv = c(argv, "--rates")
  if(time)    argv = c(argv, "--time")
  if(debug)   argv = c(argv, "--debug")

  argv = c(argv, "--out", out, "--chains", paste(chains), "--iter", paste(iter), 
           "--burnin", paste(burnin), "--seed", paste(seed), 
           "--phi-prior", paste(phiPrior), "--alpha-prior", paste(alphaPrior), 
           "--delta-prior", paste(deltaPrior), "--sigma-c0", paste(sigmaC0),
           "--d0", paste(d0), "--a-tau", paste(aTau), "--a-alpha", paste(aAlpha), 
           "--a-delta", paste(aDelta), "--b-tau", paste(bTau), "--b-alpha", paste(bAlpha), 
           "--b-delta", paste(bDelta), "--gamma-phi", paste(gammaPhi), "--gamma-alpha",
           paste(gammaAlpha), "--gamma-delta", paste(gammaDelta), "--sigma-phi0",
           paste(sigmaPhi0), "--sigma-alpha0", paste(sigmaAlpha0), "--sigma-delta0",
           paste(sigmaDelta0))
 
  if(!is.null(sigmaC))     argv = c(argv, "--sigma-c", paste(sigmaC))
  if(!is.null(d))          argv = c(argv, "--d", paste(d))
  if(!is.null(thetaPhi))   argv = c(argv, "--theta-phi", paste(thetaPhi))
  if(!is.null(thetaAlpha)) argv = c(argv, "--theta-alpha", paste(thetaAlpha))
  if(!is.null(thetaDelta)) argv = c(argv, "--theta-delta", paste(thetaDelta))
  if(!is.null(sigmaPhi))   argv = c(argv, "--sigma-phi", paste(sigmaPhi))
  if(!is.null(sigmaAlpha)) argv = c(argv, "--sigma-alpha", paste(sigmaAlpha))
  if(!is.null(sigmaDelta)) argv = c(argv, "--sigma-delta", paste(sigmaDelta))
  if(!is.null(piAlpha))    argv = c(argv, "--pi-alpha", paste(piAlpha))
  if(!is.null(piDelta))    argv = c(argv, "--pi-delta", paste(piDelta))

  argc = length(argv)
  
  library.dynam("heterosis", "heterosis", lib.loc = NULL)

# SEGFAULTS:
# .Call("mcmc", as.integer(argc), as.character(argv), PACKAGE = "heterosis")

# THIS ONE AVOIDS ALL THE PROBLEMS OF CALLING C FROM R:
  argv[1] = paste(.libPaths(), "/heterosis/doc/heterosis-mcmc", sep="")
  cmd = paste(argv, collapse=" ")  
  system(cmd)

  if(is.matrix(data) || is.array(data) || is.data.frame(data))
    file.remove("data.txt")

  if(is.numeric(group))
    file.remove("group.txt")

  return(out)
}

test_heterosis_mcmc = function(){
  data(exampleData)
  data(exampleGroup)
  heterosis_mcmc(data = exampleData, group = exampleGroup, hyper = T, dic = T)
  gelmanFactors("out")
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