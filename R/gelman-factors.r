oneParm = function(dir, parmName, parmNum, outFile, nfiles){
  library(coda, quietly = T)

  l = mcmc.list()
  for(i in 1:nfiles){
    cmd = paste("cut -d ' ' -f", parmNum, " ", dir, "chain", i, ".txt", sep="")
    con = pipe(cmd)
    l[[i]] = mcmc(as.numeric(scan(con, what = character(), quiet = T)[-1]))
    close(con)
  }

  gelman = paste(round(unlist(gelman.diag(l, autoburnin = F, multivariate = F)),3), collapse = " ")
  write(paste(parmName, gelman, collapse = " "), outFile, append = T, sep = " ")
}

readHeader = function(dir){
  con = file(paste(dir, "chain1.txt", sep = ""), "r") 
  h = readLines(con, 1)
  close(con)
  strsplit(h, split = " ")[[1]]
}

oneDir = function(dir, outFile){
  if(substr(dir, nchar(dir), nchar(dir)) != "/")
    dir = paste(dir, "/", sep = "")

  h = readHeader(dir)
  nfiles = length(list.files(dir)) 

  if(!nfiles)
    return

  for(i in 1:length(h))
    oneParm(dir, h[i], i, outFile, nfiles)
}

gelmanFactors = function(outFile = "../out/diagnostics/gelman-factors.txt"){


  if(!file.exists("../out/")){
    print("ERROR: no mcmc output to work with.")
    return(NULL);
  }

  if(!file.exists("../out/diagnostics"))
    dir.create("../out/diagnostics")
    
  write("parameter gelman-point-est 95%-upper-bd", outFile, append = T, sep = " ")

  dirs = paste("../out/", c("hyper/", "parms/"), sep = "")

  for(dir in dirs)
    if(file.exists(dir))
      oneDir(dir, outFile)
}

gelmanFactors()