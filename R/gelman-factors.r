oneParm = function(dir, parmName, parmNum, outCon, nfiles){
  library(coda, quietly = T)

  l = mcmc.list()
  for(i in 1:nfiles){
    cmd = paste("cut -d ' ' -f", parmNum, " ", dir, "chain", i - 1, ".txt", sep="")
    con = pipe(cmd)
    l[[i]] = mcmc(as.numeric(scan(con, what = character(), quiet = T)[-1]))
    close(con)
  }

  gelman = paste(round(unlist(gelman.diag(l, autoburnin = F, multivariate = F)),3), collapse = " ")
  write(paste(parmName, gelman, collapse = " "), outCon, append = T, sep = " ")
}

readHeader = function(dir){
  con = file(paste(dir, "chain0.txt", sep = ""), "r") 
  h = readLines(con, 1)
  close(con)
  strsplit(h, split = " ")[[1]]
}

getOutCon = function(outFile){
  if(!file.exists(outFile))
    file.create(outFile)

  outCon = file(outFile, "w")
}

oneDir = function(dir, outCon){
  if(substr(dir, nchar(dir), nchar(dir)) != "/")
    dir = paste(dir, "/", sep = "")

  h = readHeader(dir)

  con = pipe("ls | wc -l")
  nfiles = scan(con, quiet = T)
  close(con)

  if(!nfiles)
    return

  for(i in 1:length(h))
    oneParm(dir, h[i], i, outCon, nfiles)
}

gelmanFactors = function(){
  outFile = "../out/diagnostics/gelman-factors.txt"  

  outCon = getOutCon(outFile)
  write("parameter gelman-point-est 95%-upper-bd", outCon, append = T, sep = " ")

  dirs = paste("../out/", c("hyper/", "parms/"), sep = "")

  for(dir in dirs)
    if(file.exists(dir))
      oneDir(dir, outCon)

  close(outCon)
}
