d = 100
tau = 1
nu = 1
G = 10000
gam = 1/sqrt(rgamma(G, shape = nu/2, rate = nu * tau^2/2))

lp = function(nu){
 # nu[nu <= 0] = 0.0001
  - lgamma(nu/2) + nu/2 * log(nu*tau^2/2) - nu* mean(log(gam) + tau^2/2 * 1/gam^2)
}

p = function(nu){
  exp(lp(nu))
}

wd = 0.5
bd = 10

f = Vectorize(function(x){
  i = x / wd
  if(i <= 1)
    p(wd)
  else
    max(p(floor(i) * wd), p(ceiling(i) * wd))
}, "x")

ys = f(1:(bd / wd)*wd + wd/2)
cy = cumsum(ys)
sampleF = function(){
  u = runif(1, 0, max(cy))

  for(i in length(ys):1)
    if(cy[i] < u)
      break

  m = min(ys[i], ys[i + 1])  
  M = max(ys[i], ys[i + 1])

  (runif(1) + i - 1) * wd 
}

sampleFs = function(n){
  ret = c()
  for(i in 1:n)
    ret[i] = sampleF()
  ret
}

ss = sampleFs(10000)

hist(ss,  freq = F, col = "red")
curve(f, 0.001, 10, col = "blue", add = T, lwd = 2)
xsp = 1:10000/1000
ysp = p(xsp)
lines(x =xsp, y = ysp, lwd = 2)

legend("topright", fill = c("black", "blue"), legend = c("Const * exp(h(nu))", "Approximation"))

