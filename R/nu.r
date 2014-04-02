d = 100
tau = 1
nu = 1
G = 10000
gam = 1/sqrt(rgamma(G, shape = nu/2, rate = nu * tau^2/2))

lp = function(nu){
  - lgamma(nu/2) + nu/2 * log(nu*tau^2/2) - nu* mean(log(gam) + tau^2/2 * 1/gam^2)
}

p = function(nu){
  exp(lp(nu))
}

# insert new gamma approximation here

f = function(nu){
  dgamma(nu, shape = 3/2, rate = -1/2 - log(tau^2)/2 + mean(log(gam) + tau^2/2 * 1/gam^2))
}

curve(f, 0.001, 10, col = "blue", lwd = 2)
xsp = 1:10000/1000
ysp = p(xsp)
lines(x =xsp, y = ysp, lwd = 2)

legend("topright", legend = c("full conditional", "approximation"), fill = c("black", "blue"))