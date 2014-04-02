d = 1000
tau = 1
nu = 1
G = 10000
gam = 1/sqrt(rgamma(G, shape = nu/2, rate = nu * tau^2/2))

lp = function(nu){
  - lgamma(nu/2) + nu/2 * log(nu*tau^2/2) - nu* mean(log(gam) + tau^2/2 * 1/gam^2)
}

up = function(nu){
  exp(lp(nu))
}

library(pracma)
x = seq(0.01, 10*d, by = .1)
area = trapz(x = x, up(x))

p = function(x){
  up(x) / area
}

f = function(nu){
  dgamma(nu, shape = 3/2, rate = -1/2 - log(tau^2)/2 + mean(log(gam) + tau^2/2 * 1/gam^2))
}

curve(f, 0.001, 10, col = "blue", lwd = 2)
xsp = 1:(10*d)/10
ysp = p(xsp)
lines(x =xsp, y = ysp, lwd = 2)

legend("topright", legend = c("full conditional", "approximation"), fill = c("black", "blue"))

# x0 = 1000; f(x0)/p(x0)