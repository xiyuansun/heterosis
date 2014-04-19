d = 1000
u = 8
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
x = seq(0.001, d, by = .01)
area = trapz(x = x, up(x))

p = function(x){
  up(x) / area
}

#md = x[which.max(lp(x))]
rate = -1/2 - log(tau^2)/2 + mean(log(gam) + tau^2/2 * 1/gam^2) 
shape = 3/2 #md * rate + 1

rate = rate * 0.99

mn = shape/rate
v = shape/rate^2

f = function(nu){
  dgamma(nu, shape = shape, rate = rate)
}

curve(p, 0.001, u, col = "black", lwd = 2)
xsp = seq(0.001, u, by = .01)
ysp = f(xsp)
lines(x =xsp, y = ysp, lwd = 2, col = "blue")

legend("topright", legend = c("full conditional", "approximation"), fill = c("black", "blue"))

x0 = 1000; f(x0)/p(x0)