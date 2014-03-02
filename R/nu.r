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

x = 1:10000/1000
m = sum(x * p(x))/sum(p(x))
v = sum((x-m)^2 * p(x))/sum(p(x))

shape = m^2/v
rate = m/v

f = function(x){dgamma(x, shape = shape, rate = rate) /2}
curve(f, 0, 10, col = "blue")
xs = 1:10000/1000
ys = p(xs)/max(p(xs)) * .15
lines(x =xs, y = ys)
legend("topright", fill = c("black", "blue"), legend = c("Const * exp(h(nu))", "Gamma approximation"))