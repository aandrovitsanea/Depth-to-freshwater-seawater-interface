# Depth to freshwater-seawater interface for different inland flux with 0.5, 0.2 and 0.1 $m^2/s$ and a recharge rate of 100 mm/year

H  = 100 # m
L  = 1000 #m
a  = 0.1*L
va = 1E-6
Dc  = va*a
e  = 40 #
es = e*(1-a/H)^(1/6)
K  = 10 # m/d
Qf = -0.5 # sm per d
Qf1= -0.25
Qf2= -0.1
N  = 0.001 # m per d
x  = seq(0,1000,1)
hs = sqrt(((Qf *x-N*L*x+N*(x^2/2))*(2*e^2))/(-(1+e)*K))
y  = hs
plot(y~x , ylim = rev(range(y)) ,
     lwd=4 , type="l" , bty="n" , ylab="depth below sealevel" , col=rgb(0.2,0.4,0.6,0.8) )
curve(sqrt(((Qf *x-N*L*x+N*(x^2/2))*(2*e^2))/(-(1+e)*K)), add = TRUE, col="blue")
curve(sqrt(((Qf1*x-N*L*x+N*(x^2/2))*(2*e^2))/(-(1+e)*K)), add = TRUE, col="red")
curve(sqrt(((Qf2*x-N*L*x+N*(x^2/2))*(2*e^2))/(-(1+e)*K)), add = TRUE, col="green")


# Depth to freshwater-seawater interface for a confined aquifer with a flux of $0.5 m^2/d$, aquifer thickness of 100 m and a hydraulic conductivity of 1 m/d.

H  = 100 # m aquifer thickness at coast
e  = 40 #
K  = 1 # m/d hydraulic conductivity
L  = 100 # m
Qf = 0.5 # sm per d
h  = seq(0,100,1)
x =  K/(e*Qf)*(h^2/2)
y  = h
plot(y~x , ylim = rev(range(y)) ,
     lwd=4 , type="l" , bty="n" , ylab="depth below ground water level (m)" , col=rgb(0.2,0.4,0.6,0.8) )
     

# Depth to transition zone for an unconfined aquifer and a confined aquifer

# pumping_unconfined
rs = 1.024
rh = 1.000
e  = (rs-rh)/rh #
H  = 100 # aquifer thickness at coast
K  = 1 # m/s
K1 = 5
K2 = 0.5
Qw = 1 # qm/s
Qw2= 25
Qx0= 1 # qm/s
p  = pi
x  = seq(1,500,1)
mu = Qw/(Qx0*x)
lambda = 2*(1-mu/p)^(1/2)+mu/p*log((1-(1-mu/p)^(1/2))/(1+(1-mu/p)^(1/2)))
y = sqrt((lambda*Qx0*x)/(K*e))


plot(y~x,ylim = rev(range(y)),lwd=4 , type="l" , bty="n" , ylab="depth below ground water level (m)")
curve(sqrt(((2*(1-(Qw/(Qx0*x))/pi)^(1/2)+(Qw/(Qx0*x))/pi*log((1-(1-(Qw/(Qx0*x))/pi)^(1/2))/(1+(1-(Qw/(Qx0*x))/pi)^(1/2))))*Qx0*x)/(K1*e)), add = TRUE, col="blue")
curve(sqrt(((2*(1-(Qw/(Qx0*x))/pi)^(1/2)+(Qw/(Qx0*x))/pi*log((1-(1-(Qw/(Qx0*x))/pi)^(1/2))/(1+(1-(Qw/(Qx0*x))/pi)^(1/2))))*Qx0*x)/(K2*e)), add = TRUE, col="green")
curve(sqrt(((2*(1-(Qw2/(Qx0*x))/pi)^(1/2)+(Qw2/(Qx0*x))/pi*log((1-(1-(Qw2/(Qx0*x))/pi)^(1/2))/(1+(1-(Qw2/(Qx0*x))/pi)^(1/2))))*Qx0*x)/(K*e)), add = TRUE, col="red")

# pumping_confined
rs = 1.024
rh = 1.000 
e  = (rs/rh^2)*(rs-rh) #
a  = 0.1*10
H  = 100 # aquifer thickness at coast
K  = 1 # m/s
K1  = 2 
Qw = 1 # qm/s
Qx0= 2 # qm/s
p  = pi
x  = seq(1,500,1)
mu = Qw/(Qx0*x)
lambda = 2*(1-mu/p)^(1/2)+mu/p*log((1-(1-mu/p)^(1/2))/(1+(1-mu/p)^(1/2)))
es = e*(1-(a/H)^(1/6))
y = sqrt((lambda*Qx0*x)/(K*es))
plot(y~x,ylim = rev(range(y)),lwd=4 , type="l" , bty="n" , ylab="depth below ground water level (m)")
curve(sqrt(((2*(1-(Qw/(Qx0*x))/pi)^(1/2)+(Qw/(Qx0*x))/pi*log((1-(1-(Qw/(Qx0*x))/pi)^(1/2))/(1+(1-(Qw/(Qx0*x))/pi)^(1/2))))*Qx0*x)/(K*e)), add = TRUE, col="blue")
