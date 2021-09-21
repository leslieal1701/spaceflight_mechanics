function xdot = TwoBodyMS(t,x)
%Spacecraft + Moon
m1=  73.48e21;% mass of moon in kg
G = 6.67408e-11; %gravitational constant
mu = (G*m1) / (1000)^3; %km^3/s^2

r = x(1:3); %position vector of cassini relative to saturn
v = x(4:6); %velocity vector of cassini relative to saturn

rho = norm(r); 

a = (-mu/rho^3)*r; %acceleration of cassini relative to saturn
xdot = [v;a];