function xdot = TwoBody(t,x)
%Cassini + Saturn
m2 =2150; %mass of cassini orbiter in kg 
m1 = 568.5e24; %mass of Saturn in kg
G = 6.67408e-11; %gravitational constant
mu = (G*(m1+m2)) / (1000)^3; %km^3/s^2

r = x(1:3); %position vector of cassini relative to saturn
v = x(4:6); %velocity vector of cassini relative to saturn

rho = norm(r); 

a = (-mu/rho^3)*r; %acceleration of cassini relative to saturn
xdot = [v;a];