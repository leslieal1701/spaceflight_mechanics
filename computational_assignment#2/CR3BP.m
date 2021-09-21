function xdot = CR3BP(t,x)

mu_earth = 3.986e5; %km^3/s^2
mu_moon =  0.0123*mu_earth;
mu = mu_earth + mu_moon;
r_em = 384399;%distance from earth to moon
omega = sqrt(mu/(r_em)^3); %angular velocity
%-----------------------------------------------

r_SB3= x(1:3); %position of spacecraft relative to barycenter
v_SB3= x(4:6);%velocity of spacecaft relative to barycenter
r = x;

rEarthCR3BP = [(-mu_moon/mu)*r_em;0;0];
rMoonCR3BP = [(mu_earth/mu)*r_em;0;0];

r1 = r_SB3 - rEarthCR3BP; % spacecraft relative to earth in CR3BP frame
r2 = r_SB3 - rMoonCR3BP;
mu1 = mu_earth;
mu2 = mu_moon;
p1 = mu1/mu;
p2 = mu2/mu;
r12 = r_em;

a1 = (-mu1/norm(r1)^3).*[r(1)+(p2*r_em);r(2);r(3)];
a2 = (-mu2/norm(r2)^3).*[r(1)-(p1*r12);r(2);r(3)];
a3 = (-2).*[r(5)*(-omega);r(4)*omega;0];
a4 = (-1)*[r(1)*((-1)*(omega^2));r(2)*((-1)*(omega^2));0];
neta = a1+a2+a3+a4;
xdot = [v_SB3;neta];
