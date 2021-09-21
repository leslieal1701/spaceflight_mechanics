function xdot = NBodyEMB(t,x)
% SPACECRAFT + EARTH + MOON

%%Part 4: Propagate position and velocity of Cassini forward in time using n-body equation
mu_earth = 3.986e5; %km^3/s^2
mu_moon =  0.0123*mu_earth;
mu = mu_earth + mu_moon;
r_em = 384399;%distance from earth to moon
omega = sqrt(mu/(r_em)^3); %angular velocity
r_EM = r_em*[cos(omega*t);sin(omega*t);0];
r_EB = ( (-mu_moon/mu)*r_em ).*[cos(omega*t);sin(omega*t);0];% position of earth relative to barycenter
r_MB = ( (mu_earth/mu)*r_em).*[cos(omega*t);sin(omega*t);0]; %position of moon relative to barycenter
%------------------------------------------
%r12: radius from earth to spacecraft
%r32: radius from moon to spacecraft

%-----------------------------------------------

r_SB= x(1:3); %position of spacecraft relative to barycenter
v_BS= x(4:6);%velocity of spacecaft relative to barycenter
r12 = r_SB- r_EB; %position of spacecraft relative to Earth
r32 = r_SB - r_MB ;
rho12 = norm(r12);
rho32 = norm(r32);

mu1 = mu_earth;
mu3 = mu_moon;
term1 = (mu1/rho12^3)*r12;
term3 = (mu3/rho32^3)*r32;


a = -(term1+term3); %acceleration of spacecraft relative to barycenter
xdot = [v_BS;a];
