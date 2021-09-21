function xdot = NBodyEC(t,x)
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
r_SM = [8500; 0; 0];% position of spacecraft relative to moon
%------------------------------------------
%r12: radius from earth to spacecraft
%r32: radius from moon to spacecraft
%r31: radius from moon to Earth

%-----------------------------------------------

r_SE= x(1:3); %position of spacecraft relative to Earth center
v_SE= x(4:6);%velocity of spacecaft relative to Earth center
r_MS = r_SE - r_EM; %position of moon relative to spacecraft
r_ME = -r_EM;%r_SE - r_SM; %position of moon relative to Earth

r32 = r_MS;
rho32 = norm(r32);

r31 = r_ME;
rho31 = norm(r31);

mu1 = mu_earth;
mu3 = mu_moon;
term1 = r32/rho32^3;
term2 = r31/rho31^3;

r12 = r_SE;
rho12 = norm(r12);


a = (-(mu1/rho12^3)*r12) -(mu3*(term1-term2)); %acceleration of spacecraft relative to barycenter
xdot = [v_SE;a];
