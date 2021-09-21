function xdot = nBody3(t,x)
%CASSINI + SATURN + TITAN + SUN
%Part d: Propagate position and velocity of Cassini forward in time using n-body equation
m1 =568.5e24; %mass of Saturn in kg;
m2 =2150; %mass of cassini orbiter in kg 
m3 = 1.345e23;% mass of Titan in kg
m4 = 1.989e30; %mass of the Sun in kg
G = 6.67408e-11; %gravitational constant
mu = (G*(m1+m2)) / (1000)^3; %mu Cassini/Saturn in km^3/s^2
mu3 = (G*m3) / (1000)^3;%mu for Titan
mu4 = (G*m4) / (1000)^3; %mu for Sun

%------------------------------------------
%r32: radius from titan to cassini
%r31: radius from titan to saturn
%r42: radius from sun to cassini
%r41: radius from sun  to saturn

[state31, lt31] = cspice_spkezr('SATURN', t , 'J2000', 'LT', 'TITAN'); %GOOD
r31 = state31(1:3); %radius from titan to saturn

r32 = r31 + x(1:3); %radius from titan to cassini

[state41, lt41] = cspice_spkezr('SATURN',t , 'J2000', 'LT', 'SUN'); %GOOD
r41 = state41(1:3); %radius from sun  to saturn

r42 = r41 + x(1:3); %radius from sun to cassini

%-----------------------------------------------

r = x(1:3);%radius from cassini to saturn
v = x(4:6); %acceleration of cassini relative to saturn
rho = norm(r);

rho32 = norm(r32);
rho31 = norm(r31);
r_hat32 = r32 / rho32;
r_hat31 = r31 / rho31;
termP1 = mu3 *( (r_hat32/(rho32)^2) - (r_hat31/(rho31)^2) );%first permutation term - titan

rho42 = norm(r42);
rho41= norm(r41);
r_hat42 = r42 / rho42;
r_hat41 = r41 / rho41;
termP2= mu4 *( (r_hat42/(rho42)^2) - (r_hat41/(rho41)^2) ); %second permutation term - sun


a = ((-mu/rho^3)*r) - (termP1 + termP2); %acceleration of cassini relative to saturn
xdot = [v;a];