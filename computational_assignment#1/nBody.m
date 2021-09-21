function xdot = nBody(t,x)
% SATURN + CASSINI  + TITAN

%%Part d: Propagate position and velocity of Cassini forward in time using n-body equation
m1 =568.5e24; %mass of Saturn in kg
m2 =2150; %mass of cassini orbiter in kg 
m3 = 1.345e23;% mass of Titan in kg
G = 6.67408e-11; %gravitational constant
mu = (G*(m1+m2)) / (1000)^3; %mu Cassini/Saturn in km^3/s^2
mu3 = (G*m3) / (1000)^3;%mu for Titan
%------------------------------------------
%r32: radius from titan to cassini
%r31: radius from titan to saturn
%r12: radius from saturn to cassini -- x(1:3)
[state31, lt31] = cspice_spkezr('SATURN', t , 'J2000', 'LT', 'TITAN'); 
r31 = state31(1:3); %radius from titan to saturn

r32 = r31 + x(1:3); %radius from titan to cassini

%-----------------------------------------------

r = x(1:3); %radius from saturn to cassini
v = x(4:6);%velocity  of cassini 
rho = norm(r);

rho32 = norm(r32);
rho31 = norm(r31);
r_hat32 = r32 / rho32;
r_hat31 = r31 / rho31;

termP1 = mu3 *( (r_hat32/(rho32)^2) - (r_hat31/(rho31)^2) ); %Permutation term


a = ( (-mu/rho^3)*r ) - termP1; %acceleration of cassini relative to saturn
xdot = [v;a];

