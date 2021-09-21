% MANE 4100 Computational Assignment #2
% By: LESLIE ALEMAN

%% Clear Workspace
clear all
close all
clc

%% Part 1
mu_earth = 3.986e5; %km^3/s^2
mu_moon =  0.0123*mu_earth;
r_em = 384399; %distance from earth to moon km

%Eqs to compute the position and velocity of the Earth-Moon system relative to the
%barycenter expressed in the INERTIAL FRAME

%Patched conics method for lunar trajectory
%Earth plays the role of the Sun
mass_moon =  73.48e21;% mass of moon in kg
mass_earth = 5.974e24; %mass of earth in kg
R_moon = 384.4e3; %semi-major axis of moon orbit in km
r_SOI_moon = R_moon*(mass_moon/mass_earth)^(2/5); %Size of moon's sphere of influence

%% Part 2
r_SM = [8500;0;0]; %position of spacecraft relative to moon in km
v_SM = [166.9;969.1;0]*10^-3; %velocity of spacecraft relative to moon in km/s
v_SM_norm = norm(v_SM);
%-------
E = (v_SM_norm)^2 / 2 - mu_moon/r_SM(1); %Energy of moon-centered orbit
a = -mu_moon / (2*E); %semi-major axis of moon-centered orbit
T_circular = ((2*pi) / sqrt(mu_moon))*(a)^(3/2);%period of moon-centered orbit

time = linspace(0,T_circular ,1000);
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,x] = ode45(@TwoBodyMS, time,[r_SM;v_SM],options);
r = x(:,1:3);%spacecraft relative to moon
figure(1)
plot(r(:,1),r(:,2)), axis equal, grid on, hold on
plot(0,0,'.k','MarkerSize',20)
xlabel('MCI Position in i Direction (km)')
ylabel('MCI Position in j Direction (km)')
title('Two-Body Orbit: Spacecraft + Moon')
legend('Satellite Orbit', 'Moon')
hold off

%% Part 3

%Compute the position and inertial velocity of the spacecraft at t0 relative to the Earth-Moon barycenter

mu = mu_earth + mu_moon;
omega = sqrt(mu/(r_em)^3); %angular velocity
t0=0;
r_EB = ( (-mu_moon/mu)*r_em )*[cos(omega*t0);sin(omega*t0);0];% position of earth relative to barycenter
r_MB = ( (mu_earth/mu)*r_em )*[cos(omega*t0);sin(omega*t0);0]; %position of moon relative to barycenter

r_SB = r_SM + r_MB; %position of spacecraft relative to barycenter
v_SB = v_SM + cross([0;0;omega],r_MB);%velocity of spacecraft relative to barycenter

r_SE = r_SM+ [r_em;0;0]; %position of spacecraft relative to Earth center
v_SE = v_SM + cross([0;0;omega],[r_em;0;0]); %velocity of spacecraft relative to Earth center


%% Part 4
%Propagate postion and velocity of spacecraft using the N-Body in
%inertially aligned frame centered at the earth-moon barycenter

tf =2*T_circular;
time4 = linspace(t0,tf,5000)';

options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[tB4,xB4] = ode45(@NBodyEMB, time4,[r_SB;v_SB],options); %Spacecraft + Earth + Moon in EMB frame
rB4 = xB4(:,1:3);
figure(2)
plot(rB4(:,1),rB4(:,2),'r'), axis equal, grid on, hold on
title('3 Body Orbit about Earth-Moon Barycenter Inertial Frame')
xlabel('Position in i direction relative to Barycenter(km)')
ylabel('Position in j direction relative to Barycenter(km)')
plot(0,0,'.k','MarkerSize',20)
hold off
legend('Satellite Orbit','Barycenter')

%% Part 5
%Integrate the equations of motion in an inertially aligned frame centered at Earth

options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[tB5,xB5] = ode45(@NBodyEC, time4,[r_SE;v_SE],options); %Spacecraft + Earth + Moon in EC frame
rB5 = xB5(:,1:3); %position of spacecraft relative to Earth center
figure(3)
plot(rB5(:,1),rB5(:,2),'r'), axis equal, grid on, hold on
title('3 Body Orbit about Earth Centered Inertial Frame')
xlabel('Position in i direction relative to Earth Center(km)')
ylabel('Position in j direction relative to Earth Center(km)')
plot(0,0,'.k','MarkerSize',20)
hold off 
legend('Satellite Orbit','Earth Center')

%% Part 6
% Transform the results of #5 to the E-M barycenter to perform the comparison.
r_EBt = ( (-mu_moon/mu)*r_em ).*[cos(omega.*time4), sin(omega.*time4),zeros(1, length(time4))']; %postion of earth relative to barycenter at each time

%origin shift from Earth centered inertial frame to E-M barycenter frame 
EC2EMB_shift = rB5(:,1:3) + r_EBt;
error = rB4(:,1:3)-EC2EMB_shift;
figure(4)
plot(time4,error(:,1),'.r','MarkerSize',10), grid on
hold on
plot(time4,error(:,2),'.m','MarkerSize',10), grid on
hold on
plot(time4,error(:,3),'.k','MarkerSize',10), grid on
hold off
xlabel('Time(s)')
ylabel('Position Error(km)')
legend('Error Along i Direction','Error Along j Direction','Error Along k Direction','Location','northwest')
title('Difference Between N-Body Solutions')
%{
figure(5)
plot(rB4(:,1),rB4(:,2),'r'), axis equal, grid on, hold on
plot(rB5(:,1),rB5(:,2),'g'), axis equal, grid on, hold on
plot(EC2EMB_shift(:,1),EC2EMB_shift(:,2),'k'), axis equal, grid on
hold off
legend('Part 4','Part 5','Shifted')
%}
%% Part 7
% Transform the three-body result from #5  to be centered about the Moon. In a
%Moon-centered inertial (MCI) frame
r_EMt = r_em.*[cos(omega.*time4), sin(omega.*time4),zeros(1, length(time4))']; %distance from earth to moon at each time
r_SMvec = [8500.*ones(1, length(time4))',zeros(1, length(time4))',zeros(1, length(time4))'];
%origin shift from Earth centered inertial frame to Moon-centered inertial frame

EC2MCI_shift = rB5(:,1:3) - r_EMt;
%-----------------------

x = 0;
y = 0;
th = 0:pi/50:2*pi;
xunit = r_SOI_moon * cos(th) + x;
yunit = r_SOI_moon * sin(th) + y;
figure(6)
plot(xunit, yunit),axis equal, grid on, hold on % plot sphere of influence
plot(0,0,'.k','MarkerSize',20),axis equal, grid on, hold on %moon
plot(r(:,1),r(:,2)),axis equal, grid on, hold on %2 body
plot(EC2MCI_shift(:,1),EC2MCI_shift(:,2),'g'), axis equal, grid on %3 body
hold off
l1 = 'Moon SOI';
l2= 'Moon';
l3 = '2-Body Trajectory';
l4 = '3-Body Trajectory';
legend(l1,l2,l3,l4);
xlabel('Position in i Direction(km)')
ylabel('Position in j Direction(km)')
title('Two-Body and Three-Body Motion About Moon-Centered Inertial Frame')

%% Part 8
%Compute the spacecraft position and velocity at t0 in the rotating CR3BP frame

r0_CR3BP = r_SB;%position of spacecraft relative to barycenter
v0_CR3BP= v_SB - cross([0;0;omega],r_SB); %velocity of spacecraft relative to barycenter

%% Part 9
% Integrate the EOM in the rotating CR3BP frame. 

options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[tB9,xB9] = ode45(@CR3BP, time4,[r0_CR3BP;v0_CR3BP],options); 

rB9 = xB9(:,1:3); 
figure(7)
plot(rB9(:,1),rB9(:,2),'r'), axis equal, grid on,hold on
plot(r_MB(1),0,'.k','MarkerSize',50), axis equal, grid on
hold off
title('Circular Restricted 3-Body Motion')
xlabel('Position in i direction(km)')
ylabel('Position in j direction(km)')
legend('Satellite Orbit','Moon')

%% Part 10
% Transform the results in an inertially aligned frame from part#4
%into the CR3BP frame 

%First, transform from CR3BP to inertial frame with transformation matrix


for i=1:size(rB4)
    phi = omega*time4(i);
   
    rotMatrix= [cos(phi),sin(phi),0; 
                -sin(phi),cos(phi),0;
                0,0,1];
            
    trans = rotMatrix*rB4(i,:)';
    P4_pos_CR3BP(i,:) = trans';
end



error10 = rB9(:,1:3)-P4_pos_CR3BP ;
figure(8)
plot(time4,error10(:,1),'.r','MarkerSize',10), grid on
hold on
plot(time4,error10(:,2),'.m','MarkerSize',10), grid on
hold on
plot(time4,error10(:,3),'.k','MarkerSize',10), grid on
hold off
xlabel('Time(s)')
ylabel('Position Error(km)')
legend('Error Along i Direction','Error Along j Direction','Error Along k Direction','Location','northwest')
title('Difference Between Circular Restricted 3 Body Problem Results')

%% Part 11

x9 = r0_CR3BP(1);
y9 = r0_CR3BP(2);
%r12 = r_em
p1 = mu_earth/mu;
p2 = mu_moon/mu; 
r1 = (sqrt( (x9+p2*r_em)^2+y9^2 ));
r2 = (sqrt( (x9-p1*r_em)^2+y9^2 ));

C = .5*norm(v0_CR3BP)^2 - .5*omega^2*(r0_CR3BP(1)^2 +r0_CR3BP(2)^2) - (mu_earth/r1) - (mu_moon/r2);


V= @(xv,yv) (omega^2*(xv^2+yv^2) + 2*(mu_earth/(sqrt( (xv+p2*r_em)^2+yv^2 ))) + 2*(mu_moon/(sqrt( (xv-p1*r_em)^2+yv^2 ))) + 2*C);
figure(9)
fimplicit(V,[-0.5*(10^6) 0.5*(10^6) -0.5*(10^6) 0.5*(10^6)]), hold on
plot(rB9(:,1),rB9(:,2),'r'), axis equal, grid on, hold on
plot(r_em,0,'.k','MarkerSize',10), hold on
hold off
xlabel('Position in i Direction')
ylabel('Postion in j Direction')
legend('Forbidden Region','Satellite Trajectory', 'Moon')
title('Spacecraft Trajectory in the CR3BP Frame ')

%% Part 12
time12 = linspace(0,7776000,5000)';

options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[tB12,xB12] = ode45(@CR3BP, time12,[r0_CR3BP;v0_CR3BP],options); %Spacecraft + Earth + Moon in EC 

rB12= xB12(:,1:3); %position of spacecraft relative to Earth center
figure(10)
fimplicit(V,[-0.5*(10^6) 0.5*(10^6) -0.5*(10^6) 0.5*(10^6)]), hold on
plot(rB12(:,1),rB12(:,2),'r'), axis equal, grid on,hold on
plot(r_em,0,'.k','MarkerSize',20), hold on
plot(0,0,'.g','MarkerSize',10)
hold off

xlabel('Position in i direction(km)')
ylabel('Postion in j direction (km)')
legend('Forbidden Region','Satellite Trajectory', 'Moon','Earth')
title('Spacecraft Trajectory in the CR3BP Frame (90 Days)')
