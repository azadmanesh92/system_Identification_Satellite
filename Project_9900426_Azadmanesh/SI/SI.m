%%%%%%%%%%%%%%%%%%%% System Identification %%%%%%%%%%
%%%% 4-April-2021

clc
clear
close all

global mu W mc G ma RB 

mc = 150;      % spacecraft's mass(kg)
ma = 6.687e15; % 433 Eros mass (kg)
RB = 16;       % reference radius of Eros(km)
i = 10.28*pi/180;     % inclination


hs = 218e6;       % height from sun (km)
Rs = 696342;      % sun radius(km)
ra = hs+Rs;       % radius from the sun to asteroid
G = 6.6742e-20;   % universal gravitational constant(km^3/kg.s^2)
ms = 1.989e30;    % sun mass (kg)
mu = G*(ma+ms);   % gravitational parameter(km^3/s^2)
W = sqrt(mu/ra^3); %#rad/s


%#initial condition
x0 = [-5 0 -15]; dx0 = [0 40e-3 0];
ra0 = [-ra*cos(i) 0 ra*sin(i)];
int_parameter = [1/130 1 1e-7 1e-14];


%%% ode solving
options = odeset('RelTol',1e-3,'AbsTol',1e-5);
int = [x0 dx0 ra0 0 -sqrt(mu/ra) 0 int_parameter 20*ones(1,4) 1*ones(1,4) ones(1,4) 10*ones(1,4) x0 dx0];

        [t,x]=ode45(@fun_SI,[0 500],int);


%# plots        
%%% estimation of spacecraft mass        
figure
plot(t,1./x(:,13),t,mc*ones(numel(t),1))
xlabel('t(s)'),ylabel('mc(kg)'),grid on
legend('estimation','desired')

%%% estimation of angular velocity       
figure
plot(t,x(:,15),t,W*ones(numel(t),1))
xlabel('t(s)'),ylabel('W(rad/s)'),grid on
legend('estimation','desired')

%%% relative distance       
figure
plot(t,sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2)),hold on
plot(t,sqrt(x(:,33).^2+x(:,34).^2+x(:,35).^2))
xlabel('t(s)'),ylabel('R(km)'),grid on
legend('estimation','desired')

%%% relative velocity        
figure
plot(t,sqrt(x(:,4).^2+x(:,5).^2+x(:,6).^2)),hold on
plot(t,sqrt(x(:,36).^2+x(:,37).^2+x(:,38).^2))
xlabel('t(s)'),grid on
h = ylabel('$ \dot{R}(km/s) $');
set(h,'Interpreter','latex')
legend('estimation','desired')
   
%% gravitational potential fitting with a 2-order system        
for i = 1:numel(t)
[U(:,i)] = gravitational_potential(t(i),x(i,:));
end

delta = 1e-5*sin(0.2*pi.*t);

figure
plot(t,sqrt((U(1,:)+delta').^2+(U(2,:)+delta').^2+(U(3,:)+delta').^2)),hold on

f=fit(t,U(1,:)'+delta,'poly2');
f=fit(t,U(2,:)'+delta,'poly2');
f=fit(t,U(3,:)'+delta,'poly2');


W12 = (7.253e-11  *t.^2  -1.067e-07 .*t -5.709e-06);
W22 = (1.51e-12 *t.^2 -1.66e-09 .*t + 6.81e-07 )  ;
W32 = ( 2.151e-10*t.^2 -3.006e-07.*t -1.894e-05 ) ;

plot(t,sqrt(W12.^2+W22.^2+W32.^2)),grid on
xlabel('t(s)'),legend('$ U+\delta $','$\hat{U}+\hat{\delta}$','Interpreter','latex')
