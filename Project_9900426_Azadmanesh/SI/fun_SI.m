function dx = fun_SI(t,x)
dx = zeros(38,1);

global mu W mc


    [U] = gravitational_potential(t,x);


%%% sliding mode control
landa = diag([.01 .01 .01]);k = [100 100 100]; epsilon = 0.8;
e = x(1:3);de = x(4:6);s = de+landa*e;
Fx = mc*(-2*W*x(6)-U(1)-landa(1)*x(4))-k(1)*satlins(s(1)/epsilon);
Fy = mc*(W^2*x(2)-U(2)-landa(2)*x(5))-k(2)*satlins(s(2)/epsilon);
Fz = mc*(2*W*x(4)-3*W^2*x(3)-U(3)-landa(3)*x(6))-k(3)*satlins(s(3)/epsilon);

%%%%%%% %%%%% spacecraft dynamics system
delta = 1e-5*sin(0.2*pi*t);  % disturbance



d2x = 2*W*x(6)+(Fx/mc+U(1)+delta);
d2y = -W^2*x(2)+(Fy/mc+U(2)+delta);
d2z = -2*W*x(4)+3*W^2*x(3)+Fz/mc+U(3)+delta;

dx(1)=x(4);    %x(1)=x
dx(2)=x(5);    %x(2)=y
dx(3)=x(6);    %x(3)=z
dx(4)=d2x;     %x(4)=xdot
dx(5)=d2y;     %x(5)=ydot
dx(6)=d2z;     %x(6)=zdot

%%%%%%% asteroid's motion relative to sun
ra = sqrt(x(7)^2+x(8)^2+x(9)^2);
d2xa = -mu*x(7)/(ra^3);
d2ya = -mu*x(8)/(ra^3);
d2za = -mu*x(9)/(ra^3);
dx(7) = x(10);
dx(8) = x(11);
dx(9) = x(12);
dx(10) = d2xa;
dx(11) = d2ya;  
dx(12) = d2za;  

%%%%%%%%%%%%%%%%%%%%

W11 = Fx;
W21 = Fy;
W31 = Fz;

%%%% gravitational potential fitting with a 2-order system        

W12 = (7.253e-11  *t^2  -1.067e-07 *t -5.709e-06);
W22 = (1.51e-12 *t^2 -1.66e-09 *t + 6.81e-07 )  ;
W32 = ( 2.151e-10*t^2 -3.006e-07*t -1.894e-05 ) ;


W13 = 2*x(6);
W23 = 0;
W33 = -2*x(4);

W14 = 0;
W24 = -x(2);
W34 = 3*x(3);


We = [W11 W12 W13 W14;W21 W22 W23 W24;W31 W32 W33 W34];   % signal vector

% phi = [1/mc 1 W W^2]';     % system's parameters

phi_hat = x(13:16);   
tau_tilda =  [d2x;d2y;d2z]-We*phi_hat;

P = [x(17) x(21) x(25) x(29);...
    x(18) x(22) x(26) x(30);...
    x(19) x(23) x(27) x(31);...
    x(20) x(24) x(28) x(32)];% Kalman Filter matrix

% landa0 = diag([.3 .01 .1 .1]); k0 = 1;
landa0 = diag([.5 0.6 .2 .2]); k0 = 1;
la = landa0*(1-norm(P)/k0);  %#forgetting factor
dp = la*P-P*(We')*We*P;
dphi_hat = P*We'*tau_tilda;

dx(13:16) = dphi_hat;
dx(17:32) = dp; 


%%%%%%%%%%%%%%%%%%%%%%%%% Estimation of system
landa = diag([.01 .01 .01]);k = [10 10 10]; epsilon = 0.8;
e = x(33:35);de = x(36:38);s = de+landa*e;
Fx = (1/x(13))*(-2*x(15)*x(36)-W12-landa(1)*x(36))-k(1)*satlins(s(1)/epsilon);
Fy = (1/x(13))*(x(15)^2*x(34)-W22-landa(2)*x(37))-k(2)*satlins(s(2)/epsilon);
Fz = (1/x(13))*(2*x(15)*x(36)-3*x(15)^2*x(35)-W32-landa(3)*x(38))-k(3)*satlins(s(3)/epsilon);


d2x_hat = 2*x(15)*x(38)+(Fx*x(13)+W12);
d2y_hat = -x(15)^2*x(34)+(Fy*x(13)+W22);
d2z_hat = -2*x(15)*x(36)+3*x(15)^2*x(35)+Fz*x(13)+W32;

dx(33)=x(36);    %x(1)=x
dx(34)=x(37);    %x(2)=y
dx(35)=x(38);    %x(3)=z
dx(36)=d2x_hat;     %x(4)=xdot
dx(37)=d2y_hat;     %x(5)=ydot
dx(38)=d2z_hat;     %x(6)=zdot




end

