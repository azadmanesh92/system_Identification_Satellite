function [U] = gravitational_potential(~,x)

global RB G ma 
r = sqrt(x(1)^2+x(2)^2+x(3)^2)+RB;
phi = atan2(sqrt(x(1)^2+x(2)^2),x(3));  %# latitude angle
landa = atan2(x(2),x(1));   %# Longitudinal angle

%# Legendre coefficients
P00 = 1; P10 = sin(phi); P11 = cos(phi);
P20 = 0.5*(3*sin(phi)^2-1); P21 = 3*sin(phi)*cos(phi); P22 = 3*cos(phi)^2*sin(phi);
P30 = 0.5*(5*sin(phi)^3-3*sin(phi)); P31 = 0.5*cos(phi)*(15*sin(phi)^2-3); P32 = 15*cos(phi)^2*sin(phi);
P_sinphi = [P00 P10 P11;P20 P21 P22;P30 P31 P32];

%# Stokes' coefficients
C = [0 0.001175 -0.000348;-0.052851 0.000102 0.083203; -0.001747 0.004083 0.002129];
S = [0 0 0.000088;0 0.000012 -0.028023; 0 0.003404 -0.000836];


    U_bar = zeros(1,1);
%# gravitational potential function
for n = 2:3
    for m = 1:n
        U_bar = U_bar +(RB/(r))^n*P_sinphi(n,m)+C(n,m)*cos(m*landa)+S(n,m)*sin(m*landa);
    end
end

        U = (G*ma/r)*(1+U_bar);


ux = U*sin(phi)*cos(landa);
uy = U*sin(phi)*sin(landa);
uz = U*cos(phi);
U = [ux uy uz];


end


