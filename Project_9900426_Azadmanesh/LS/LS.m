clear all
close all;
clc;
Ts=2;
N_samples=4000;

u=zeros(1,N_samples);
u=[ones(1,N_samples/10) zeros(1,N_samples/10) ones(1,N_samples/10) zeros(1,N_samples/10) ones(1,N_samples/10) ones(1,N_samples/10) zeros(1,N_samples/10) ones(1,N_samples/10) zeros(1,N_samples/10) ones(1,N_samples/10)];
u=2*u+3*randn(1,N_samples);
noise=0.0001*randn(1,N_samples);
num=[0.1317 0.6007 0.3582 0.03051];

den=[1 0.4749 0.496 0.1807 0.09072];
% Sampling
y(1)=0;
y(2)=0;
y(3)=0;
y(4)=0;
for temp=5:N_samples
y(temp)=-den(2)*y(temp-1)-den(3)*y(temp-2)-den(4)*y(temp-3)-den(5)*y(temp-4)+num(1)*u(temp-1)+num(2)*u(temp-2) +num(3)*u(temp-3)+num(4)*u(temp-4)+noise(temp);
end

%Identification LS
FFI=[0 0 0 0 0 0 0 0;0 0 0 0 u(1) 0 0 0;0 0 0 0 u(2) u(1) 0 0;0 0 0 0 u(3) u(2) u(1) 0];
Y=[y(1);y(2);y(3);y(4)];
for temp=5:N_samples
 FFI=[FFI;[-y(temp-1) -y(temp-2) -y(temp-3) -y(temp-4) u(temp-1) u(temp-2) u(temp-3) u(temp-4)]];
 Y=[Y;y(temp)];
end
theta_hat=(FFI'*FFI)^-1*FFI'*Y
x=(FFI'*FFI);
w=cond(x)
figure(1)
plot(0:Ts:Ts*(N_samples-1),u)
xlabel('sample time')
ylabel('sample')
title('input')
figure(2)
plot(0:Ts:Ts*(N_samples-1),y)
xlabel('sample time')
ylabel('sample')
title('output')