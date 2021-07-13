clear;
clc;
close all;

N=1000;
T=0.1;
Mc=150;
R=16e3;
Ma=6.687e15;
G=6.67e-11;
w=1.12e-7;

u1=5+2*randn(1,N);
u2=5+2*randn(1,N);
u3=5+2*randn(1,N);

for i=1:N
   u1(i)=u1(i)+0.1*i;
   u2(i)=u2(i)+0.1*i;
   u3(i)=u3(i)+0.1*i;
end
    
noise1=0.0001*randn(1,N);
noise2=0.0001*randn(1,N);
noise3=0.0001*randn(1,N);


x1=zeros(1,N);
x2=zeros(1,N);
x3=zeros(1,N);
x4=zeros(1,N);
x5=zeros(1,N);
x6=zeros(1,N);

for i=1:N-1
    x1(i+1)=x1(i)+T*x2(i);
    x2(i+1)=x2(i)+2*T*w*x6(i)+(1/Mc)*T*u1(i)+T*G*Mc*Ma/(R+x1(i))^2+noise1(i);
    x3(i+1)=x3(i)+T*x4(i);
    x4(i+1)=x4(i)-T*w^2*x3(i)+(1/Mc)*T*u2(i)+T*G*Mc*Ma/(R+x3(i))^2+noise2(i);
    x5(i+1)=x5(i)+T*x6(i);
    x6(i+1)=x6(i)-2*T*w*x2(i)+3*T*w^2*x5(i)+(1/Mc)*T*u3(i)+T*G*Mc*Ma/(R+x5(i))^2+noise3(i); 
    
end    

U=[u1' u2' u3'];
Y=[x1' x3' x5'];