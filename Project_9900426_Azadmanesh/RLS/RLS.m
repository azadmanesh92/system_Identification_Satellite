clear all
close all;
clc;
num=[0.1317 0.6007 0.3582 0.03051];
den=[1 0.4749 0.496 0.1807 0.09072];
N_samples=100;
Theta1=den(2)*ones(N_samples);
Theta2=den(3)*ones(N_samples);
Theta3=den(4)*ones(N_samples);
Theta4=den(5)*ones(N_samples);
Theta5=num(1)*ones(N_samples);
Theta6=num(2)*ones(N_samples);
Theta7=num(3)*ones(N_samples);
Theta8=num(4)*ones(N_samples);
noise=0.003*randn(1,N_samples);

u=zeros(1,N_samples);
u=[ones(1,N_samples/10) zeros(1,N_samples/10) ones(1,N_samples/10) zeros(1,N_samples/10) ones(1,N_samples/10) ones(1,N_samples/10) zeros(1,N_samples/10) ones(1,N_samples/10) zeros(1,N_samples/10) ones(1,N_samples/10)];


y(1)=0;
y(2)=0;
y(3)=0;
y(4)=0;
for temp=5:N_samples
y(temp)=-den(2)*y(temp-1)-den(3)*y(temp-2)-den(4)*y(temp-3)-den(5)*y(temp-4)+num(1)*u(temp-1)+num(2)*u(temp-2)+num(3)*u(temp-3)+num(4)*u(temp-4)+noise(temp);
end
P_pre=100000*eye(8);
theta_hat(1:8,1:4)=0;
for temp=5:N_samples
ffiT=[-y(temp-1) -y(temp-2) -y(temp-3) -y(temp-4) u(temp-1) u(temp-2) u(temp-3) u(temp-4)];
P=P_pre-P_pre*ffiT'*(eye(1)+ffiT*P_pre*ffiT')^(-1)*ffiT*P_pre;
K(1:8,temp)=P*ffiT';
theta_hat(1:8,temp)=theta_hat(1:8,temp-1)+K(1:8,temp)*(y(temp)-ffiT*theta_hat(1:8,temp-1));
P_pre=P;
end

figure(1)
plot(theta_hat(1,:),'r','LineWidth',3)
hold on
plot(Theta1,'g--','LineWidth',2)
legend('Estimated Theta 1','Desired Theta 1','Location','best')
title('Theta 1')

figure(2)
plot(theta_hat(2,:),'r','LineWidth',2)
hold on
plot(Theta2,'g--','LineWidth',2)
legend('Estimated Theta 2','Desired Theta 2','Location','best')
title('Theta 2')

figure(3)
plot(theta_hat(3,:),'r','LineWidth',2)
hold on
plot(Theta3,'g--','LineWidth',2)
legend('Estimated Theta 3','Desired Theta 3','Location','best')
title('Theta 3')

figure(4)
plot(theta_hat(4,:),'r','LineWidth',2)
hold on
plot(Theta4,'g--','LineWidth',2)
legend('Estimated Theta 4','Desired Theta 4','Location','best')
title('Theta 4')

figure(5)
plot(theta_hat(5,:),'r','LineWidth',2)
hold on
plot(Theta5,'g--','LineWidth',2)
legend('Estimated Theta 5','Desired Theta 5','Location','best')
title('Theta 5')

figure(6)
plot(theta_hat(6,:),'r','LineWidth',2)
hold on
plot(Theta6,'g--','LineWidth',2)
legend('Estimated Theta 6','Desired Theta 6','Location','best')
title('Theta 6')

figure(7)
plot(theta_hat(7,:),'r','LineWidth',2)
hold on
plot(Theta7,'g--','LineWidth',2)
legend('Estimated Theta 7','Desired Theta 7','Location','best')
title('Theta 7')

figure(8)
plot(theta_hat(8,:),'r','LineWidth',2)
hold on
plot(Theta8,'g--','LineWidth',2)
legend('Estimated Theta 8','Desired Theta 8','Location','best')
title('Theta 8')

figure(9)
plot(0:2:2*(N_samples-1),u)
xlabel('sample time')
ylabel('sample')
title('input')
figure(10)
plot(0:2:2*(N_samples-1),y)
xlabel('sample time')
ylabel('sample')
title('output')


