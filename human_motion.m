function human_motion
clear;clc;format compact;
t = 0:0.01:2*pi;
r1 = 0.05;
r2 = 0.1;
w = 2*pi/2;
len = length(t);
xh = r1*cos(w*t);
yh = r2*sin(w*t);
dxh = -r1*w*sin(w*t);
dyh =  r2*w*cos(w*t);
ddxh = -r1*w*w*cos(w*t);
ddyh = r2*w*w*sin(w*t);
qh = [xh ; yh ; zeros(1,len)];
dqh = [dxh ; dyh; zeros(1,len)];
ddqh = [ddxh ; ddyh; zeros(1,len)];

% subplot(131);
% plot(xh, yh, 'b' , 'Linewidth', 1); axis equal;
% subplot(132);
% plot(dxh, dyh, 'b' , 'Linewidth', 1); axis equal;
% subplot(133);
% plot(ddxh, ddyh, 'b' , 'Linewidth', 1); axis equal;

ddq = [];
tau = [];
lambda = [];
alpha = 40;
for i = 1:len
    if i == 1
[ddq_, tau_, lambda_] = myQP(ddqh(1,i), ddqh(2,i), ddqh(3,i), ...
    0, 0, 0, -35, 0);
    else
[ddq_, tau_, lambda_] = myQP(ddqh(1,i), ddqh(2,i), ddqh(3,i), ...
    ddqh(1,i)-ddqh(1,i-1), ddqh(2,i)-ddqh(2,i-1), ddqh(3,i)-ddqh(3,i-1), lambda(2,i-1), alpha);
    end
ddq = [ddq, ddq_];
tau = [tau, tau_];
lambda = [lambda, lambda_];
end

wd = 1.5;
start = 1;
subplot(221); 
plot(start:len, tau(1,:), 'r' , 'Linewidth', wd); hold on;
plot(start:len, tau(2,:), 'b' , 'Linewidth', wd);
plot(start:len, tau(3,:), 'k' , 'Linewidth', wd);
title('SuperLimb tau');
subplot(222); 
plot(start:len, tau(4,:), 'k' , 'Linewidth', wd); hold on;
plot(start:len, tau(5,:), 'k' , 'Linewidth', wd);
plot(start:len, tau(6,:), 'b' , 'Linewidth', wd);
title('Human tau');
subplot(223);
plot(start:len, lambda(1,:), 'r' ,'Linewidth', wd); hold on;
plot(start:len, lambda(2,:), 'b' ,'Linewidth', wd);
title('Supporting force');
subplot(224);
plot(start:len, ddq(1,:), 'r' ,'Linewidth', wd); hold on;
plot(start:len, ddq(2,:), 'b' ,'Linewidth', wd);
title('ddqh');
end
