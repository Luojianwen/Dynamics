function superlimb_opt
clear;clc;format compact;
global mh Ih m1 m2 m3;
global Sk Skc Sh;
global l0 l1 l2 l3 l11 l22 l33;
mh = 0.0001; Ih = 0.0001; m1 = 2.5; m2 = 1.2; m3 = 0.1;
l0 = 0.01; l1 = 0.2; l2 = 0.2; l3 = 0.05; l11 = 0.05; l22 = 0.05; l33 = 0.005;
Sk = [eye(2), zeros(2,4)]; Sh = [zeros(3),eye(3)];
Skc = [zeros(4,2), eye(4,4)];
dt = 0.001;
t = 0:dt:1*pi;
r1 = 0.001; %0.001
r2 = 0.01; %0.01
w = 2*pi/1;
len = length(t);
[xh0, yh0, thetah0] = init;
xh = xh0 + r1*(cos(w*t)-1);
yh = yh0 + r2*sin(w*t);
thetah = zeros(1,len);
dxh = -r1*w*sin(w*t);
dyh =  r2*w*cos(w*t);
% dthetah = zeros(1,len);
ddxh = -r1*w*w*cos(w*t);
ddyh = -r2*w*w*sin(w*t);
ddthetah = zeros(1,len);
ddqh = [ddxh;ddyh;ddthetah];

q0 = [45*pi/180, 90*pi/180, -45*pi/180, xh0, yh0, thetah0];
[A, h, Jc, Jhs, dot_Jhs] = kindynModel(q0, zeros(1,6));
q = [q0',zeros(6,len-1)];
dq = zeros(6,len);
ddq = zeros(6,len);
tau = zeros(6,len);
lambda = zeros(2,len);

use_QP = 0;
deepBlue = [0 0 139]/255;
maroon = [207, 33, 71]/255;
if use_QP == 0
    color = deepBlue;
else
    color = maroon;
end
for i = 1:len-1
    if i == 1
        qs_tmp = [45*pi/180, 90*pi/180, -45*pi/180]' + Jhs^-1 * [xh(i+1)-xh(i), yh(i+1)-yh(i), 0]';
        qh = [xh0, yh0, thetah0]';
        q(:,i+1) = [qs_tmp; qh];
        dqs_tmp = pinv(Jhs) * [dxh(i), dyh(i), 0]';
        dq(:,i+1) = [dqs_tmp; [dxh(i), dyh(i), 0]'];
    else
        [A, h, Jc, Jhs, dot_Jhs] = kindynModel(q(:,i)', dq(:,i)');
        qs_tmp = q(1:3,i) + 1 * pinv(Jhs) * [xh(i+1)-xh(i), yh(i+1)-yh(i), 0]';
        qh_tmp = [xh(i), yh(i), thetah(i)]';
        q(:,i+1) = [qs_tmp; qh_tmp];
        dqs_tmp = Jhs^-1 * [dxh(i), dyh(i), 0]';
        dq(:,i+1) = [dqs_tmp; [dxh(i), dyh(i), 0]'];
        ddq(1:3,i) = pinv(Jhs) * (ddqh(:,i) - dot_Jhs * dqs_tmp);
        ddq(4:6,i) = ddqh(:,i);
        
        [Q, R] = qr(Jc');
        SQ = Skc*Q';
        S_bar = dynamicalInv(A, SQ)*(SQ);
        Nkc = eye(6) - S_bar;
        R_bar = R(1:2,1:2)^-1*Sk*Q'*Nkc;
        Aqh = (A*ddq(:,i) + h);

        lambda_b = R_bar*Aqh;
        tau_ini = [-0.6, 57, 0.8]';
        W = Sh*Jc'; b = tau_ini - W*lambda_b;
        Q = diag([1e5,1e5,1e5]); H = W'*Q*W; H = (H+H')/2; 
        f = b'*Q*W;
%         if abs(Aqh - tau(:,i) - Jc'*lambda(:,i))>1e3
%             flag = 0
%         end
        if use_QP == 1
            sol = quadprog(H, f);
            tau(:,i) = Jc'*(lambda_b - sol);
        else
            tau(:,i) = Jc'*(lambda_b - [0,-20]');
        end
        lambda(:,i) = R_bar*(Aqh-tau(:,i));
    end
end


x = zeros(1,len);
y = zeros(1,len);
for i = 1:len
    [x(i), y(i)] = FK(q(:,i));
end

wd = 2;

% figure(1)
% plot(xh, yh, 'Color', color, 'Linewidth', wd); axis equal; hold on; xlabel('x_h (m)'); ylabel('y_h (m)');
n = 6;
% figure(2)
% subplot(n,1,1);
% plot(t, 180*q(1,:)/pi, 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; ylabel('qs1 (degree)');
% subplot(n,1,2);
% plot(t, 180*q(2,:)/pi, 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; ylabel('qs2 (degree)');
% subplot(n,1,3);
% plot(t, 180*q(3,:)/pi, 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; ylabel('qs3 (degree)');
% subplot(n,1,4);
% plot(t, q(4,:), 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; ylabel('xh (m)');
% % plot(t, x, 'r' , 'Linewidth', wd);
% subplot(n,1,5);
% plot(t, q(5,:), 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; ylabel('yh (m)');
% % plot(t, y, 'r' , 'Linewidth', wd);
% subplot(n,1,6);
% plot(t, 180*q(6,:)/pi, 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; xlabel('time (s)'); ylabel('thetah (degree)');

% figure(3)
% n = 6;
% subplot(n,1,1);
% plot(t, dq(1,:), 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; ylabel('dqs1 (rad/s)');
% subplot(n,1,2);
% plot(t, dq(2,:), 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; ylabel('dqs2 (rad/s)');
% subplot(n,1,3);
% plot(t, dq(3,:), 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; ylabel('dqs3 (rad/s)');
% subplot(n,1,4);
% plot(t, dq(4,:), 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; ylabel('dxh (m/s)');
% subplot(n,1,5);
% plot(t, dq(5,:), 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; ylabel('dyh (m/s)');
% subplot(n,1,6);
% plot(t, dq(6,:), 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; xlabel('time (s)'); ylabel('dthetah (rad/s)');

figure(4)
n = 6;
subplot(n,1,1);plot(t, ddq(1,:), 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; ylabel('ddqs1 (rad/s^2)');
subplot(n,1,2);plot(t, ddq(2,:), 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; ylabel('ddqs2 (rad/s^2)');
subplot(n,1,3);plot(t, ddq(3,:), 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; ylabel('ddqs3 (rad/s^2)');
subplot(n,1,4);plot(t, ddq(4,:), 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; ylabel('ddxh (m/s^2)');
subplot(n,1,5);plot(t, ddq(5,:), 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; ylabel('ddyh (m/s^2)');
subplot(n,1,6);plot(t, ddq(6,:), 'Color', color, 'Linewidth', wd); xlim([0,3]); hold on; xlabel('time (s)'); ylabel('ddthetah (rad/s^2)');




% figure(5)
% tau(:,1) = tau(:,2);tau(:,len) = tau(:,len-1);
% n = 6;
% subplot(n,1,1);plot(t, tau(1,:), 'Color', color, 'Linewidth', wd); hold on; 
% xlim([0,3]);
% title('q_{s1} torques');
% subplot(n,1,2);plot(t, tau(2,:), 'Color', color, 'Linewidth', wd); hold on; 
% xlim([0,3]);
% title('q_{s2} torques');
% subplot(n,1,3);plot(t, tau(3,:), 'Color', color, 'Linewidth', wd); hold on; 
% xlim([0,3]);
% title('q_{s3} torques');
% subplot(n,1,4);plot(t, tau(4,:), 'Color', color, 'Linewidth', wd); hold on; 
% xlim([0,3]);
% title('x_h torques');
% subplot(n,1,5);plot(t, tau(5,:), 'Color', color, 'Linewidth', wd); hold on; 
% xlim([0,3]);
% title('y_h torques');
% subplot(n,1,6);plot(t, tau(6,:), 'Color', color, 'Linewidth', wd); hold on; 
% xlim([0,3]);
% title('theta_h torques');
% 
% figure(6)
% lambda(:,1) = lambda(:,2);lambda(:,len) = lambda(:,len-1);
% n = 2;
% subplot(n,1,1);plot(t, lambda(1,:), 'Color', color, 'Linewidth', wd); hold on;title('\lambda_x'); 
% xlim([0,3]); ylim([-0.3,0.3]);
% subplot(n,1,2);plot(t, lambda(2,:), 'Color', color, 'Linewidth', wd); hold on;title('\lambda_y'); 
% xlim([0,3]); ylim([-21,-18]);

end

function [xh0, yh0, thetah0] = init
    global l1 l2 l3;
    q = [45*pi/180, 90*pi/180, -45*pi/180];
    xh0 = -( l1*cos(pi/4)+l2*cos(pi/4+q(2))+l3*cos(pi/4+q(2)+q(3)) );
    yh0 = -( l1*sin(pi/4)+l2*sin(pi/4+q(2))+l3*sin(pi/4+q(2)+q(3)) );
    thetah0 = pi/4 - q(1);
end

function [x, y] = FK(q)
    global l1 l2 l3;
    x = -( l1*cos(pi/4)+l2*cos(pi/4+q(2))+l3*cos(pi/4+q(2)+q(3)) );
    y = -( l1*sin(pi/4)+l2*sin(pi/4+q(2))+l3*sin(pi/4+q(2)+q(3)) );
end

function W_bar = dynamicalInv(A, W)
% W_mn for m > n
    A_inv = A^-1;
    W_bar = A_inv*W'*(W*A_inv*W')^-1;
end




