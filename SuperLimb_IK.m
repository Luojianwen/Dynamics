function SuperLimb_IK
clear;clc;format compact;
global l1 l2 l3
l1 = 0.1;
l2 = 0.15;
l3 = 0.04;
theta10=45*pi/180;
theta20=90*pi/180;
theta30=-45*pi/180;
thetah0=0*pi/180;
[xh0, yh0]=FK(theta10, theta20, theta30, thetah0)
x_upper = xh0;
x_lower = xh0;
y_upper = yh0;
y_lower = yh0;
for theta2=45*pi/180:0.1:135*pi/180
    for theta3 = -70*pi/180:0.1:-20*pi/180
        [x_tmp, y_tmp] = FK(45*pi/180, theta2, theta3, 0);
        x_upper = max([x_upper, x_tmp]);
        x_lower = min([x_lower, x_tmp]);
        y_upper = max([y_upper, y_tmp]);
        y_lower = min([y_lower, y_tmp]);
    end
end
x_upper - x_lower
y_upper - y_lower
x_upper
x_lower
y_upper
y_lower
% x_upper =  0.118905680009684
% x_lower = -0.084391483851682
% y_upper = -0.046040348036566
% y_lower = -0.260662538624126

len = 1000;
t = linspace(0,4*6.28,len);
theta2 = theta20 + zeros(len);
theta3 = theta30 + zeros(len);
d_theta = zeros(2,len);
xhd = xh0 + 0.035*(cos(t)-1);  % 0.035 -0.066 ~ 0
yhd = yh0 + 0.023*sin(t);      % -0.217 -0.02 ~ 0.02
xha = zeros(len);
yha = zeros(len);
for i=1:len-1
    Jhs = computeJhs(theta2(i), theta3(i));
    d_theta_ = (Jhs^-1)*[xhd(i+1)-xhd(i), yhd(i+1)-yhd(i), 0]';
    d_theta(1,i) = d_theta_(2);
    d_theta(2,i) = d_theta_(3);
%     if(norm(d_theta_(2))>0.01)
%         180*theta2(i)/pi
%         180*theta3(i)/pi
%         Jhs^-1
%     end
    theta2(i+1) = theta2(i) + d_theta(1,i);
    theta3(i+1) = theta3(i) + d_theta(2,i);
    [xha(i), yha(i)] = FK(theta10, theta2(i), theta3(i), thetah0);
end
figure(1);
plot(t(1:len-1), 180*theta2(1:len-1)/pi, 'b','Linewidth',1); hold on;
plot(t(1:len-1), 180*theta3(1:len-1)/pi, 'r','Linewidth',1); hold off;
figure(2);
plot(xhd, yhd, 'b','Linewidth',1); axis equal; hold on;
plot(xha(1:len-1), yha(1:len-1), 'r','Linewidth',1); hold off;grid on;
figure(3);
plot(t(1:len-1), d_theta(1,1:len-1), 'b','Linewidth',0.5); hold on;
plot(t(1:len-1), d_theta(2,1:len-1), 'r-','Linewidth',1.5); hold off;
end

function Jhs = computeJhs(theta2, theta3)
global l2 l3
    s2  = sin(theta2+pi/4);
    s23 = sin(theta2+theta3+pi/4); 
    c2  = cos(theta2+pi/4);
    c23 = cos(theta2+theta3+pi/4);   
    Jhs = [...
        0,  l2*s2+l3*s23,  l3*s23;...
        0, -l2*c2-l3*c23, -l3*c23;...
        -1,0,0
    ];
end

function [xh, yh]= FK(theta1, theta2, theta3, thetah)
global l1 l2 l3
    s1h   = sin(theta1+thetah);
    s12h  = sin(theta1+theta2+thetah);
    s123h = sin(theta1+theta2+theta3+thetah);     
    c1h   = cos(theta1+thetah);
    c12h  = cos(theta1+theta2+thetah);
    c123h = cos(theta1+theta2+theta3+thetah);
    xh = -(l1*c1h+l2*c12h+l3*c123h);
    yh = -(l1*s1h+l2*s12h+l3*s123h);
end
