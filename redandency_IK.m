function redandency_IK
clear;clc;format compact;
the10 = 15*pi/180;
the20 = 45*pi/180;
the30 = 45*pi/180;
the40 = 45*pi/180;
[x0, y0] = FK(the10, the20, the30, the40);
len = 1000;
t = linspace(0,10*pi,len);
the1 = the10 + zeros(len);
the2 = the20 + zeros(len);
the3 = the30 + zeros(len);
the4 = the40 + zeros(len);
r1 = 0.1;
r2 = 0.2;
x = x0 - r1 + r1*cos(t);
y = y0 + r2*sin(t);
xa = zeros(len);
ya = zeros(len);
for i = 1:len-1
    J = getJacobian(the1(i), the2(i), the3(i), the4(i));
    dx = x(i+1) - x(i);
    dy = y(i+1) - y(i);
    dthe = pinv(J) * [dx, dy]';
    the1(i+1) = the1(i) + dthe(1);
    the2(i+1) = the2(i) + dthe(2);
    the3(i+1) = the3(i) + dthe(3);
    the4(i+1) = the4(i) + dthe(4);
    [xa(i), ya(i)] = FK(the1(i), the2(i), the3(i), the4(i));
end

plot(x, y, 'b-','Linewidth',1.5); axis equal; hold on;
plot(xa(1:len-1), ya(1:len-1), 'r','Linewidth',0.5); hold off;
end

function J = getJacobian(the1, the2, the3, the4)
    s1    = sin(the1);
    s12   = sin(the1+the2);
    s123  = sin(the1+the2+the3);
    s1234 = sin(the1+the2+the3+the4);
    c1    = cos(the1);
    c12   = cos(the1+the2);
    c123  = cos(the1+the2+the3);
    c1234 = cos(the1+the2+the3+the4);
    
    J = [...
        -s1-s12-s123-s1234, -s12-s123-s1234, -s123-s1234, -s1234;...
         c1+c12+c123+c1234,  c12+c123+c1234,  c123+c1234,  c1234;
    ];
end

function [x, y] = FK(the1, the2, the3, the4)
    s1    = sin(the1);
    s12   = sin(the1+the2);
    s123  = sin(the1+the2+the3);
    s1234 = sin(the1+the2+the3+the4);
    c1    = cos(the1);
    c12   = cos(the1+the2);
    c123  = cos(the1+the2+the3);
    c1234 = cos(the1+the2+the3+the4);
    
    x = c1 + c12 + c123 +c1234;
    y = s1 + s12 + s123 +s1234;
end
