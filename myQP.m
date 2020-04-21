function [ddq, tau, lambda] = myQP(ddqh1, ddqh2, ddqh3, delta_ddqh1, delta_ddqh2, delta_ddqh3, lambda_y, alpha)
% clear;clc;format compact;
% Q = eye(3);
% fun = @(x) x^2 + 10;
% X = fmincon(fun, 0,[],[])
global mh Ih m1 m2 m3;
global l0 l1 l2 l3 l11 l22 l33;
global Sk Skplus Sh;
global tau0 dtau lambda_min lambda_max;
mh = 0.0001; Ih = 0.0001; m1 = 2.5; m2 = 1.2; m3 = 0.1;
l0 = 0.01; l1 = 0.2; l2 = 0.2; l3 = 0.05; l11 = 0.05; l22 = 0.05; l33 = 0.005;
Sk = [eye(2), zeros(2,4)]; Sh = [zeros(3,3),eye(3),zeros(3,6)];
Skplus = [zeros(4,2), eye(4,4)];
tau0 = 0; dtau = 30;

q = [60*pi/180, 120*pi/180, -60*pi/180, 0, 0, 0];
dotq = [0.0 0.0 0.0 0.0 0.0 0.0];
[A, h, Jc] = computeAh(q, dotq);
[Q, R] = qr(Jc');

SQ = Skplus*Q';
S_bar = dynamicalInv(A, SQ)*(SQ);
Nkplus = eye(6) - S_bar;
R_bar = R(1:2,1:2)^-1*Sk*Q'*Nkplus;
W1 = [S_bar*A, Nkplus];
W2 = [R_bar*A, -R_bar];
W = [W1;-W1;W2;-W2];
H = [1000*eye(6), zeros(6,6);zeros(6,6), 2000*eye(6)]; f=zeros(12,1);

lambda_ylim = alpha*(Jc*A^-1*Jc')^-1*Jc*[0,0,0,delta_ddqh1, delta_ddqh2, delta_ddqh3]';
lambda_min_bound = 1;
lambda_min = [-lambda_min_bound, -100]'; lambda_max = [lambda_min_bound,-35]';
lambda_max(2) = lambda_y + lambda_ylim(2);

b = [-S_bar*h+tau0+dtau; S_bar*h-tau0+dtau; -R_bar*h+lambda_max; R_bar*h-lambda_min];

% tic
x = quadprog(H, f, W, b, Sh, [ddqh1, ddqh2, ddqh3]');
ddq = x(1:6);
tau = (W1 * x + S_bar*h);
lambda = (W2 * x + R_bar*h);
% toc

% ddq = x(1:6)'
% delta_tau = x(7:12)'
% lambda = (W2 * x + R_bar*h)'
% tau = (W1 * x + S_bar*h)'
% (m1+m2+m3+mh)*9.8
% A*ddq'+h - tau' - Jc'*lambda'

% wd = 1;
% figure(1); 
% plot(1:len, tau(1,:), 'r' , 'Linewidth', wd); hold on;
% plot(1:len, tau(2,:), 'b' , 'Linewidth', wd);
% plot(1:len, tau(3,:), 'k' , 'Linewidth', wd);
% title('SuperLimb tau');
% figure(2); 
% plot(1:len, tau(4,:), 'r' , 'Linewidth', wd); hold on;
% plot(1:len, tau(5,:), 'b' , 'Linewidth', wd);
% plot(1:len, tau(6,:), 'k' , 'Linewidth', wd);
% title('Human tau');
% figure(3)
% plot(1:len, lambda(1,:), 'r' ,'Linewidth', wd); hold on;
% plot(1:len, lambda(2,:), 'b' ,'Linewidth', wd);
end

function W_bar = dynamicalInv(A, W)
% W_mn for m > n
    A_inv = A^-1;
    W_bar = A_inv*W'*(W*A_inv*W')^-1;
end

function [A, h, Jc] = computeAh(t, t_dot)
% t is the generalized coordinates
% t = [theta1, theta2, theta3, xh, yh, thetah]^T;
% t_dot = [dtheta1, dtheta2, dtheta3, dxh, dyh, dthetah]^T;
global mh Ih m1 m2 m3;
global l0 l1 l2 l3 l11 l22 l33;
sh    = sin(t(6));
s1h   = sin(t(1)+t(6));
s12h  = sin(t(1)+t(2)+t(6));
s123h = sin(t(1)+t(2)+t(3)+t(6));
ch    = cos(t(6));
c1h   = cos(t(1)+t(6));
c12h  = cos(t(1)+t(2)+t(6));
c123h = cos(t(1)+t(2)+t(3)+t(6));

    Ah = diag([0,0,0,mh,mh,Ih]);
    J1 = [...
        -l11*s1h 0 0 1 0 -l0*sh-l11*s1h;
         l11*c1h 0 0 0 1  l0*ch+l11*c1h
    ];
    J2 = [...
        -l1*s1h-l22*s12h -l22*s12h 0 1 0 -l0*sh-l1*s1h-l22*s12h;
         l1*c1h+l22*c12h  l22*c12h 0 0 1  l0*ch+l1*c1h+l22*c12h
    ];
    J3 = [...
        -l1*s1h-l2*s12h-l33*s123h -l2*s12h-l33*s123h -l33*s123h 1 0 -l0*sh-l1*s1h-l2*s12h-l33*s123h;
         l1*c1h+l2*c12h+l33*c123h  l2*c12h+l33*c123h  l33*c123h 0 1  l0*ch+l1*c1h+l2*c12h+l33*c123h
    ];
    A = Ah + J1'* diag([m1,m1])*J1 + J2'*diag([m2,m2])*J2 + J3'*diag([m3,m3])*J3;
    Jc = [...
        -l1*s1h-l2*s12h-l3*s123h -l2*s12h-l3*s123h -l3*s123h 1 0 -l0*sh-l1*s1h-l2*s12h-l3*s123h;
         l1*c1h+l2*c12h+l3*c123h  l2*c12h+l3*c123h  l3*c123h 0 1  l0*ch+l1*c1h+l2*c12h+l3*c123h
    ];

    dt1 = t_dot(1); dt2 = t_dot(2); dt3 = t_dot(3); dt4 = t_dot(4); dt5 = t_dot(5); dth = t_dot(6);
    J1_dot = [...
        -l11*c1h*(dt1+dth) 0 0 0 0 -l0*ch*dth-l11*c1h*(dt1+dth);
        -l11*s1h*(dt1+dth) 0 0 0 0 -l0*sh*dth-l11*s1h*(dt1+dth);
    ];

    J2_dot = [...
        -l1*c1h*(dt1+dth)-l22*c12h*(dt1+dt2+dth), -l22*c12h*(dt1+dt2+dth) 0 0 0 -l0*ch*dth-l1*c1h*(dt1+dth)-l22*c12h*(dt1+dt2+dth);
        -l1*s1h*(dt1+dth)-l22*s12h*(dt1+dt2+dth), -l22*s12h*(dt1+dt2+dth) 0 0 0 -l0*sh*dth-l1*s1h*(dt1+dth)-l22*s12h*(dt1+dt2+dth);
    ];
    J3_dot = [...
        -l1*c1h*(dt1+dth)-l2*c12h*(dt1+dt2+dth)-l33*c123h*(dt1+dt2+dt3+dth), ...
        -l2*c12h*(dt1+dt2+dth)-l33*c123h*(dt1+dt2+dt3+dth), -l33*c123h*(dt1+dt2+dt3+dth) 0 0 ...
        -l0*ch*dth-l1*c1h*(dt1+dth)-l2*c12h*(dt1+dt2+dth)-l33*c123h*(dt1+dt2+dt3+dth);
        -l1*s1h*(dt1+dth)-l2*s12h*(dt1+dt2+dth)-l33*s123h*(dt1+dt2+dt3+dth), ...
        -l2*s12h*(dt1+dt2+dth)-l33*s123h*(dt1+dt2+dt3+dth), -l33*s123h*(dt1+dt2+dt3+dth) 0 0 ...
        -l0*sh*dth-l1*s1h*(dt1+dth)-l2*s12h*(dt1+dt2+dth)-l33*s123h*(dt1+dt2+dt3+dth);
    ];

    A_dot = m1*(J1_dot'*J1 + J1'*J1_dot) + m2*(J2_dot'*J2 + J2'*J2_dot) + m3*(J3_dot'*J3 + J3'*J3_dot);
    
    J1t1 = [-l11*c1h 0 0 0 0 -l11*c1h; -l11*s1h 0 0 0 0 -l11*s1h];
    J1t2 = [0 0 0 0 0 0;0 0 0 0 0 0];
    J1t3 = [0 0 0 0 0 0;0 0 0 0 0 0];
    J1t4 = [0 0 0 0 0 0;0 0 0 0 0 0];
    J1t5 = [0 0 0 0 0 0;0 0 0 0 0 0];
    J1t6 = [-l11*c1h 0 0 0 0 -l0*ch-l11*c1h; -l11*s1h 0 0 0 0 -l0*sh-l11*s1h];
    
    J2t1 = [-l1*c1h-l22*c12h, -l22*c12h 0 0 0 -l1*c1h-l22*c12h; -l1*s1h-l22*s12h, -l22*s12h 0 0 0 -l1*s1h-l22*s12h];
    J2t2 = [-l22*c12h, -l22*c12h 0 0 0 -l22*c12h; -l22*s12h, -l22*s12h 0 0 0 -l22*s12h];
    J2t3 = [0 0 0 0 0 0;0 0 0 0 0 0];
    J2t4 = [0 0 0 0 0 0;0 0 0 0 0 0];
    J2t5 = [0 0 0 0 0 0;0 0 0 0 0 0];
    J2t6 = [-l1*c1h-l22*c12h, -l22*c12h 0 0 0 -l0*ch-l1*c1h-l22*c12h; ...
        -l1*s1h-l22*s12h, -l22*s12h 0 0 0 -l0*sh-l1*s1h-l22*s12h];
    
    J3t1 = [...
        -l1*c1h-l2*c12h-l33*c123h, -l2*c12h-l33*c123h, -l33*c123h 0 0 -l1*c1h-l2*c12h-l33*c123h;
        -l1*s1h-l2*s12h-l33*s123h, -l2*s12h-l33*s123h, -l33*s123h 0 0 -l1*s1h-l2*s12h-l33*c123h];
    J3t2 = [...
        -l2*c12h-l33*c123h, -l2*c12h-l33*c123h, -l33*c123h 0 0 -l2*c12h-l33*c123h;
        -l2*s12h-l33*s123h, -l2*s12h-l33*s123h, -l33*s123h 0 0 -l2*s12h-l33*c123h];
    J3t3 = [...
        -l33*c123h, -l33*c123h, -l33*c123h 0 0 -l33*c123h;
        -l33*s123h, -l33*s123h, -l33*s123h 0 0 -l33*c123h];
    J3t4 = [0 0 0 0 0 0;0 0 0 0 0 0];
    J3t5 = [0 0 0 0 0 0;0 0 0 0 0 0];
    J3t6 = [...
        -l1*c1h-l2*c12h-l33*c123h, -l2*c12h-l33*c123h, -l33*c123h 0 0 -l0*ch-l1*c1h-l2*c12h-l33*c123h;
        -l1*s1h-l2*s12h-l33*s123h, -l2*s12h-l33*s123h, -l33*s123h 0 0 -l0*sh-l1*s1h-l2*s12h-l33*c123h];
    
    dAt1 = m1*(J1t1'*J1+J1'*J1t1) + m2*(J2t1'*J2+J2'*J2t1) + m3*(J3t1'*J3+J3'*J3t1);
    dAt2 = m1*(J1t2'*J1+J1'*J1t2) + m2*(J2t2'*J2+J2'*J2t2) + m3*(J3t2'*J3+J3'*J3t2);
    dAt3 = m1*(J1t3'*J1+J1'*J1t3) + m2*(J2t3'*J2+J2'*J2t3) + m3*(J3t3'*J3+J3'*J3t3);
    dAt4 = m1*(J1t4'*J1+J1'*J1t4) + m2*(J2t4'*J2+J2'*J2t4) + m3*(J3t4'*J3+J3'*J3t4);
    dAt5 = m1*(J1t5'*J1+J1'*J1t5) + m2*(J2t5'*J2+J2'*J2t5) + m3*(J3t5'*J3+J3'*J3t5);
    dAt6 = m1*(J1t6'*J1+J1'*J1t6) + m2*(J2t6'*J2+J2'*J2t6) + m3*(J3t6'*J3+J3'*J3t6);
    
    g = 9.8*[...
        (m1*l11+m2*l1+m3*l1)*c1h+(m2*l22+m3*l2)*c12h+m3*l33*c123h;
        (m2*l22+m3*l2)*c12h+m3*l33*c123h;
        m3*l33*c123h;
        0;
        (m1+m2+m3+mh);
        (m1+m2+m3)*l0*ch+(m1*l11+m2*l1+m3*l1)*c1h+(m2*l22+m3*l2)*s12h+m3*l33*c123h;
    ];
    
    h = A_dot*t_dot' ...
        - 0.5*[t_dot*dAt1*t_dot', t_dot*dAt2*t_dot', t_dot*dAt3*t_dot', ...
        t_dot*dAt4*t_dot', t_dot*dAt5*t_dot', t_dot*dAt6*t_dot']' + g;
end


