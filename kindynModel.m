function [A, h, Jc, Jhs, dot_Jhs] = kindynModel(q, dot_q)
% q is the generalized coordinates
% q = [theta1, theta2, theta3, xh, yh, thetah]^T;
% dot_q = [dtheta1, dtheta2, dtheta3, dxh, dyh, dthetah]^T;
global mh Ih m1 m2 m3;
global l0 l1 l2 l3 l11 l22 l33;

% mh = 0.0001; Ih = 0.0001; m1 = 2.5; m2 = 1.2; m3 = 0.1;
% l0 = 0.01; l1 = 0.2; l2 = 0.2; l3 = 0.05; l11 = 0.05; l22 = 0.05; l33 = 0.005;

sh    = sin(q(6));
s1h   = sin(q(1)+q(6));
s12h  = sin(q(1)+q(2)+q(6));
s123h = sin(q(1)+q(2)+q(3)+q(6));
ch    = cos(q(6));
c1h   = cos(q(1)+q(6));
c12h  = cos(q(1)+q(2)+q(6));
c123h = cos(q(1)+q(2)+q(3)+q(6));

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

    dt1 = dot_q(1); dt2 = dot_q(2); dt3 = dot_q(3); dt4 = dot_q(4); dt5 = dot_q(5); dth = dot_q(6);
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
    
    h = A_dot*dot_q' ...
        - 0.5*[dot_q*dAt1*dot_q', dot_q*dAt2*dot_q', dot_q*dAt3*dot_q', ...
        dot_q*dAt4*dot_q', dot_q*dAt5*dot_q', dot_q*dAt6*dot_q']' + g;
    
    c2pi4  = cos(q(2)+pi/4);
    c23pi4 = cos(q(2)+q(3)+pi/4);
    s2pi4  = sin(q(2)+pi/4);
    s23pi4 = sin(q(2)+q(3)+pi/4);
    
    Jhs = [...
        0,  l2*s2pi4+l3*s23pi4,  l3*s23pi4;
        0, -l2*c2pi4-l3*c23pi4, -l3*c23pi4;
       -1,                   0,          0;
    ];
    dot_Jhs = [...
        0, l2*c2pi4+2*l3*c23pi4, 2*l3*c23pi4;
        0, l2*s2pi4+2*l3*s23pi4, 2*l3*s23pi4;
        0,                    0,           0;
    ];
end
