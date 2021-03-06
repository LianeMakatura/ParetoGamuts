function deriv = BikeRocker_deriv(in1,z1)
%BIKEROCKER_DERIV
%    DERIV = BIKEROCKER_DERIV(IN1,Z1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    09-Aug-2021 18:16:55

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
t2 = x3.*5.0;
t3 = t2+5.0;
t5 = z1.*4.9e1;
t4 = t5-9.0;
t6 = t5+4.1e1;
t20 = x1.*2.0e1;
t7 = t20+2.5e1;
t8 = t4.^2;
t9 = t8+1.0e4;
t10 = sqrt(t9);
t11 = x2.*7.0;
t12 = t11+7.0;
t13 = t6.^2;
t14 = t13+4.444444444444444e3;
t15 = sqrt(t14);
t16 = sqrt(1.3e1);
t17 = z1.*4.802e3;
t18 = x1.*8.0e2;
t19 = t18+1.0e3;
t21 = t7.^2;
t22 = 1.0./t3;
t23 = t10.*7.0;
t24 = t15.*7.0;
t25 = t16.*(3.5e2./3.0);
t26 = t23+t24+t25;
t27 = t10.*t12;
t28 = t12.*t15;
t29 = t12.*t16.*(5.0e1./3.0);
t30 = t21+t27+t28+t29;
t31 = 1.0./t30.^2;
t32 = t17-8.82e2;
t33 = 1.0./sqrt(t9);
t34 = (t12.*t32.*t33)./2.0;
t35 = t17+4.018e3;
t36 = 1.0./sqrt(t14);
t37 = (t12.*t35.*t36)./2.0;
t38 = t34+t37;
t39 = 1.0./t30;
deriv = reshape([(t3.*t19)./5.3111e4,t10.*t19.*t22.*t31.*(-1.176470588235294e2),(t3.*t26)./5.3111e4,t10.*t22.*t26.*t31.*(-1.176470588235294e2),t21.*9.414245636497147e-5+t10.*t12.*9.414245636497147e-5+t12.*t15.*9.414245636497147e-5+t12.*t16.*1.569040939416191e-3,1.0./t3.^2.*t10.*t39.*(-5.882352941176471e2),(t3.*t38)./5.3111e4,t10.*t22.*t31.*t38.*(-1.176470588235294e2)+t22.*t32.*t33.*t39.*(1.0e3./1.7e1)],[2,4]);
