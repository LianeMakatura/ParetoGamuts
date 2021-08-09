function hess = BikeStay_Mass_hess(in1,z1)
%BIKESTAY_MASS_HESS
%    HESS = BIKESTAY_MASS_HESS(IN1,Z1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    09-Aug-2021 18:33:04

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
t4 = x1.*6.0e1;
t5 = z1.*9.8e1;
t2 = -t4+t5+2.32e2;
t8 = x1.*3.0e1;
t3 = t8-1.3e2;
t6 = t2.^2;
t7 = t6./4.0;
t9 = t3.^2;
t10 = t9./4.0;
t11 = t7+t10;
t12 = x3.*1.0e1;
t13 = t12+6.0;
t14 = x2.*1.0e1;
t15 = t14+1.0e1;
t18 = x1.*2.25e3;
t19 = z1.*2.94e3;
t16 = -t18+t19+8.91e3;
t17 = 1.0./sqrt(t11);
t20 = 1.0./t11.^(3.0./2.0);
t21 = z1.*4.802e3;
t27 = x1.*2.94e3;
t22 = t21-t27+1.1368e4;
t23 = x1.*(7.2e1./6.29e2);
t24 = t23-(t15.*t16.*t17)./1.5725e4+1.2e2./6.29e2;
t25 = sqrt(t11);
t26 = t25.*1.271860095389507e-3;
t28 = (t13.*t15.*t16.*t20.*t22)./3.145e5;
t29 = t28-t13.*t15.*t17.*1.869634340222576e-2;
t30 = (t13.*t17.*t22)./1.5725e4;
t31 = (t15.*t17.*t22)./1.5725e4;
hess = reshape([x3.*(7.2e1./6.29e2)+t13.*t15.*t17.*(9.0./6.29e2)-(t13.*t15.*t16.^2.*t20)./3.145e5+6.868044515103339e-2,t13.*t16.*t17.*(-6.359300476947536e-5),t24,t29,t13.*t16.*t17.*(-6.359300476947536e-5),0.0,t26,t30,t24,t26,0.0,t31,t29,t30,t31,t13.*t15.*t17.*3.053736089030207e-2-(t13.*t15.*t20.*t22.^2)./3.145e5],[4,4]);