function f1 = f1_mass_constrained_z_0(x,z)
%F1_MASS_CONSTRAINED_Z_0
%    F1 = F1_MASS_CONSTRAINED_Z_0(X1,X2)
%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-May-2020 23:57:32
% Includes scripted file modifications (LBracket.fixFileInputs) to match system specs

% unpack the block vars into the variable names used by the function
x1 = x(1, :).';
x2 = x(2, :).';

% Original content from Symbolic Math Toolbox
t2 = x1.*4.246881205779131e-1;
t3 = t2+2.831254137186087e-1;
t4 = x1.*4.238395925819106e-1;
t5 = t4+2.825597283879404e-1;
t6 = x1.*1.415627068593044e-1;
t7 = x1.*2.587201358060298e-1;
t8 = 1.0./t3;
t9 = 1.0./t5;
t10 = t4+3.533410818175926e-1;
t14 = t3.*t9.*t10;
t11 = t2-t14+1.124854816216236e-1;
t12 = (t5.*t8.*t11)./6.0;
t13 = t7+t12+1.606831982990778e-1;
t15 = x2.*(2.0./5.0);
t16 = t15+1.0./1.0e1;
t17 = x1.*3.412798641939702e-1;
t18 = t5.*t8.*t11.*(5.0./6.0);
t19 = t17+t18+2.393168017009222e-1;
t20 = t16.*t19;
t21 = x1.*8.255972838794041e-2;
t22 = t12+t21+1.140242801166704e-1;
t23 = t16.*t22;
t24 = t20+t23;
t25 = conj(x1);
t26 = t6+1.593514922337206e-2;
t27 = t25.*4.246881205779131e-1;
t28 = t25.*4.238395925819106e-1;
t29 = t28+2.825597283879404e-1;
t30 = t27+2.831254137186087e-1;
t31 = t25.*2.831254137186087e-1;
t32 = t31+1.671902644952367e-1;
t33 = -t12+t17+2.393168017009222e-1;
t34 = t16.*t33;
t35 = t12+t21+4.324292668701824e-2;
t36 = t25.*1.415627068593044e-1;
t37 = t36+1.593514922337206e-2;
t38 = t16.*t35;
t39 = t34+t38;
t40 = t6+1.159351492233721e-1;
t41 = conj(x2);
t42 = t41.*(2.0./5.0);
t43 = t42+1.0./1.0e1;
t44 = t13.*t16;
t45 = t34+t44;
t46 = t36+1.159351492233721e-1;
t47 = t20+t44;
t48 = x1.*2.831254137186087e-1;
t49 = t25.*3.412798641939702e-1;
t50 = 1.0./t30;
t51 = 1.0./t29;
t52 = t28+3.533410818175926e-1;
t58 = t30.*t51.*t52;
t53 = t27-t58+1.124854816216236e-1;
t54 = t48+1.671902644952367e-1;
t55 = t16.*t54;
t56 = t16.*t40;
t57 = t55+t56;
t59 = t29.*t50.*t53.*(5.0./6.0);
t60 = t49+t59+2.393168017009222e-1;
t61 = t16.*t26;
t65 = (t29.*t50.*t53)./6.0;
t62 = t49-t65+2.393168017009222e-1;
t63 = t48+9.655033239825157e-2;
t64 = t16.*t63;
t66 = t61+t64;
t67 = t25.*8.255972838794041e-2;
t68 = t65+t67+4.324292668701824e-2;
t69 = t13.*t26;
t70 = t19.*t40;
t71 = t56-t61;
t72 = t31+9.655033239825157e-2;
t73 = t65+t67+1.140242801166704e-1;
t74 = t25.*2.587201358060298e-1;
t75 = t65+t74+1.606831982990778e-1;
f1 = t43.*(t69+t70).*(-1.0./6.0)+(t24.*t32)./3.0+(t24.*t37)./6.0+(t37.*t39)./6.0-(t37.*t45)./6.0-(t37.*t47)./6.0-(t39.*t46)./6.0+(t45.*t46)./6.0+(t46.*t47)./6.0-(t57.*t60)./6.0+(t57.*t62)./6.0-(t60.*t66)./6.0+(t60.*t71)./3.0+(t62.*t71)./6.0-(t66.*t68)./3.0-(t68.*t71)./6.0-(t32.*(t20-t34))./6.0+(t72.*(t23-t38))./6.0+(t68.*(t55-t56))./6.0+(t73.*(t55-t64))./6.0+(t75.*(t56-t61))./6.0+(t43.*(t13.*t40+t33.*t40))./3.0+(t43.*(t22.*t26-t19.*t63))./2.0-(t43.*(t35.*t40-t33.*t54))./2.0+(t43.*(t22.*t54-t35.*t63))./6.0+(t32.*(t34-t16.*t35))./6.0-(t43.*(t69+t19.*t26))./6.0-(t43.*(t69-t13.*t40))./2.0-(t43.*(t70-t26.*t33))./6.0-(t16.*t19.*t37)./3.0+t13.*t16.*t46.*(2.0./3.0)+(t16.*t22.*t37)./3.0-t16.*t26.*t75.*(2.0./3.0)+(t16.*t40.*t62)./3.0-(t16.*t35.*t72)./3.0-(t16.*t60.*t63)./3.0+(t16.*t54.*t73)./3.0;

