function f1 = f1_mass_constrained_z_0_2(x,z)
%F1_MASS_CONSTRAINED_Z_0_2
%    F1 = F1_MASS_CONSTRAINED_Z_0_2(X1,X2)
%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-May-2020 23:57:57
% Includes scripted file modifications (LBracket.fixFileInputs) to match system specs

% unpack the block vars into the variable names used by the function
x1 = x(1, :).';
x2 = x(2, :).';

% Original content from Symbolic Math Toolbox
t2 = x1.*4.856217119292648e-1;
t3 = t2+3.237478079528432e-1;
t4 = x1.*3.523798417941783e-1;
t5 = t4+2.349198945294522e-1;
t6 = 1.0./t3;
t7 = t4+3.15856846517663e-1;
t8 = 1.0./t5;
t13 = t3.*t7.*t8;
t9 = t2-t13+1.650178343204802e-1;
t10 = x1.*3.237478079528432e-1;
t11 = t10+1.92286867573956e-1;
t12 = x1.*3.491989452945221e-2;
t14 = (t5.*t6.*t9)./6.0;
t15 = t12+t14+9.79043768826634e-3;
t16 = conj(x2);
t17 = t16.*(2.0./5.0);
t18 = t17+1.0./1.0e1;
t19 = x1.*1.618739039764216e-1;
t20 = t19+1.314609403788872e-1;
t21 = x2.*(2.0./5.0);
t22 = t21+1.0./1.0e1;
t23 = conj(x1);
t24 = t23.*4.856217119292648e-1;
t25 = t23.*3.523798417941783e-1;
t26 = t25+2.349198945294522e-1;
t27 = t24+3.237478079528432e-1;
t28 = x1.*3.174599472647261e-1;
t29 = t23.*1.618739039764216e-1;
t30 = t29+1.314609403788872e-1;
t31 = x1.*2.825400527352739e-1;
t32 = t14+t31+1.748705431588141e-1;
t33 = t22.*t32;
t34 = -t14+t28+2.251294568411859e-1;
t35 = t22.*t34;
t36 = t33+t35;
t37 = t29+3.146094037888723e-2;
t38 = t5.*t6.*t9.*(5.0./6.0);
t39 = t28+t38+2.251294568411859e-1;
t40 = t22.*t39;
t41 = t33+t40;
t42 = t19+3.146094037888723e-2;
t43 = t20.*t22;
t44 = 1.0./t27;
t45 = 1.0./t26;
t46 = t25+3.15856846517663e-1;
t52 = t27.*t45.*t46;
t47 = t24-t52+1.650178343204802e-1;
t48 = (t26.*t44.*t47)./6.0;
t49 = t11.*t22;
t50 = t43+t49;
t51 = t23.*3.174599472647261e-1;
t53 = t10+1.335568939415929e-1;
t54 = t12+t14+9.072738967647714e-2;
t55 = t23.*3.237478079528432e-1;
t56 = t55+1.92286867573956e-1;
t57 = t23.*3.491989452945221e-2;
t58 = t26.*t44.*t47.*(5.0./6.0);
t59 = t51+t58+2.251294568411859e-1;
t60 = t22.*t42;
t61 = -t48+t51+2.251294568411859e-1;
t62 = t43-t60;
t63 = t22.*t54;
t64 = t40+t63;
t65 = t20.*t32;
t66 = t20.*t39;
t67 = t22.*t53;
t68 = t60+t67;
t69 = t48+t57+9.79043768826634e-3;
t70 = t32.*t42;
t71 = t15.*t22;
t72 = t35+t71;
t73 = t48+t57+9.072738967647714e-2;
t74 = t55+1.335568939415929e-1;
t75 = t23.*2.825400527352739e-1;
t76 = t48+t75+1.748705431588141e-1;
f1 = t18.*(t66+t70).*(-1.0./6.0)+(t30.*t36)./6.0+(t30.*t41)./6.0-(t36.*t37)./6.0-(t37.*t41)./6.0+(t37.*t64)./6.0-(t30.*t72)./6.0+(t37.*t72)./6.0-(t50.*t59)./6.0+(t50.*t61)./6.0+(t56.*t64)./3.0+(t59.*t62)./3.0+(t61.*t62)./6.0-(t59.*t68)./6.0-(t62.*t69)./6.0-(t68.*t69)./3.0+(t56.*(t35-t40))./6.0-(t69.*(t43-t49))./6.0+(t74.*(t63-t71))./6.0+(t18.*(t11.*t54-t15.*t53))./6.0-(t18.*(t39.*t53-t42.*t54))./2.0-(t18.*(t15.*t20-t11.*(t28-(t5.*t6.*t9)./6.0+2.251294568411859e-1)))./2.0+(t56.*(t35-t15.*t22))./6.0+(t18.*(t65+t20.*t34))./3.0+(t18.*(t65-t32.*t42))./2.0-(t18.*(t66-t34.*t42))./6.0-(t18.*(t70+t39.*t42))./6.0+(t76.*(t43-t22.*t42))./6.0+(t73.*(t49-t22.*t53))./6.0+t22.*t30.*t32.*(2.0./3.0)-(t22.*t37.*t39)./3.0+(t20.*t22.*t61)./3.0+(t11.*t22.*t73)./3.0-(t15.*t22.*t74)./3.0+(t22.*t37.*t54)./3.0-(t22.*t53.*t59)./3.0-t22.*t42.*t76.*(2.0./3.0);

