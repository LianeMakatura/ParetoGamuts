function f1 = f1_mass_constrained_z_1(x,z)
%F1_MASS_CONSTRAINED_Z_1
%    F1 = F1_MASS_CONSTRAINED_Z_1(X1,X2)
%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-May-2020 23:59:20
% Includes scripted file modifications (LBracket.fixFileInputs) to match system specs

% unpack the block vars into the variable names used by the function
x1 = x(1, :).';
x2 = x(2, :).';

% Original content from Symbolic Math Toolbox
t2 = x1.*5.99999700000025e-1;
t3 = t2+3.999998000000167e-1;
t4 = x1.*5.999998999999756e-4;
t5 = t4+3.999999333333171e-4;
t6 = x1.*1.999999000000083e-1;
t7 = t6+1.666499333361167e-1;
t8 = 1.0./t3;
t9 = t4+1.003999499333375e-1;
t10 = 1.0./t5;
t31 = t3.*t9.*t10;
t11 = t2-t31+2.998998000166833e-1;
t12 = (t5.*t8.*t11)./6.0;
t13 = x2.*(2.0./5.0);
t14 = t13+1.0./1.0e1;
t15 = conj(x1);
t16 = t15.*5.99999700000025e-1;
t17 = t15.*5.999998999999756e-4;
t18 = t17+3.999999333333171e-4;
t19 = t16+3.999998000000167e-1;
t20 = t7.*t14;
t21 = x1.*3.999998000000167e-1;
t22 = t21+2.333498666639e-1;
t23 = t14.*t22;
t24 = t20+t23;
t25 = t15.*2.001999999666667e-1;
t26 = 1.0./t19;
t27 = 1.0./t18;
t28 = t17+1.003999499333375e-1;
t33 = t19.*t27.*t28;
t29 = t16-t33+2.998998000166833e-1;
t30 = x1.*2.001999999666667e-1;
t32 = -t12+t30+1.501333249777785e-1;
t34 = t15.*1.996000000666667e-1;
t35 = (t18.*t26.*t29)./6.0;
t36 = t6+6.664993333611667e-2;
t37 = x1.*3.998000000333333e-1;
t38 = t12+t37+2.498666750222215e-1;
t39 = conj(x2);
t40 = t39.*(2.0./5.0);
t41 = t40+1.0./1.0e1;
t42 = x1.*1.996000000666667e-1;
t43 = -t12+t42+1.497333250444451e-1;
t44 = t15.*1.999999000000083e-1;
t45 = t5.*t8.*t11.*(5.0./6.0);
t46 = t30+t45+1.501333249777785e-1;
t47 = -t12+t42+4.973337504444098e-2;
t48 = t21+2.332498666805667e-1;
t49 = t44+6.664993333611667e-2;
t50 = t14.*t32;
t51 = t14.*t43;
t52 = t36.*t38;
t54 = t14.*t36;
t53 = t20-t54;
t55 = t18.*t26.*t29.*(5.0./6.0);
t56 = t25+t55+1.501333249777785e-1;
t57 = t7.*t46;
t58 = t34-t35+1.497333250444451e-1;
t59 = t44+1.666499333361167e-1;
t60 = t14.*t46;
t61 = t14.*t38;
t62 = t14.*t48;
t63 = t54+t62;
t64 = t15.*3.999998000000167e-1;
t65 = t64+2.333498666639e-1;
t66 = t60-t14.*t47;
t67 = t50+t61;
t68 = t60+t61;
t69 = t7.*t38;
t70 = t15.*3.998000000333333e-1;
t71 = t35+t70+2.498666750222215e-1;
t72 = t64+2.332498666805667e-1;
t73 = t25-t35+1.501333249777785e-1;
f1 = t41.*(t52+t57).*(-1.0./6.0)+(t65.*(t50+t51))./6.0-(t24.*t56)./6.0+(t53.*t56)./3.0+(t53.*t58)./6.0+(t49.*t66)./6.0-(t49.*t67)./6.0-(t49.*t68)./6.0-(t56.*t63)./6.0+(t58.*t63)./3.0+(t53.*t71)./6.0+(t53.*t73)./6.0+(t59.*t67)./6.0+(t59.*t68)./6.0+(t65.*t66)./3.0-((t23-t14.*t48).*(t34-(t18.*t26.*t29)./6.0+4.973337504444098e-2))./6.0+(t24.*(t25-(t18.*t26.*t29)./6.0+1.501333249777785e-1))./6.0+(t58.*(t20-t23))./6.0+(t49.*(t50-t51))./6.0-(t41.*(t52-t69))./2.0+(t65.*(t50-t60))./6.0+(t41.*(t7.*t43+t22.*t32))./2.0-(t41.*(t22.*t47-t43.*t48))./6.0-(t41.*(t36.*t47+t46.*t48))./2.0+(t41.*(t69+t7.*t32))./3.0-(t41.*(t57-t32.*t36))./6.0-(t59.*(t50-t14.*t43))./6.0-(t41.*(t52+t36.*t46))./6.0+(t72.*(t51-t14.*t47))./6.0-(t14.*t22.*(t34-t35+4.973337504444098e-2))./3.0+(t7.*t14.*t73)./3.0-(t14.*t46.*t49)./3.0-(t14.*t47.*t49)./3.0+t14.*t38.*t59.*(2.0./3.0)-(t14.*t48.*t56)./3.0-t14.*t36.*t71.*(2.0./3.0)+(t14.*t43.*t72)./3.0;

