function f1 = f1_mass_constrained_z_0_1(x,z)
%F1_MASS_CONSTRAINED_Z_0_1
%    F1 = F1_MASS_CONSTRAINED_Z_0_1(X1,X2)
%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-May-2020 23:57:47
% Includes scripted file modifications (LBracket.fixFileInputs) to match system specs

% unpack the block vars into the variable names used by the function
x1 = x(1, :).';
x2 = x(2, :).';

% Original content from Symbolic Math Toolbox
t2 = x1.*4.565551683920277e-1;
t3 = x1.*3.893037094795364e-1;
t4 = t3+2.595358063196909e-1;
t5 = t2+3.043701122613518e-1;
t6 = x2.*(2.0./5.0);
t7 = t6+1.0./1.0e1;
t8 = 1.0./t5;
t9 = 1.0./t4;
t10 = t3+3.356283343850289e-1;
t18 = t5.*t9.*t10;
t11 = t2-t18+1.394861606814291e-1;
t12 = conj(x1);
t13 = t12.*4.565551683920277e-1;
t14 = t12.*3.893037094795364e-1;
t15 = t14+2.595358063196909e-1;
t16 = t13+3.043701122613518e-1;
t17 = x1.*5.953580631969091e-2;
t19 = (t4.*t8.*t11)./6.0;
t20 = t17+t19+1.031009776009089e-1;
t21 = t7.*t20;
t22 = x1.*3.297679031598455e-1;
t23 = t4.*t8.*t11.*(5.0./6.0);
t24 = t22+t23+2.3252735678412e-1;
t25 = t7.*t24;
t26 = t21+t25;
t27 = t12.*1.521850561306759e-1;
t28 = t27+2.397604549046348e-2;
t29 = x1.*2.702320968401545e-1;
t30 = t19+t29+1.6747264321588e-1;
t31 = t7.*t30;
t32 = -t19+t22+2.3252735678412e-1;
t33 = t7.*t32;
t34 = t31+t33;
t35 = t25+t31;
t36 = t27+1.239760454904635e-1;
t37 = t12.*3.043701122613518e-1;
t38 = t37+1.803940667708883e-1;
t39 = t17+t19+2.700844953557095e-2;
t40 = t7.*t39;
t41 = t33+t40;
t42 = x1.*3.043701122613518e-1;
t43 = x1.*1.521850561306759e-1;
t44 = t12.*3.297679031598455e-1;
t45 = 1.0./t16;
t46 = t14+3.356283343850289e-1;
t47 = 1.0./t15;
t54 = t16.*t46.*t47;
t48 = t13-t54+1.394861606814291e-1;
t49 = t42+1.803940667708883e-1;
t50 = t7.*t49;
t51 = t43+1.239760454904635e-1;
t52 = t7.*t51;
t53 = t50+t52;
t55 = t15.*t45.*t48.*(5.0./6.0);
t56 = t44+t55+2.3252735678412e-1;
t57 = t43+2.397604549046348e-2;
t58 = conj(x2);
t59 = t58.*(2.0./5.0);
t60 = t59+1.0./1.0e1;
t61 = t30.*t57;
t62 = t7.*t57;
t66 = (t15.*t45.*t48)./6.0;
t63 = t44-t66+2.3252735678412e-1;
t64 = t42+1.155101151909656e-1;
t65 = t7.*t64;
t67 = t62+t65;
t68 = t12.*5.953580631969091e-2;
t69 = t24.*t51;
t70 = t30.*t51;
t71 = t66+t68+2.700844953557095e-2;
t72 = t52-t62;
t73 = t66+t68+1.031009776009089e-1;
t74 = t12.*2.702320968401545e-1;
t75 = t66+t74+1.6747264321588e-1;
t76 = t37+1.155101151909656e-1;
f1 = t60.*(t61+t69).*(-1.0./6.0)+(t26.*t28)./6.0-(t28.*t34)./6.0-(t28.*t35)./6.0+(t26.*t38)./3.0+(t28.*t41)./6.0+(t34.*t36)./6.0+(t35.*t36)./6.0-(t36.*t41)./6.0-(t53.*t56)./6.0+(t53.*t63)./6.0-(t56.*t67)./6.0+(t56.*t72)./3.0+(t63.*t72)./6.0-(t67.*t71)./3.0-(t71.*t72)./6.0+(t72.*t75)./6.0-(t38.*(t25-t33))./6.0+(t76.*(t21-t40))./6.0+(t71.*(t50-t52))./6.0+(t73.*(t50-t65))./6.0-(t60.*(t61-t70))./2.0+(t60.*(t20.*t57-t24.*t64))./2.0+(t60.*(t32.*t49-t39.*t51))./2.0+(t60.*(t20.*t49-t39.*t64))./6.0+(t38.*(t33-t7.*t39))./6.0-(t60.*(t61+t24.*t57))./6.0+(t60.*(t70+t32.*t51))./3.0-(t60.*(t69-t32.*t57))./6.0+(t7.*t20.*t28)./3.0-(t7.*t24.*t28)./3.0+t7.*t30.*t36.*(2.0./3.0)+(t7.*t51.*t63)./3.0-(t7.*t39.*t76)./3.0-(t7.*t56.*t64)./3.0+(t7.*t49.*t73)./3.0-t7.*t57.*t75.*(2.0./3.0);

