function f1 = f1_mass_constrained_z_0_3(x,z)
%F1_MASS_CONSTRAINED_Z_0_3
%    F1 = F1_MASS_CONSTRAINED_Z_0_3(X1,X2)
%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-May-2020 23:58:09
% Includes scripted file modifications (LBracket.fixFileInputs) to match system specs

% unpack the block vars into the variable names used by the function
x1 = x(1, :).';
x2 = x(2, :).';

% Original content from Symbolic Math Toolbox
t2 = x1.*5.117094573379158e-1;
t3 = x1.*3.132944801156505e-1;
t4 = t3+2.08862986743767e-1;
t5 = t2+3.411396382252772e-1;
t6 = 1.0./t5;
t7 = 1.0./t4;
t8 = t3+2.941478963000863e-1;
t11 = t5.*t7.*t8;
t9 = t2-t11+1.889238915393354e-1;
t10 = x1.*1.705698191126386e-1;
t12 = (t4.*t6.*t9)./6.0;
t13 = t10+3.834392162743544e-2;
t14 = x1.*3.044314933718835e-1;
t15 = conj(x2);
t16 = t15.*(2.0./5.0);
t17 = t16+1.0./1.0e1;
t18 = x1.*2.955685066281165e-1;
t19 = t12+t18+1.828315194926911e-1;
t20 = x2.*(2.0./5.0);
t21 = t20+1.0./1.0e1;
t22 = conj(x1);
t23 = t22.*1.705698191126386e-1;
t24 = t23+1.383439216274354e-1;
t25 = t19.*t21;
t26 = t4.*t6.*t9.*(5.0./6.0);
t27 = t14+t26+2.171684805073089e-1;
t28 = x1.*3.411396382252772e-1;
t29 = t22.*5.117094573379158e-1;
t30 = t22.*3.132944801156505e-1;
t31 = t30+2.08862986743767e-1;
t32 = t29+3.411396382252772e-1;
t33 = -t12+t14+2.171684805073089e-1;
t34 = t21.*t33;
t35 = x1.*8.862986743767021e-3;
t36 = t12+t35-8.305493763541869e-3;
t37 = t13.*t19;
t38 = t10+1.383439216274354e-1;
t39 = t13.*t21;
t40 = 1.0./t32;
t41 = 1.0./t31;
t42 = t30+2.941478963000863e-1;
t52 = t32.*t41.*t42;
t43 = t29-t52+1.889238915393354e-1;
t44 = (t31.*t40.*t43)./6.0;
t45 = t21.*t36;
t46 = t21.*t27;
t47 = t12+t35+7.697941579277743e-2;
t48 = t23+3.834392162743544e-2;
t49 = t22.*3.411396382252772e-1;
t50 = t28+1.505799699119e-1;
t53 = t21.*t38;
t51 = t39-t53;
t54 = t22.*3.044314933718835e-1;
t55 = t34+t45;
t56 = t27.*t38;
t57 = t49+2.027957165978418e-1;
t58 = t22.*8.862986743767021e-3;
t59 = t44+t58-8.305493763541869e-3;
t60 = t28+2.027957165978418e-1;
t61 = t21.*t50;
t62 = t21.*t47;
t63 = t46+t62;
t64 = -t44+t54+2.171684805073089e-1;
t65 = t21.*t60;
t66 = t53+t65;
t67 = t31.*t40.*t43.*(5.0./6.0);
t68 = t54+t67+2.171684805073089e-1;
t69 = t25+t34;
t70 = t25+t46;
t71 = t39+t61;
t72 = t49+1.505799699119e-1;
t73 = t44+t58+7.697941579277743e-2;
t74 = t22.*2.955685066281165e-1;
t75 = t44+t74+1.828315194926911e-1;
f1 = t17.*(t37+t56).*(-1.0./6.0)-(t24.*t55)./6.0+(t24.*t69)./6.0+(t24.*t70)./6.0+(t48.*t55)./6.0+(t51.*t59)./6.0+(t48.*t63)./6.0-(t51.*t64)./6.0-(t48.*t69)./6.0-(t48.*t70)./6.0-(t51.*t68)./3.0+(t57.*t63)./3.0-(t51.*t75)./6.0-(t59.*t71)./3.0+(t64.*t66)./6.0-(t66.*t68)./6.0-(t68.*t71)./6.0+(t57.*(t34-t45))./6.0+(t57.*(t34-t46))./6.0-(t59.*(t53-t65))./6.0+(t17.*(t60.*(t14-(t4.*t6.*t9)./6.0+2.171684805073089e-1)-t36.*t38))./2.0+(t17.*(t19.*t38+t33.*t38))./3.0+(t17.*(t13.*t47-t27.*t50))./2.0-(t17.*(t36.*t50-t47.*t60))./6.0-(t17.*(t37+t13.*t27))./6.0-(t17.*(t37-t19.*t38))./2.0-(t17.*(t56-t13.*t33))./6.0-(t72.*(t45-t21.*t47))./6.0-(t73.*(t61-t21.*t60))./6.0+t19.*t21.*t24.*(2.0./3.0)-(t21.*t27.*t48)./3.0-t13.*t21.*t75.*(2.0./3.0)+(t21.*t47.*t48)./3.0+(t21.*t38.*t64)./3.0-(t21.*t36.*t72)./3.0-(t21.*t50.*t68)./3.0+(t21.*t60.*t73)./3.0;
