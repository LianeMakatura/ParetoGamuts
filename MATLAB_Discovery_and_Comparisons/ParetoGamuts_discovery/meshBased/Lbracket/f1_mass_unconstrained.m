function f1 = f1_mass_unconstrained(x,z)
%F1_MASS_UNCONSTRAINED
%    F1 = F1_MASS_UNCONSTRAINED(X1,X2,Z1)
%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-May-2020 23:59:40
% Includes scripted file modifications (LBracket.fixFileInputs) to match system specs

% unpack the block vars into the variable names used by the function
x1 = x(1, :).';
x2 = x(2, :).';
z1 = z(1, :).';

% Original content from Symbolic Math Toolbox
t2 = conj(z1);
t3 = t2.*7.833981633974483e-1;
t4 = t3+7.863981633974483e-1;
t5 = sin(t4);
t6 = cos(t4);
t7 = x1.*(3.0./5.0);
t8 = t7+2.0./5.0;
t9 = x2.*(2.0./5.0);
t10 = t9+1.0./1.0e1;
t11 = x1./5.0;
t12 = 1.0./t5;
t13 = t6./1.0e1;
t14 = 1.0./t6;
t15 = t5./1.0e1;
t16 = t6.*t8;
t17 = t15+t16;
t18 = t5.*t14.*t17;
t23 = t5.*t8;
t19 = t13+t18-t23+1.0./1.0e1;
t20 = (t6.*t12.*t19)./6.0;
t21 = z1.*7.833981633974483e-1;
t22 = t21+7.863981633974483e-1;
t32 = t5./1.2e1;
t33 = t6.*t8.*(2.0./3.0);
t24 = t11+t20-t32-t33+2.0./1.5e1;
t25 = t5./6.0e1;
t26 = cos(t22);
t27 = sin(t22);
t28 = conj(x1);
t29 = t28.*(3.0./5.0);
t30 = t29+2.0./5.0;
t31 = t10.*t24;
t34 = (t6.*t8)./3.0;
t43 = t6.*t12.*t19.*(5.0./6.0);
t35 = t11+t25+t34-t43+2.0./1.5e1;
t44 = t10.*t35;
t36 = t31-t44;
t37 = t26./6.0e1;
t41 = (t27.*t30)./3.0;
t38 = t37-t41+1.0./1.5e1;
t39 = t11+t20+t25+t34+2.0./1.5e1;
t40 = t10.*t39;
t49 = x1.*(2.0./5.0);
t42 = t20+t25+t34-t49-4.0./1.5e1;
t45 = t11+t20+t25-t33+2.0./1.5e1;
t46 = t27.*t30.*(2.0./3.0);
t47 = t37+t46-1.0./3.0e1;
t48 = t10.*t45;
t52 = t10.*t42;
t50 = t40-t52;
t51 = -t37+t41+1.0./3.0e1;
t53 = t6./1.2e1;
t54 = t5.*t8.*(2.0./3.0);
t55 = t6./6.0e1;
t56 = t28./5.0;
t57 = 1.0./t27;
t58 = t26./1.0e1;
t59 = 1.0./t26;
t60 = t27./1.0e1;
t61 = t26.*t30;
t62 = t60+t61;
t63 = t27.*t59.*t62;
t72 = t27.*t30;
t64 = t58+t63-t72+1.0./1.0e1;
t65 = (t26.*t57.*t64)./6.0;
t66 = t53-t54+1.0./3.0e1;
t67 = t10.*t66;
t73 = (t5.*t8)./3.0;
t68 = t55-t73+1.0./1.5e1;
t69 = t10.*t68;
t70 = t67+t69;
t71 = t27./6.0e1;
t74 = (t26.*t30)./3.0;
t75 = -t55+t73+1.0./3.0e1;
t76 = t10.*t75;
t77 = t69+t76;
t83 = t26.*t30.*(2.0./3.0);
t78 = t56+t65+t71-t83+2.0./1.5e1;
t85 = t26.*t57.*t64.*(5.0./6.0);
t79 = t56+t71+t74-t85+2.0./1.5e1;
t80 = t54+t55-1.0./3.0e1;
t81 = t10.*t80;
t82 = t56+t65+t71+t74+2.0./1.5e1;
t84 = t76+t81;
t86 = conj(x2);
t87 = t86.*(2.0./5.0);
t88 = t87+1.0./1.0e1;
t89 = t42.*t68;
t90 = t35.*t75;
t91 = t42.*t75;
t92 = t26./1.2e1;
t93 = t28.*(-2.0./5.0)+t65+t71+t74-4.0./1.5e1;
f1 = (t47.*(t40+t48))./6.0-(t88.*(t89+t90))./6.0-(t88.*(t89+t91))./2.0+((t31-t10.*(t11+t20+t25-t6.*t8.*(2.0./3.0)+2.0./1.5e1)).*(t92-t27.*t30.*(2.0./3.0)+1.0./3.0e1))./6.0+(t36.*t38)./6.0-(t36.*t47)./3.0+(t38.*t50)./6.0+(t50.*t51)./6.0-(t70.*t78)./3.0+(t70.*t79)./6.0+(t77.*t78)./6.0+(t77.*t79)./3.0+(t77.*t82)./6.0-(t79.*t84)./6.0+(t82.*t84)./6.0-(t77.*t93)./6.0+(t47.*(t40-t44))./6.0-(t51.*(t40-t48))./6.0+(t51.*(t44-t52))./6.0+(t78.*(t76-t81))./6.0-((t81+t10.*(t53-t5.*t8.*(2.0./3.0)+1.0./3.0e1)).*(t27.*(-1.0./1.2e1)+t56+t65-t26.*t30.*(2.0./3.0)+2.0./1.5e1))./6.0+(t88.*(t24.*t68+t35.*t66))./2.0+(t88.*(t35.*t68-t42.*t68))./6.0-(t88.*(t24.*t80+t45.*t66))./6.0+(t88.*(t39.*t80+t45.*t75))./2.0-(t38.*(t40-t10.*t45))./6.0+(t38.*(t44-t10.*t42))./6.0-(t88.*(t90+t39.*t68))./6.0-(t88.*(t91-t39.*t75))./3.0-(t10.*t80.*(t27.*(-1.0./1.2e1)+t56+t65-t83+2.0./1.5e1))./3.0+(t10.*t24.*t38)./3.0+(t10.*t35.*t38)./3.0-t10.*t42.*t51.*(2.0./3.0)+(t10.*t66.*t79)./3.0+(t10.*t75.*t82)./3.0-t10.*t68.*t93.*(2.0./3.0)-(t10.*t45.*(-t46+t92+1.0./3.0e1))./3.0;
