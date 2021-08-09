function f2 = f2_rcom_constrained_z_0_2(x,z)
%F2_RCOM_CONSTRAINED_Z_0_2
%    F2 = F2_RCOM_CONSTRAINED_Z_0_2(X1,X2)
%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-May-2020 23:57:57
% Includes scripted file modifications (LBracket.fixFileInputs) to match system specs

% unpack the block vars into the variable names used by the function
x1 = x(1, :).';
x2 = x(2, :).';

% Original content from Symbolic Math Toolbox
t2 = x2.*(2.0./5.0)-9.0./1.0e1;
t3 = x1.*1.618739039764216e-1-8.685390596211128e-1;
t4 = x1.*4.856217119292648e-1;
t5 = t4+3.237478079528432e-1;
t6 = x1.*3.523798417941783e-1;
t7 = t6+2.349198945294522e-1;
t8 = x1.*(-3.174599472647261e-1)+(t7.*(t4-(t5.*(t6+3.15856846517663e-1))./t7+1.650178343204802e-1))./(t5.*6.0)+7.748705431588141e-1;
f2 = t2.^2./3.0+t3.^2./3.0+t8.^2./3.0;
