function f2 = f2_rcom_constrained_z_0_5(x,z)
%F2_RCOM_CONSTRAINED_Z_0_5
%    F2 = F2_RCOM_CONSTRAINED_Z_0_5(X1,X2)
%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-May-2020 23:58:29
% Includes scripted file modifications (LBracket.fixFileInputs) to match system specs

% unpack the block vars into the variable names used by the function
x1 = x(1, :).';
x2 = x(2, :).';

% Original content from Symbolic Math Toolbox
t2 = x1.*5.54327719506772e-1;
t3 = x1.*2.296100594190539e-1;
t4 = t2+3.695518130045147e-1;
t5 = t3+1.530733729460359e-1;
t6 = x1.*(-2.76536686473018e-1)+(t5.*(t2-(t4.*(t3+2.454613261971646e-1))./t5+2.312834697680057e-1))./(t4.*6.0)+8.002442168094666e-1;
t7 = x2.*(2.0./5.0)-9.0./1.0e1;
t8 = x1.*1.847759065022573e-1-8.498607862045799e-1;
f2 = t6.^2./3.0+t7.^2./3.0+t8.^2./3.0;

