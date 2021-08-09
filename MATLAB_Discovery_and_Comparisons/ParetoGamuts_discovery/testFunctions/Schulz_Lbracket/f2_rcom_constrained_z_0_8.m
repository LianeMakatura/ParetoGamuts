function f2 = f2_rcom_constrained_z_0_8(x,z)
%F2_RCOM_CONSTRAINED_Z_0_8
%    F2 = F2_RCOM_CONSTRAINED_Z_0_8(X1,X2)
%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    10-May-2020 23:59:00
% Includes scripted file modifications (LBracket.fixFileInputs) to match system specs

% unpack the block vars into the variable names used by the function
x1 = x(1, :).';
x2 = x(2, :).';

% Original content from Symbolic Math Toolbox
t2 = x2.*(2.0./5.0)-9.0./1.0e1;
t3 = x1.*5.925565812827096e-1;
t4 = x1.*9.421622991049698e-2;
t5 = t4+6.281081994033132e-2;
t6 = t3+3.950377208551397e-1;
t7 = x1.*(-2.314054099701657e-1)+(t5.*(t3-(t6.*(t4+1.615702501541162e-1))./t5+2.793350158700569e-1))./(t6.*6.0)+8.292698216509254e-1;
t8 = x1.*1.975188604275699e-1-8.376045438791339e-1;
f2 = t2.^2./3.0+t7.^2./3.0+t8.^2./3.0;
