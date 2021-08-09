function fHandle = Turbine_Mass(x,z)
%TURBINE_MASS
%    FHANDLE = TURBINE_MASS(X1,X2,X3)
%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    09-Aug-2021 18:15:10
% Includes scripted file modifications (Turbine.fixFileInputs) to match system specs

% unpack the block vars into the variable names used by the function
x1 = x(1, :).';
x2 = x(2, :).';
x3 = x(3, :).';
z1 = z(1, :).';

% Original content from Symbolic Math Toolbox
fHandle = ((x2.*9.0+6.0).*(x1.*6.0e1+4.0e1).*8.4678717818753e-3)./(pi.*(x3.*1.5e1+5.0))-1.078162836318169e-2;
