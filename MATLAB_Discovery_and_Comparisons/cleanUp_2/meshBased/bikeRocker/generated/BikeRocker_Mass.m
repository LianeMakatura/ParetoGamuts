function fHandle = BikeRocker_Mass(in1,z1)
%BIKEROCKER_MASS
%    FHANDLE = BIKEROCKER_MASS(IN1,Z1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    26-Jan-2021 11:16:12

%Automatically generated by BikeRocker.createProblem()
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
t2 = x1.*2.0e1+2.5e1;
t6 = z1.*4.9e1;
t3 = t6-9.0;
t4 = x2.*7.0;
t5 = t4+7.0;
t7 = t6+4.1e1;
fHandle = ((x3.*5.0+5.0).*(sqrt(1.3e1).*t5.*(5.0e1./3.0)+t5.*sqrt(t3.^2+1.0e4)+t5.*sqrt(t7.^2+4.444444444444444e3)+t2.^2))./5.3111e4-2.192577808740186e-1;
