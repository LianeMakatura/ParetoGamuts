function fHandle = BikeStay_Mass(in1,z1)
%BIKESTAY_MASS
%    FHANDLE = BIKESTAY_MASS(IN1,Z1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    27-Jan-2021 19:43:24

%Automatically generated by BikeStay.createProblem()
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
t4 = x1.*3.0e1;
t2 = t4+5.0e1;
t3 = x1.*-6.0e1+z1.*9.8e1+2.32e2;
t5 = t4-1.3e2;
t6 = x3.*1.0e1;
t7 = t6+6.0;
fHandle = (t2.^2.*t7)./1.5725e5+(t7.*(x2.*1.0e1+1.0e1).*sqrt(t3.^2./4.0+t5.^2./4.0))./7.8625e4-9.4e1./6.29e2;
