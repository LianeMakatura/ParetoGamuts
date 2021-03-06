function fHandle = ZDT6_f1(in1,z1)
%ZDT6_F1
%    FHANDLE = ZDT6_F1(IN1,Z1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    09-Aug-2021 17:31:01

%Automatically generated by ZDT6.createProblem()
x1 = in1(1,:);
t2 = sin(x1.*pi.*6.0);
t3 = t2.^2;
t4 = t3.^2;
fHandle = -t3.*t4.*exp(x1.*-4.0)+1.0;
