function fHandle = ZDT2_f2(in1,z1)
%ZDT2_F2
%    FHANDLE = ZDT2_F2(IN1,Z1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    09-Aug-2021 17:27:15

%Automatically generated by ZDT2.createProblem()
x1 = in1(1,:);
x2 = in1(2,:);
t2 = x2.*(9.0./2.0);
t3 = z1.*(9.0./2.0);
t4 = t2+t3+1.0;
fHandle = t4.*(1.0./t4.^2.*x1.^2-1.0).*(-1.0./1.0e1);
