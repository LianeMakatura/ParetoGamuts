function hess = ZDT4_f1_hess(in1,z1)
%ZDT4_F1_HESS
%    HESS = ZDT4_F1_HESS(IN1,Z1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    27-Jan-2021 18:30:23

x1 = in1(1,:);
x2 = in1(2,:);
hess = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,3]);
