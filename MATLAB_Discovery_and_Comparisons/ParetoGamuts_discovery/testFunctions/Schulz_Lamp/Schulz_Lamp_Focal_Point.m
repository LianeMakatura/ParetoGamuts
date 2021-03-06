function fHandle = Schulz_Lamp_Focal_Point(in1)
%SCHULZ_LAMP_FOCAL_POINT
%    FHANDLE =SCHULZ_LAMP_FOCAL_POINT(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    20-May-2020 20:39:24

%Automatically generated by Schulz_Lamp.createProblem()
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
x7 = in1(7,:);
x8 = in1(8,:);
x9 = in1(9,:);
x10 = in1(10,:);
x11 = in1(11,:);
x12 = in1(12,:);
x13 = in1(13,:);
x14 = in1(14,:);
x15 = in1(15,:);
x16 = in1(16,:);
x17 = in1(17,:);
x18 = in1(18,:);
x19 = in1(19,:);
x20 = in1(20,:);
x21 = in1(21,:);
t2 = x1.*2.0;
t3 = x2.*2.0;
t4 = x3.*2.0;
fHandle = sqrt(abs(t2+x7.*2.0+x16.*(3.0./5.0)-1.1e+1).^2+abs(t3+x8.*2.0+x17.*(3.0./5.0)-3.8e+1./5.0).^2+abs(t4+x9.*2.0+x18-2.0).^2)./6.0e+1+sqrt(abs(t2+x4.*2.0+x13.*(3.0./5.0)-3.8e+1./5.0).^2+abs(t3+x5.*2.0+x14.*(3.0./5.0)-3.8e+1./5.0).^2+abs(t4+x6.*2.0+x15-2.0).^2)./6.0e+1+sqrt(abs(t3+x11.*2.0+x20.*(3.0./5.0)-1.1e+1).^2+abs(t2+x10.*2.0+x19.*(3.0./5.0)-3.8e+1./5.0).^2+abs(t4+x12.*2.0+x21-2.0).^2)./6.0e+1;
