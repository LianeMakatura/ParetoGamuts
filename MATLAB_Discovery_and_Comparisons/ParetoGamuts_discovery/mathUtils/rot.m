function R= rot(V,theta)
%This function returns the 3D rotation matrix about an arbitary vector V
%passing the origin.
%
%This may serves as an extension to the matlab built in funciton :
%rotx, roty and rotz.
%This can also be extended for V not going through origin by simple
%translation
%
%Input:
%V: 3 by 1 or 1 by 3 arbitrary rotation vector, it does not need to be
%normalized
%
%theta:The rotation angle is positive if the rotation is in the counter-clockwise 
%direction when viewed by an observer looking along the V towards the origin.
%
%Output:
%R: 3 by 3 rotation matrix
%To use R on a column vector v: R*v
%To use R on a row vector w: w*R'
%
% Glenn Murray's step by step approach
%
% (1) Translate space so that the rotation axis passes through the origin.
% (2) Rotate space about the z axis so that the rotation axis lies in the xz plane.
% (3) Rotate space about the y axis so that the rotation axis lies along the z axis.
% (4) Perform the desired rotation by ? about the z axis.
% (5) Apply the inverse of step (3).
% (6) Apply the inverse of step (2).
% (7) Apply the inverse of step (1).
%
% Glenn Murray's online article :
%http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
%
%Example:
%Rotate pi/3 about x-axis
% rot([1,0,0],pi/3)
% =
% 
%     1.0000         0         0
%          0    0.5000   -0.8660
%          0    0.8660    0.5000
%
%By Tianxiao Jiang
%Unversity of Houston
%jtxinnocence@gmail.com

V=V/norm(V);

R=[V(1)^2+(1-V(1)^2)*cos(theta),             V(1)*V(2)*(1-cos(theta))-V(3)*sin(theta), V(1)*V(3)*(1-cos(theta))+V(2)*sin(theta);
   V(1)*V(2)*(1-cos(theta))+V(3)*sin(theta), V(2)^2+(1-V(2)^2)*cos(theta),             V(2)*V(3)*(1-cos(theta))-V(1)*sin(theta);
   V(1)*V(3)*(1-cos(theta))-V(2)*sin(theta), V(2)*V(3)*(1-cos(theta))+V(1)*sin(theta), V(3)^2+(1-V(3)^2)*cos(theta)];

end
