function [li,phii,psii] = SPMinverseKinematicsCode(x,y,z,alpha,beta,gama)
%SPMINVERSEKINEMATICSCODE Summary of this function goes here
%   Detailed explanation goes here
li = zeros(6,1);
phii = zeros(6,1);
psii = zeros(6,1);

%Constants
rb = 0.175;
rt = 0.135;
gammab = 0.2985;
gammat = 0.6573;
% x = 0.1; y = 0.1; z = 0.7;
% alpha = 20*pi/180; beta = 15*pi/180; gama = 30*pi/180;

%Vectors and Matrices
O1=[0;0;0];
O2=[x;y;z];
p = O2-O1;
b = zeros(3,6);
t = zeros(3,6);
l = zeros(3,6);
l_mod = zeros(1,6);
d = zeros(3,6);
Rx = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
Ry = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
Rz = [cos(gama) -sin(gama) 0; sin(gama) cos(gama) 0; 0 0 1];
R = Rx*Ry*Rz; %First about x-axis, then about new y-axis, then about new z-axis

%Position coordinates
b(:,1)=[rb*cos(-gammab);rb*sin(-gammab);0];
b(:,2)=[rb*cos(gammab);rb*sin(gammab);0];
b(:,3)=[rb*cos((2*pi/3)-gammab);rb*sin((2*pi/3)-gammab);0];
b(:,4)=[rb*cos((2*pi/3)+gammab);rb*sin((2*pi/3)+gammab);0];
b(:,5)=[rb*cos((4*pi/3)-gammab);rb*sin((4*pi/3)-gammab);0];
b(:,6)=[rb*cos((4*pi/3)+gammab);rb*sin((4*pi/3)+gammab);0];

t(:,1)=[rt*cos(-gammat);rt*sin(-gammat);0];
t(:,2)=[rt*cos(gammat);rt*sin(gammat);0];
t(:,3)=[rt*cos((2*pi/3)-gammat);rt*sin((2*pi/3)-gammat);0];
t(:,4)=[rt*cos((2*pi/3)+gammat);rt*sin((2*pi/3)+gammat);0];
t(:,5)=[rt*cos((4*pi/3)-gammat);rt*sin((4*pi/3)-gammat);0];
t(:,6)=[rt*cos((4*pi/3)+gammat);rt*sin((4*pi/3)+gammat);0];

for i = 1:6
    d(:,i) = p + (R*t(:,i));
    l(:,i) = p + (R*t(:,i)) - b(:,i);
    l_mod(1,i) = sqrt(dot(l(:,i),l(:,i)));
end
li = transpose(l_mod);
for i = 1:6
    psii(i,1) = asin(-(l(2,i)/li(i,1)));
    phii(i,1) = asin(l(1,i)/(li(i,1)*cos(psii(i,1))));
end

end

