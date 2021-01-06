function Y=Adjoint(H)
% To find the Adjoint matrix for H.
% H=[R,b;zeros(1,3),1] 
% where R belongs to SO(3), b to 3x1 vector
% The corresponding se(3) is [w;v].
% - made by Jin Seob Kim.

R=H(1:3,1:3);
b=H(1:3,4);
B=[0,-b(3),b(2);b(3),0,-b(1);-b(2),b(1),0];

Y=[R,zeros(3);B*R,R];

