function v=vect_se3(V)
% This function gives the "right" dual vector
% of given homogeneous transformation matrix.
% - made by Jin Seob Kim.

Vin = zeros(4);
Vin = real(V);
v=[Vin(3,2);-Vin(3,1);Vin(2,1);Vin(1,4);Vin(2,4);Vin(3,4)];
