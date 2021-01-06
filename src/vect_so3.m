function v=vect_so3(V)
% This function gives the "right" dual vector
% of given rotation matrix.
% - made by Jin Seob Kim.
Vin = zeros(3);
Vin = real(V);
v=[Vin(3,2);-Vin(3,1);Vin(2,1)];



