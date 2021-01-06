function v=vect(V)
% This function gives the "right" dual vector
% of given homogeneous transformation matrix.
% - made by Jin Seob Kim.

if size(V,1)==4
   v=[V(3,2);-V(3,1);V(2,1);V(1,4);V(2,4);V(3,4)];
elseif size(V,1)==3
   v=[V(3,2);-V(3,1);V(2,1)];
else 
   disp('Try again!!');
end



