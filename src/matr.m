function S = matr(xi)
% make so(3) element from 3x1 vector or
% make se(3) element from 6x1 vector
% - written by Jin Seob Kim

N = size(xi,1);

if N == 1
    S = [0,-xi;xi,0];
    
elseif N == 3
    S = zeros(3);
    
    S(1,2) = -xi(3);
    S(1,3) = xi(2);
    S(2,1) = xi(3);
    S(2,3) = -xi(1);
    S(3,1) = -xi(2);
    S(3,2) = xi(1);
    
elseif N == 6
    S = zeros(4);

    w = xi(1:3);
    v = xi(4:6);

    S(1,2) = -w(3);
    S(1,3) = w(2);
    S(2,1) = w(3);
    S(2,3) = -w(1);
    S(3,1) = -w(2);
    S(3,2) = w(1);
    S(1,4) = v(1);
    S(2,4) = v(2);
    S(3,4) = v(3);
else
    disp('input is not of correct form')
end
   