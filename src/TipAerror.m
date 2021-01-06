function [total_er, A_er, tip_error] = TipAerror(pmat_total)
%
% computes the error from area difference and tip position difference.
% only single bend insertion case.
% used for ideal needle insertion (i.e., kappa_0(s))
% pmat_total: cell array
% 
% - written by Jin Seob Kim

Ncase = length(pmat_total);
ds = 0.5;

%% =====================o===================================================
%
%                            5 insertion lengths
%
%==========================================================================
if Ncase == 5
    
L1 = 90; N1 = L1/ds+1; s1 = linspace(0,L1,N1);
L2 = 105; N2 = L2/ds+1; s2 = linspace(0,L2,N2);
L3 = 120; N3 = L3/ds+1; s3 = linspace(0,L3,N3);
L4 = 150; N4 = L4/ds+1; s4 = linspace(0,L4,N4);
L5 = 180; N5 = L5/ds+1; s5 = linspace(0,L5,N5);

pmat1 = pmat_total{1};
pmat2 = pmat_total{2};
pmat3 = pmat_total{3};
pmat4 = pmat_total{4};
pmat5 = pmat_total{5};

%% error using area of a parameterized curve
x1 = pmat1(3,:); y1 = pmat1(2,:);
x2 = pmat2(3,:); y2 = pmat2(2,:);
x3 = pmat3(3,:); y3 = pmat3(2,:);
x4 = pmat4(3,:); y4 = pmat4(2,:);
x5 = pmat5(3,:); y5 = pmat5(2,:);

% indices for L's
ixL1 = find(s2 == L1);
ixL2 = find(s3 == L2);
ixL3 = find(s4 == L3);
ixL4 = find(s5 == L4);

%>>> curve area error calculation
A_error = zeros(1,4);

% L = L1 case
A_error(1) = abs(Acurve_param(-y2(1:ixL1),x2(1:ixL1)) - Acurve_param(-y1,x1)) ...
    + abs(Acurve_param(-y3(1:ixL1),x3(1:ixL1)) - Acurve_param(-y1,x1)) ...
    + abs(Acurve_param(-y4(1:ixL1),x4(1:ixL1)) - Acurve_param(-y1,x1)) ...
    + abs(Acurve_param(-y5(1:ixL1),x5(1:ixL1)) - Acurve_param(-y1,x1)) ...
    + abs(Acurve_param(-y3(1:ixL1),x3(1:ixL1)) - Acurve_param(-y2(1:ixL1),x2(1:ixL1))) ...
    + abs(Acurve_param(-y4(1:ixL1),x4(1:ixL1)) - Acurve_param(-y2(1:ixL1),x2(1:ixL1))) ...
    + abs(Acurve_param(-y5(1:ixL1),x5(1:ixL1)) - Acurve_param(-y2(1:ixL1),x2(1:ixL1))) ...
    + abs(Acurve_param(-y4(1:ixL1),x4(1:ixL1)) - Acurve_param(-y3(1:ixL1),x3(1:ixL1))) ...
    + abs(Acurve_param(-y5(1:ixL1),x5(1:ixL1)) - Acurve_param(-y3(1:ixL1),x3(1:ixL1))) ...
    + abs(Acurve_param(-y5(1:ixL1),x5(1:ixL1)) - Acurve_param(-y4(1:ixL1),x4(1:ixL1)));
A_error(1) = A_error(1)/10;

% L = L2 case
A_error(2) = abs(Acurve_param(-y3(1:ixL2),x3(1:ixL2)) - Acurve_param(-y2,x2)) ...
    + abs(Acurve_param(-y4(1:ixL2),x4(1:ixL2)) - Acurve_param(-y2,x2)) ...
    + abs(Acurve_param(-y5(1:ixL2),x5(1:ixL2)) - Acurve_param(-y2,x2)) ...
    + abs(Acurve_param(-y4(1:ixL2),x4(1:ixL2)) - Acurve_param(-y3(1:ixL2),x3(1:ixL2))) ...
    + abs(Acurve_param(-y5(1:ixL2),x5(1:ixL2)) - Acurve_param(-y3(1:ixL2),x3(1:ixL2))) ...
    + abs(Acurve_param(-y5(1:ixL2),x5(1:ixL2)) - Acurve_param(-y4(1:ixL2),x4(1:ixL2)));
A_error(2) = A_error(2)/6;

% L = L3 case
A_error(3) = abs(Acurve_param(-y4(1:ixL3),x4(1:ixL3)) - Acurve_param(-y3,x3)) ...
    + abs(Acurve_param(-y5(1:ixL3),x5(1:ixL3)) - Acurve_param(-y3,x3)) ...
    + abs(Acurve_param(-y5(1:ixL3),x5(1:ixL3)) - Acurve_param(-y4(1:ixL3),x4(1:ixL3)));
A_error(3) = A_error(3)/3;

% L = L4 case
A_error(4) = abs(Acurve_param(-y5(1:ixL4),x5(1:ixL4)) - Acurve_param(-y4,x4));

% curve area error
A_er = sum(A_error)/4;

%>>> tip position error calculation
tip_error = zeros(1,4);

% L = L1 case
tip_error(1) = sqrt((x1(end) - x2(ixL1))^2 + (y1(end) - y2(ixL1))^2) ...
    + sqrt((x1(end) - x3(ixL1))^2 + (y1(end) - y3(ixL1))^2) ...
    + sqrt((x1(end) - x4(ixL1))^2 + (y1(end) - y4(ixL1))^2) ...
    + sqrt((x1(end) - x5(ixL1))^2 + (y1(end) - y5(ixL1))^2) ...
    + sqrt((x2(ixL1) - x3(ixL1))^2 + (y2(ixL1) - y3(ixL1))^2) ...
    + sqrt((x2(ixL1) - x4(ixL1))^2 + (y2(ixL1) - y4(ixL1))^2) ...
    + sqrt((x2(ixL1) - x5(ixL1))^2 + (y2(ixL1) - y5(ixL1))^2) ...
    + sqrt((x3(ixL1) - x4(ixL1))^2 + (y3(ixL1) - y4(ixL1))^2) ...
    + sqrt((x3(ixL1) - x5(ixL1))^2 + (y3(ixL1) - y5(ixL1))^2) ...
    + sqrt((x4(ixL1) - x5(ixL1))^2 + (y4(ixL1) - y5(ixL1))^2);
tip_error(1) = tip_error(1)/10*L1;

% L = L2 case
tip_error(2) = sqrt((x2(end) - x3(ixL2))^2 + (y2(end) - y3(ixL2))^2) ...
    + sqrt((x2(end) - x4(ixL2))^2 + (y2(end) - y4(ixL2))^2) ...
    + sqrt((x2(end) - x5(ixL2))^2 + (y2(end) - y5(ixL2))^2) ...
    + sqrt((x3(ixL2) - x4(ixL2))^2 + (y3(ixL2) - y4(ixL2))^2) ...
    + sqrt((x3(ixL2) - x5(ixL2))^2 + (y3(ixL2) - y5(ixL2))^2) ...
    + sqrt((x4(ixL2) - x5(ixL2))^2 + (y4(ixL2) - y5(ixL2))^2);
tip_error(2) = tip_error(2)/6*L2;

% L = L3 case
tip_error(3) = sqrt((x3(end) - x4(ixL3))^2 + (y3(end) - y4(ixL3))^2) ...
    + sqrt((x3(end) - x5(ixL3))^2 + (y3(end) - y5(ixL3))^2) ...
    + sqrt((x4(ixL3) - x5(ixL3))^2 + (y4(ixL3) - y5(ixL3))^2);
tip_error(3) = tip_error(3)/3*L3;

% L = L4 case
tip_error(4) = sqrt((x4(end) - x5(ixL4))^2 + (y4(end) - y5(ixL4))^2);
tip_error(4) = tip_error(4)*L4;

% total error
tip_er = sum(tip_error)/4;

% output
total_er = A_er + tip_er;

%% ========================================================================
%
%                            4 insertion lengths
%
%==========================================================================
elseif Ncase == 4
    
L1 = 60; N1 = L1/ds+1; s1 = linspace(0,L1,N1);
L2 = 90; N2 = L2/ds+1; s2 = linspace(0,L2,N2);
L3 = 120; N3 = L3/ds+1; s3 = linspace(0,L3,N3);
L4 = 150; N4 = L4/ds+1; s4 = linspace(0,L4,N4);

pmat1 = pmat_total{1};
pmat2 = pmat_total{2};
pmat3 = pmat_total{3};
pmat4 = pmat_total{4};

%% error using area of a parameterized curve
x1 = pmat1(3,:); y1 = pmat1(2,:);
x2 = pmat2(3,:); y2 = pmat2(2,:);
x3 = pmat3(3,:); y3 = pmat3(2,:);
x4 = pmat4(3,:); y4 = pmat4(2,:);

% indices for L's
ixL1 = find(s2 == L1);
ixL2 = find(s3 == L2);
ixL3 = find(s4 == L3);

%>>> curve area error calculation
A_error = zeros(1,3);

% L = 60 case
A_error(1) = abs(Acurve_param(-y2(1:ixL1),x2(1:ixL1)) - Acurve_param(-y1,x1)) ...
    + abs(Acurve_param(-y3(1:ixL1),x3(1:ixL1)) - Acurve_param(-y1,x1)) ...
    + abs(Acurve_param(-y4(1:ixL1),x4(1:ixL1)) - Acurve_param(-y1,x1)) ...
    + abs(Acurve_param(-y3(1:ixL1),x3(1:ixL1)) - Acurve_param(-y2(1:ixL1),x2(1:ixL1))) ...
    + abs(Acurve_param(-y4(1:ixL1),x4(1:ixL1)) - Acurve_param(-y2(1:ixL1),x2(1:ixL1))) ...
    + abs(Acurve_param(-y4(1:ixL1),x4(1:ixL1)) - Acurve_param(-y3(1:ixL1),x3(1:ixL1)));
A_error(1) = A_error(1)/6;

% L = 90 case
A_error(2) = abs(Acurve_param(-y3(1:ixL2),x3(1:ixL2)) - Acurve_param(-y2,x2)) ...
    + abs(Acurve_param(-y4(1:ixL2),x4(1:ixL2)) - Acurve_param(-y2,x2)) ...
    + abs(Acurve_param(-y4(1:ixL2),x4(1:ixL2)) - Acurve_param(-y3(1:ixL2),x3(1:ixL2)));
A_error(2) = A_error(2)/3;

% L = 120 case
A_error(3) = abs(Acurve_param(-y4(1:ixL3),x4(1:ixL3)) - Acurve_param(-y3,x3));

% curve area error
A_er = sum(A_error)/3;

%>>> tip position error calculation
tip_error = zeros(1,3);

% L = 60 case
tip_error(1) = sqrt((x1(end) - x2(ixL1))^2 + (y1(end) - y2(ixL1))^2) ...
    + sqrt((x1(end) - x3(ixL1))^2 + (y1(end) - y3(ixL1))^2) ...
    + sqrt((x1(end) - x4(ixL1))^2 + (y1(end) - y4(ixL1))^2) ...
    + sqrt((x2(ixL1) - x3(ixL1))^2 + (y2(ixL1) - y3(ixL1))^2) ...
    + sqrt((x2(ixL1) - x4(ixL1))^2 + (y2(ixL1) - y4(ixL1))^2) ...
    + sqrt((x3(ixL1) - x4(ixL1))^2 + (y3(ixL1) - y4(ixL1))^2);
tip_error(1) = tip_error(1)/6*L1;

% L = 90 case
tip_error(2) = sqrt((x2(end) - x3(ixL2))^2 + (y2(end) - y3(ixL2))^2) ...
    + sqrt((x2(end) - x4(ixL2))^2 + (y2(end) - y4(ixL2))^2) ...
    + sqrt((x3(ixL2) - x4(ixL2))^2 + (y3(ixL2) - y4(ixL2))^2);
tip_error(2) = tip_error(2)/3*L2;

% L = 120 case
tip_error(3) = sqrt((x3(end) - x4(ixL3))^2 + (y3(end) - y4(ixL3))^2);
tip_error(3) = tip_error(3)*L3;

% total error
tip_er = sum(tip_error)/3;

% output
total_er = A_er + tip_er;

%%=========================================================================
%
%             only these two cases!
%
%==========================================================================
else
    disp('Out of range for the cases')
end
