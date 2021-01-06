function [total_er, A_er, tip_er] = TipAerror_Dimitri(pmat_total, lengths)
%
% computes the error from area difference and tip position difference.
% only single bend insertion case.
% used for ideal needle insertion (i.e., kappa_0(s))
% pmat_total: cell array
% 
% - written by Dimitri Lezcano

Ncase = length(pmat_total);
ds = 0.5;

% sort the insertion scenarios so that they are increasing in length
[lengths, indices] = sort(lengths);
pmat_total = pmat_total(indices);

%% calculate the error between each of the insertions

A_error = zeros(1,length(lengths));       % area error measurements
tip_error = zeros(1,length(lengths));     % tip error measurements

for i = 1:length(lengths)-1
    num_data_pts = 0;  % for averaging

    for j = i+1:length(lengths)

        % set-up
        num_data_pts = num_data_pts + 1; % increment the data counter
        Li = lengths(i); Lj = lengths(j); % insertion lengths
        Ni = Li/ds+1; Nj = Lj/ds+1; 
        si = linspace(0, Li, Ni); sj = linspace(0, Lj, Nj);

        pmati = pmat_total{i}; pmatj = pmat_total{j};

        % error using area of a parameterized curve
        xi = pmati(1,:); yi = pmati(2,:); zi = pmati(3,:);
        xj = pmatj(1,:); yj = pmatj(2,:); zj = pmatj(3,:);

        % indices for L's
        ixLi = find(sj == Li);
        

        % interpolation ( for point correspondence )
        zij = sort(unique([zi zj(1:ixLi)]));

        xi_interp = interp1(zi, xi, zij, 'linear', 'extrap');
        xj_interp = interp1(zj, xj, zij, 'linear', 'extrap');

        yi_interp = interp1(zi, yi, zij, 'linear', 'extrap');
        yj_interp = interp1(zj, yj, zij, 'linear', 'extrap');

        % Area error
        dr_ij = vecnorm([xi_interp - xj_interp;
                         yi_interp - yj_interp]);
        dz = diff(zij);

        A_error(i) = A_error(i) + sum(dr_ij(2:end).*dz);

        % Tip Error
        tip_error(i) = tip_error(i) + Li*vecnorm([xi(end) - xj(ixLi);
                                                  yi(end) - yj(ixLi); 
                                                  zi(end) - zj(ixLi)]);

    end
    % average the errors
    A_error(i) = A_error(i)/num_data_pts;
    tip_error(i) = tip_error(i)/num_data_pts;

end

%% Perform the averaging
A_er = mean(A_error);
tip_er = mean(tip_error);

total_er = A_er + tip_er;


