function [total_er, A_er, tip_er] = TipAerror_ref_Dimitri(pmat_total, lengths, L_ref)
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
idx_ref = find(lengths == L_ref);
pmat_ref = pmat_total{idx_ref};


%% calculate the error between each of the insertions
A_error = 0;       % area error measurements
tip_error = 0;     % tip error measurements

num_data_pts = 0;  % for averaging

for j = [1:idx_ref-1, idx_ref+1:length(pmat_total)]
    % set-up
    num_data_pts = num_data_pts + 1; % increment the data counter
    Lj = lengths(j); % insertion lengths
    N_ref = L_ref/ds+1; Nj = Lj/ds+1; 
    s_ref = linspace(0, L_ref, N_ref); sj = linspace(0, Lj, Nj);

    pmatj = pmat_total{j};

    % error using area of a parameterized curve
    xi = pmat_ref(1,:); yi = pmat_ref(2,:); zi = pmat_ref(3,:);
    xj = pmatj(1,:); yj = pmatj(2,:); zj = pmatj(3,:);

    % indices for L's
    if L_ref < Lj
        ixLi = find(sj == L_ref);
        
    else
        ixLi = find(s_ref == Lj);
    
    end

    % interpolation ( for point correspondence )
    zij = sort(unique([zi(1:ixLi) zj(1:ixLi)]));

    xi_interp = interp1(zi, xi, zij, 'linear', 'extrap');
    xj_interp = interp1(zj, xj, zij, 'linear', 'extrap');

    yi_interp = interp1(zi, yi, zij, 'linear', 'extrap');
    yj_interp = interp1(zj, yj, zij, 'linear', 'extrap');

    % Area error
    dr_ij = vecnorm([xi_interp - xj_interp;
                     yi_interp - yj_interp]);
    dz = diff(zij);

    A_error = A_error + sum(dr_ij(2:end).*dz);

    % Tip Error
    tip_error = tip_error + L_ref*vecnorm([xi(end) - xj(ixLi);
                                              yi(end) - yj(ixLi); 
                                              zi(end) - zj(ixLi)]);

end
% average the errors
A_error = A_error/num_data_pts;
tip_error = tip_error/num_data_pts;


%% Perform the averaging
A_er = mean(A_error);
tip_er = mean(tip_error);

total_er = A_er + tip_er;


