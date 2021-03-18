%% gauss_bayes.m
%
% this is a function to provide a multivariate approximation of Bayesian fusion 
% of two gaussian distributions
%
% - written by: Dimitri Lezcano

function [mu3, cov3] = gauss_bayes(mu1, cov1, mu2, cov2)
    %% arguments block
    arguments
        mu1 (:,1);
        cov1 (:,:) {mustBeSquare, mustBeEqualSize(cov1, mu1, 1)};
        mu2 (:,1) {mustBeEqualSize(mu2, mu1)};
        cov2 (:,:) {mustBeSquare, mustBeEqualSize(cov1, cov2)};
    end
    %% Calculation
    inv_cov12 = inv(cov1 + cov2);
    cov3 = cov1 * inv_cov12 * cov2;
    mu3 = cov2*inv_cov12*mu1 + cov1*inv_cov12*mu2;

end