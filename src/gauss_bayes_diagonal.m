%% gauss_bayes_diagonal.m
%
% function to handle Gaussian Bayesian fusion for diagonal covariance
%
% - written by: Dimitri Lezcano

function [mu3, cov3] = gauss_bayes_diagonal(mu1, cov1, mu2, cov2)

    %% arguments block
    arguments
        mu1 (:,1);
        cov1 (:,1) {mustBeEqualSize(mu1, cov1)};
        mu2 (:,1) {mustBeEqualSize(mu2, mu1)};
        cov2 (:,1) {mustBeEqualSize(cov2, cov1)};
    end
    %% calculation
    % direct calculation
    cov3 = 1./(1./cov1 + 1./cov2);
    mu3 = cov3 .* (mu1./cov1 + mu2./cov2);
    
    % alternate method
    % [mu3, Cov3] = gauss_bayes(mu1, diag(cov1), mu2, diag(cov2));
    % cov3 = diag(Cov3);
    
end