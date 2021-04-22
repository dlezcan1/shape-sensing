%% fit_ellipse2d.m
%
% function to fit a 2-D ellipse in 3-D space
%
% - written by: Dimitri Lezcano

function [center, Sigma, R, Sig_2d] = fit_ellipse2d(pts)
    %% Arguments Block
    arguments
        pts (3,:);
    end
    
    %% Compute the center point
    center = mean(pts, 2);
    
    %% Compute the ellipse parameters
    % zero-center the points
    dpts = pts - center;
    
    % compute ellipse normal
    [u, ~, ~] = svd(dpts);
    n = u(:,end);
    n = n / norm(n); % normalize
    
    % rotate the points to align with XY plane
    [~, max_idx] = max(vecnorm(dpts,2,1));
    r1 = dpts(:,max_idx);
    r1 = r1/norm(r1);
    R = [r1, cross(n, r1), n];
    
    dpts_R = R' * dpts;
    
    % perform 2D ellipse fitting
    A = [dpts_R(1,:).^2; dpts_R(1,:).*dpts_R(2,:); dpts_R(2,:).^2]';
    b = ones(size(A,1), 1);
    sig_2d_v = A\b;
    Sig_2d = inv([sig_2d_v(1), sig_2d_v(2);
              sig_2d_v(2), sig_2d_v(3)]);
    Sig_3d = [Sig_2d, zeros(2,1); zeros(1,3)];
  
    % rerotate the center back
    Sigma = R' * Sig_3d * R;
    
end
