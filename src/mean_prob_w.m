%% mean_w_prob
% This function is to determine the mean omega based on the probability density
%
% - written by: Dimitri Lezcano


function w_mean = mean_prob_w(prob, w_i, w_j, w_k)
    % grid the omega values for averaging
    [w1_grid, w2_grid, w3_grid] = ndgrid(w_i, w_j, w_k);
    
    w_mean = zeros(3, size(prob,4));
    
    % iterate through each arclength
    for l = 1:size(prob, 4)
       prob_l = prob(:,:,:,l);
       
       % take the mean
       w1_mean = sum(w1_grid .* prob_l, 'all')/sum(prob_l, 'all');
       w2_mean = sum(w2_grid .* prob_l, 'all')/sum(prob_l, 'all');
       w3_mean = sum(w3_grid .* prob_l, 'all')/sum(prob_l, 'all');
       
       % add the mean to the data
       w_mean(:,l) = [w1_mean; w2_mean; w3_mean];
       
    end
end