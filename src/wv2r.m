function positions = wv2r(wv, L)
% Function to return a position values for each of the deformed values
%
% wv   - The angular deformation vector data
% L    - the total length of insertoin
%
% - written by: Dimitri Lezcano
% - credit to: Jin Seob Kim


    e3 = [0; 0; 1];
    positions = zeros(3,length(wv));
    ds = L / (length(positions) - 1);
    
%% Instantiations - Taken from main_singlebend_v2_tissue_90total.m
    Rmat = zeros(3,3,length(positions));
    Rmat(:,:,1) = eye(3);
    
%% Calculations - Taken from main_singlebend_v2_tissue_90total.m
    for i = 2:length(positions)   
    % orientation
        W = matr(1/2*(wv(:,i-1) + wv(:,i)));
        Rmat(:,:,i) = Rmat(:,:,i-1)*expm(ds*W);

        % position
        e3vec = squeeze(Rmat(:,3,1:i));
        if i == 2
            positions(:,i) = positions(:,i-1) + squeeze(Rmat(:,3,i))*ds;
        else
            positions(:,i) = Simpson_vec_int(e3vec,ds);
        end        
    end
end