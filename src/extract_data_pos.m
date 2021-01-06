function [RMSE,y_th,z_th] = extract_data_pos(y_data,z_data,pmat)
%
% y_data, z_data: position from image analysis
% pmat: position from the model
%
% - written by Jin Seob (Jesse) Kim

%======================================%
%%%%%%% extract theoretical data %%%%%%%
%======================================%
M = length(z_data); % data point number
N = size(pmat,2); % model point number

z = pmat(3,:)';
y = pmat(2,:)';

z_th = zeros(M,1);
y_th = zeros(M,1);
for k = 1:M
    for j = 1:N
        if abs(z(j) - z_data(k)) < 1e-7
            z_th(k) = z(j);
            y_th(k) = y(j);
        elseif (j+1 <= N) && (z(j) < z_data(k) && z(j+1) > z_data(k))
            dz1 = z_data(k) - z(j);
            dz2 = z(j+1) - z_data(k);
            dz = z(j+1) - z(j);
            z_th(k) = dz2/dz*z(j) + dz1/dz*z(j+1);
            y_th(k) = dz2/dz*y(j) + dz1/dz*y(j+1);            
        end
    end
end    

%========================================%
%%%%%%% Mean Square Error function %%%%%%%
%========================================%
RMSE = sqrt(1/M*sum((y_data - y_th).^2));
