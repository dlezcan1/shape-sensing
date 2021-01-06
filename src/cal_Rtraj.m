function [Rmat,pmat] = cal_Rtraj(wv,ds,theta0)
%
% Given an omega vector array, compute the trajectoyr R(s)
% wv: 3 x N array
%
% - written by Jin Seob (Jesse) Kim

N = size(wv,2);

Rmat = zeros(3,3,N);
Rmat(:,:,1) = Rot_x(theta0);
pmat = zeros(3,N);
for i = 2:N
    % orientation
    W = matr(1/2*(wv(:,i-1) + wv(:,i)));
    Rmat(:,:,i) = Rmat(:,:,i-1)*expm(ds*W);

    % position
    e3vec = squeeze(Rmat(:,3,1:i));
    if i == 2
        pmat(:,i) = pmat(:,i-1) + squeeze(Rmat(:,3,i))*ds;
    else
        pmat(:,i) = Simpson_vec_int(e3vec,ds);
    end        
end