function gmat = cal_pose_shear(xiv,g0,N,ds)
%
% calculate pose given \xi (extensible rod)
%
% - written by Jin Seob (Jesse) Kim

gmat = zeros(4,4,N);
gmat(:,:,1) = g0;

for i = 2:N
    xi = 1/2*(xiv(:,i-1) + xiv(:,i));
    gmat(:,:,i) = gmat(:,:,i-1)*exp_se3(ds*matr(xi));
end
