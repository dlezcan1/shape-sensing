function y = costfn_theta0(theta0,wv,pos_img,ds,scale_th)

N = size(wv,2);

% configuration
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


y_data = pos_img(:,2);
z_data = pos_img(:,1);

y = extract_data_pos(y_data,z_data,pmat);

