function Y = ad_se3(X)

% basis elements for se(3)
X_se3 = zeros(4,4,6);
X_se3(:,:,1) = [0,0,0,0;0,0,-1,0;0,1,0,0;0,0,0,0];
X_se3(:,:,2) = [0,0,1,0;0,0,0,0;-1,0,0,0;0,0,0,0];
X_se3(:,:,3) = [0,-1,0,0;1,0,0,0;0,0,0,0;0,0,0,0];
X_se3(:,:,4) = [0,0,0,1;0,0,0,0;0,0,0,0;0,0,0,0];
X_se3(:,:,5) = [0,0,0,0;0,0,0,1;0,0,0,0;0,0,0,0];
X_se3(:,:,6) = [0,0,0,0;0,0,0,0;0,0,0,1;0,0,0,0];

if size(X,2) == 4 % --> input is 4 x 4 matrix
    Y = zeros(6);
    for i=1:6
        Y(:,i) = vect(X*X_se3(:,:,i) - X_se3(:,:,i)*X);
%         eval(['Y(:,',num2str(i),')=vect(X*X',num2str(i),'-X',num2str(i),'*X);'])
    end
elseif size(X,2) == 1 % --> input is 6 x 1 vector
    Y = [matr(X(1:3)),zeros(3);matr(X(4:6)),matr(X(1:3))];
end
