function I=Simpson_vec_int(f,h)
% This function returns the integral by Simpson's rule.
% f is m x n, m is the dimension of a vector, and n is the dimension of a
% parameter.
% - written by Jin Seob Kim.

n=length(f)-1;
if n==1
    fprintf('Data has only one interval')
    return;
end
if n==2
    I=h/3*(f(:,1)+4*f(:,2)+f(:,3));
    return;
end
if n==3
    I=3/8*h*(f(:,1)+3*f(:,2)+3*f(:,3)+f(:,4));
    return;
end

I=0;
if 2*floor(n/2)~=n
    I=3/8*h*(f(:,n-2)+3*f(:,n-1)+3*f(:,n)+f(:,n+1));
    m=n-3;
else
    m=n;
end

I=I+(h/3)*(f(:,1)+4*sum(f(:,2:2:m),2)+f(:,m+1));
if m>2
    I=I+(h/3)*2*sum(f(:,3:2:m),2);
end