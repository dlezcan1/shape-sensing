function A = Acurve_param(y,x)
% 
% calculate the area defined as y(x), where the curve is parameterized as
% y(s) and x(s).
% Use the simplest form:
% A = \sum_{i=2}^n y_i (x_i - x_{i-1}) (i = 1,2, \cdots, n)
% analytical form = \int_0^L y(s) dx(s)/ds ds
%
% - written by Jin Seob Kim

yv = y(2:end);
dx = diff(x);
A = sum(yv.*dx);