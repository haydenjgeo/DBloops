function [d,r2] = ellipsoid_distance(x,y,z,p)
% Distance of points projected onto an ellipsoid.
%
% Input arguments:
% x, y, z:
%    co-ordinates of data points whose distance from the ellipsoid to measure
% p:
%    parameters of the ellipsoid expressed in implicit form

% Copyright 2011 Levente Hunyadi

[center,radii,~,R] = ellipsoid_im2ex(p);
[xp,yp,zp] = ellipsoidfit_residuals(x,y,z, center,radii,R);
d = (x-xp).^2 + (y-yp).^2 + (z-zp).^2;
r2= 1- sum((x-xp).^2 + (y-yp).^2 + (z-zp).^2)./sum((x-mean(x)).^2 + (y-mean(y)).^2 + (z-mean(z)).^2);