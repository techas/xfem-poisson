function [xi,eta] = invMap( Xe, px, py )
% [xi,eta] = invMap( Xe, px, py )
%
% Map xy-points inside an triangular element (cartesian space) 
% in xieta-points inside the reference triangle.
%
% INPUT
%   Xe      nodal coords of the element (cart space)
%   px,py   points to map (cart space)
%
% OUTPUT
%   xi,eta  mapped points (reference space)
%
x = Xe(:,1);
y = Xe(:,2);
a1 = x(2)*y(3) - x(3)*y(2);
a2 = x(3)*y(1) - x(1)*y(3);
b1 = y(2) - y(3);
b2 = y(3) - y(1);
c1 = x(3) - x(2);
c2 = x(1) - x(3);
M = [[1;1;1] Xe];
A = det( M );
xi = (a1 + b1*px + c1*py)/A;
eta = (a2 + b2*px + c2*py)/A;
