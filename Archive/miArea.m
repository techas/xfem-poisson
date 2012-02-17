function A = miArea( nodes )
% A = miArea( nodes )
%
% compute the area of a triangle
%
% INPUT
%   nodes   nodal coords of the triangle. size( nodes ) = [3,2]
%
% OUTPUT
%   A       area
%
n1 = nodes(1,:);
n2 = nodes(2,:);
n3 = nodes(3,:);
v1 = [(n3 - n2) 0];
v2 = [(n1 - n2) 0];
A1 = cross( v1, v2 );
A = A1(3) * 0.5;
