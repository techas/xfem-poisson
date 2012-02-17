function [ti,xi] = inContactWithElements( X, T, t0 );
% INCONTACTWITHELEMENTS to select the nodes and elements in contact with a 
% given group of elements
%
%  syntax [ti,xi] = inContactWithElements( X, T, t0 );
%
%  T  : connectivity matrix of the entire mesh
%  t0 : elements of the initial group
%
%  ti : elements in t0 or in contact with t0
%  xi : nodes in ti

% R. Cottereau 03/2008

x0 = unique(T(t0,:));
ti = [];
for i1=1:size(x0);
    [t1,j1] = find( T == x0(i1) );
    ti = union( ti, t1 );
end
xi = unique(T(ti,:));