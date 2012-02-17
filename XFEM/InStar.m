function in = InStar(Xs,X,Ts)
% INSTAR returns the points in Xs that are contained in the star defined by
% Ts and X
%
%  syntax: In = InStar(Xs,X,Ts)
%
% X : nodal coordinates of the mesh
% Ts: connectivity matrix of the star
% Xs: points to be checked
%
% In: vector of logicals of the same size as Xs
%
% NB: this function is a loop on the elements of Ts, calling inpolygon

% R. Cottereau 04/2008

in = zeros(size(Xs,1),1);

for i1=1:size(Ts,1)
    in = in | inpolygon( Xs(:,1), Xs(:,2), X(Ts(i1,:),1), X(Ts(i1,:),2) );
end