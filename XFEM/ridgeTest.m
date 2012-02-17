close all
clear all
clc
load ridgetest.mat

levelSet = [1 -1 1]';
% NB: with levelSet = [1 -1 1]', the ridge and the normalised rigde are the
% same !!

N = [1-X(:,1)-X(:,2) X(:,1) X(:,2)];
Nxi = ones( size(X,1), 1 ) * [-1 1 0];
Neta = ones( size(X,1), 1 ) * [-1 0 1];

R = N * abs( levelSet ) - abs( N * levelSet );
Rxi = Nxi * abs( levelSet ) - sign( N * levelSet ) .* ( Nxi * levelSet );
Reta = Neta * abs( levelSet ) - sign( N * levelSet ) .* ( Neta * levelSet );
gradR = [Rxi Reta];

NRxi = Nxi .* repmat( R, [1 3]) + N .* repmat( Rxi, [1 3]);
NReta = Neta .* repmat( R, [1 3]) + N .* repmat( Reta, [1 3]);
gradNR = [NRxi NReta];

%trisurf( T, X(:,1), X(:,2), NRxi(:,3) )

DibujaVec( X, X, T, gradNR, levelSet, 0, 'k-', 0.05 );