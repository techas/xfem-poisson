function bool = crossedByLevelSet( LS, T, tol )
% [bool,nls] = crossedByLevelSet( LS, T, tol )
% 
% returns if the element has to be enriched or not
%
% INPUT
%   LS   nodal level set values
%   T    connectivity matrix
%   tol  tolerance
%
% OUTPUT
%   bool                       true: the element has to be enriched
%                              false: std element
%   numberOfCrossingLevelSets  
%   firstCrossingLevelSet      
%

% R. Cottereau 04/2008 modified from S. Zlotnik
% OLD CALL TO FUNCTION:
%function [bool,numberOfCrossingLevelSets,firstCrossingLevelSet] = ...
%   crossedByLevelSet( LS, T, tol )
% it cannot be used as it is for several level sets

if nargin < 2
   tol = 0;
end
LST = LS(T);
if size(T,1)==1
    LST = LST';
end
% All negative level sets
N = all( LST <= tol, 2 );
% All positive level sets
P = all( LST >= -tol, 2 );
% Splited?
S = ~N & ~P ;
% The first crossing level set 
%C = find( S == 1, 1 );
%if isempty( C )
%   C = 0;
%end
% M = k, then use material k
% M = 0, then use several material or use the last material
%M = find( P == 1, 1 );
%if isempty( M )
%   M = 0;
%end
%bool = (M == 0 & C > 0) | (C < M & C > 0);
%numberOfCrossingLevelSets = length( find( S ) );
%firstCrossingLevelSet = C;
bool = S;