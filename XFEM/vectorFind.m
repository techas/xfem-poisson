function ix = vectorFind( vector, values )
% ix = vfind( vector, values )
%
% The same as find but looks for several values in a vector.
% Returns a vector of the same size of values, with the indices of the 
% values in the vector.
%
% INPUT
%   vector   vector to look on
%   values   values to look for
%
% OUTPUT
%  ix      index of values in vector.
%            vector(ix) = values;
%
% Example:
%     vec = [5 3 6 7 4 2]
%     val = [5 4]
%     ix = vfind( vec, val )
%     ans = [1 5]
%
ix = zeros( 1, length( values ) );
for I = 1:length( values )
   newix = find( vector == values(I), 1 );
   if isempty( newix )
      error( 'Index not found' )
   end
   ix(I) = newix;
end
