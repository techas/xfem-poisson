function MakeMovieGasNatural( out, X, T, time, graphiks );
% MAKEMOVIEGASNATURAL to make a movie from the output of a computation for
% the GasNatural case
%
% syntax: MakeMovieGasNatural( out, graphiks )
%
%  out: structured array containing the field 't', 'levelSet' and
%       'phreaLevel'
%  graphiks: structured array for graphical options

% R. Cottereau 04/2008

% get rid of repeated time steps
ind = find( diff( out.t ) ~= 0 );
out.t = out.t( ind );
out.levelSet = out.levelSet( :, ind );
out.waterTable = out.waterTable( ind );

% create the time vector
t = [ out.t(1) : graphiks.TimeStep : out.t(end) ];

% constants
Ne = size( out.levelSet, 1 );
Nt = length( t );

% initializations
l = zeros( Ne, Nt );

% interpolate the calculated results
for i1 = 1:Ne
    l( i1, : ) = interp1( out.t, out.levelSet( i1, : ), t );
end
wt = interp1( out.t, out.waterTable, t );

% creating the movie structure
graphiks.Type = 'Front';
graphiks.Step = 1;
graphiks.Movie = 1;
for i1 = 1:Nt
    M(i1) = RealTimeGraphics( 1, graphiks, i1, X, T, l( :, i1 ), wt(i1) );
end

% save movie on file
movie2avi(M,graphiks.MovieFileName);
