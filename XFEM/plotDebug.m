function plotDebug( Xe, Xp, Ns )
%
%
%
figure( 33 )
hold on
order = [1  2  3  1];          
%
plot( Xp(order,1), Xp(order,2), '-b' );
for I = 1:length( Xp )
   xy = ( Xp(I,:) + mean( Xp ) )/2;
   text( xy(1), xy(2), int2str( I ), ...
         'FontSize', 10, ...
         'Color', 'b' );
end
%
plot( Xe(order,1), Xe(order,2), '-k', 'linewidth', 2 );
for I = 1:length( Xe )
   xy = ( Xe(I,:) + mean( Xe ) )/2;
   text( xy(1), xy(2), int2str( I ), ...
         'FontSize', 12, ...
         'Color', 'k', ...
         'Fontweight', 'bold' );
end

if nargin == 3
   xy = Ns*Xe;
   plot( xy(:,1), xy(:,2 ), 'xr' )
end