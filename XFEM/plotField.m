function SM = plotField( X, field, xsamples, ysamples, limits )
% SM = plotField( X, field, xsamples, ysamples, limits )
% 
% plot a scalar field
%
% limits [x_lo x_up y_lo y_up]
%
if nargin < 4
   xsamples = 300;
   ysamples = 100;
end
if nargin <= 4
   x_lo = min( X(:,1) );     x_up = max( X(:,1) );  
   y_lo = min( X(:,2) );     y_up = max( X(:,2) ); 
end
if nargin == 5
   x_lo = limits(1);     x_up = limits(2);  
   y_lo = limits(3);     y_up = limits(4);   
end
xlin = linspace( x_lo, x_up, xsamples );
ylin = linspace( y_lo, y_up, ysamples );
[sX,sY] = meshgrid( xlin, ylin );
SM = griddata( X(:,1), X(:,2), field, sX, sY, 'linear', {'QJ'} );
set( gca, 'UserData', SM )
pcolor( sX, sY, SM )
shading flat
% axis equal tight
