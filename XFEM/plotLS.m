function hh = plotLS( X, levelSet, filled, xsamples, ysamples, limits  )
%
% plot the level set
%
color = [.7 .7 .7];
%color = [eye(3); 1 1 0 ; 1 0 1 ; 0 1 1 ; 1 1 1 ; 0 0 0];
if size( levelSet, 1 ) < size( levelSet, 2 )
   levelSet = levelSet';
end
if nargin < 3;
   filled = 1;
end
if nargin < 4
   xsamples = 300;
   ysamples = 300;
end
if nargin <= 6
   x_lo = min( X(:,1) );     x_up = max( X(:,1) );  
   y_lo = min( X(:,2) );     y_up = max( X(:,2) ); 
end
if nargin == 6
   x_lo = limits(1);     x_up = limits(2);  
   y_lo = limits(3);     y_up = limits(4);   
end
xlin = linspace( x_lo, x_up, xsamples );
ylin = linspace( y_lo, y_up, ysamples );
[sX,sY] = meshgrid( xlin, ylin );
hold on
hh = zeros( 1, size( levelSet, 2 ) );
for I = size( levelSet, 2 ):-1:1
   sF = griddata( X(:,1), X(:,2), levelSet(:,I), sX, sY);
   if filled
%      [C1,H1] = contourf( sX, sY, sF, [0 Inf], 'k-' );
      [C1,H1] = contourf( sX, sY, sF, [0 0], 'Tag', 'levelSetTag' );
      pp1 = get( H1, 'Children' );
%      set( pp1(2), 'FaceColor', color(I,:) )
%      set( pp1(1), 'FaceAlpha',1);
%      set( pp1(1), 'FaceAlpha',.5);
   else
%      [C,H1] = contour3( sX, sY, sF, [0 0], 'g-');
      [C,H1] = contour( sX, sY, sF, [0 0], 'k', 'Tag', 'levelSetTag' );
      set(H1,'LineWidth',2);
   end
%   hh(I) = H1;
end
% hold off


