%
% plot flux
%
plotLS( X, levelSet );
% plot( [xmin xmin xmax xmax xmin], [ymin ymax ymax ymin ymin], 'k', 'linewidth', 2 )
trimesh( T, X(:,1), X(:,2), 'color', 0.4*[1 1 1 ] )
qH = qH./max(sqrt(qH(:,1).^2+qH(:,2).^2));


onlyEnriched = 1;
if onlyEnriched 
   ngp = 4;
   ix = zeros(size(Xpg,1),1);
   nee = length(enrichedElements);
   p = nee*ngp - 1;   
   hhh=quiver( Xpg(end-p:end,1), Xpg(end-p:end,2), ...
      qH(end-p:end,1), qH(end-p:end,2), 1, 'color', 'k' );
else
   hhh=quiver( Xpg(:,1), Xpg(:,2), qH(:,1), qH(:,2), 1, 'color', 'k' );
end

