%
% plot enriched elements and nodes
%
plotLS( X, levelSet );
plot( [xmin xmin xmax xmax xmin], [ymin ymax ymax ymin ymin], 'k', 'linewidth', 2 )
trimesh( T, X(:,1), X(:,2), 'color', 0.4*[1 1 1 ] )
ix = [1:3 1];
hold on;
for eid = enrichedElements'
   Te = T(eid,:);
   Xe = X(Te,:); 
   plot( Xe(ix,1),Xe(ix,2),'-k','linewidth', 2);
   m = mean(Xe);
   txt = sprintf('%i',eid);
   text(m(1),m(2),txt,...
      'color', 'k', 'FontWeight', 'bold', 'fontsize',13, ...
      'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
end
for nid = enrichedNodes'
   txt = sprintf('%i',nid);
   text(X(nid,1),X(nid,2),txt,...
      'color', 'r', 'FontWeight', 'bold', 'fontsize',13, ...
      'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
end   
axis equal
