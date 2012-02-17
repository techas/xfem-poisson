%
% plot interface discretization
%
plotLS( X, levelSet );
plot( [xmin xmin xmax xmax xmin], [ymin ymax ymax ymin ymin], 'k', 'linewidth', 2 )
trimesh( T, X(:,1), X(:,2), 'color', 0.4*[1 1 1 ] )

for p = 1:length(polis)
   for s = 1:length(polis{p})-1
      p1 = polis{p}(s,:);
      p2 = polis{p}(s+1,:);
      m = (p1+p2)/2;
      plot( [p1(1) p2(1)],[p1(2) p2(2)],'-ob','linewidth', 2);
      if length(polis) == 1
         txt = sprintf('%i',s);
      else
         txt = sprintf('%i-%i',p,s);
      end
      text(m(1),m(2),txt, ...
         'color', 'k', 'FontWeight', 'bold', 'fontsize',11, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
            % normal to the interface
%       if exist('theNormals','var')
%          if length(polis) == 1
%             ix = Ei{p}(s);
%             ix2 = vectorFind( enrichedElements, ix );
%             xx = [m; m + theNormals(ix2,:)*0.2];
%             plot(xx(:,1),xx(:,2),'-r')
%          end
%       end

   end
end
axis equal
