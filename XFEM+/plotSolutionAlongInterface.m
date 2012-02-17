%
% plot solution along interface
%


% num of points in each segment
np = 10;
pos = [];
lastx = 0;
for p = 1:length(polis)
   for s = 1:length(polis{p})-1
      % coords of the segment
      p1 = polis{p}(s,:);
      p2 = polis{p}(s+1,:);
      pos = [linspace(p1(1),p2(1),np)' linspace(p1(2),p2(2),np)'];
      % Element data
      idelem = Ei{p}(s);
      Te = T(idelem,:);
      Lse = levelSet(Te);
      he = h(Te);
      % shape functions at pos
      Xe = X(Te,:);
      [xi,eta] = invMap( Xe, pos(:,1), pos(:,2) );
      [N,Nxi,Neta] = shapeFunction( triangle, nPoints, [xi,eta] );
      % std solution
      sol = N*he;
      % enriched solution
      Tee = vectorFind( enrichedNodes, Te );
      hee = hE(Tee);
      % Ridge
      R = N*abs(Lse) - abs(N*Lse);
      M = N.*repmat(R,1,size(N,2));
      % complete sol
      sol = sol + M*hee;
      % plot stuff
      di = pos - repmat(p1,np,1);
      di = di.^2;
      di = sqrt(sum(di,2));
      plot( di+lastx, sol, 'linewidth', 2 )
      lastx = lastx + di(end);
      % plot segments
      %          plot( p1(1), sol(1), '^r', 'MarkerFaceColor','r','MarkerSize',5)
      if length(polis{p}) < 15
         txt = sprintf('%i',s);
         text(p1(1),sol(1),txt, ...
            'color', 'r', 'FontWeight', 'bold', 'fontsize',11, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
      end
   end
end
axis tight
