%
% plot normal flux jump across the interface
%

hold on

allq = [];
% num of points in each segment
np = 10;
pos = [];
lastx = 0;
for p = 1:length(polis)
   for s = 1:length(polis{p})-1
      p1 = polis{p}(s,:);
      p2 = polis{p}(s+1,:);
      dp = p1-p2;
      pos = [linspace(p1(1),p2(1),np)' linspace(p1(2),p2(2),np)'];
      % flux
      qp = FluxosX( X, T, levelSet, h, hE, pos, 1, opts.tolerance );
      qn = FluxosX( X, T, levelSet, h, hE, pos, 0, opts.tolerance );
      % delta flux
      dq = qp - qn;
      % normal to the segment (linear levelset assumed!)
      nn = [-dp(2) dp(1)];
      nn = nn/norm(nn);
      % flux jump
      normalq = sum(dq.*repmat(nn,size(dq,1),1),2);
      % plot stuff
      di = pos - repmat(p1,np,1);
      di = di.^2;
      di = sqrt(sum(di,2));
      allq = [allq; normalq];
      plot( di+lastx, normalq, 'linewidth', 2 )
      lastx = lastx + di(end);
      if length(polis{p}) < 15
         txt = sprintf('%i',s);
         text(p1(1),0,txt, ...
            'color', 'r', 'FontWeight', 'bold', 'fontsize',11, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
      end
   end
end
% plot( [xmin xmax],[ymin ymin], '--k')
axis tight


xmin = 0;
xmax = max(di)+lastx;

plot([xmin xmax], [0 0], ':k')

max(allq)
mean(allq)