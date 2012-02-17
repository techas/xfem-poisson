%
%
%
theSide = 0;
%
hold on
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
      %qp = FluxosX( X, T, levelSet, h, hE, pos, 1, opts.tolerance );
      if analytic
         % flux - analytic 
         r0 = 0.5;
         normalq = 2*r0*ones(np,1);
      else
         qn = FluxosX( X, T, levelSet, h, hE, pos, theSide, opts.tolerance );
         % delta flux
         % dq = qp - qn;
         dq = qn;
         % normal to the segment (linear levelset assumed!)
         nn = [-dp(2) dp(1)];
         nn = nn/norm(nn);
         % flux jump
         normalq = sum(dq.*repmat(nn,size(dq,1),1),2);
      end
      % plot stuff
      di = pos - repmat(p1,np,1);
      di = di.^2;
      di = sqrt(sum(di,2));
      plot( di+lastx, normalq, 'linewidth', 2, 'color', col )
      %
      e = Ei{1}(s);
      ix = find(enrichedElements==e);
      text( mean(di)+lastx,mean(normalq), sprintf('%f', debug(ix)));
      %
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
