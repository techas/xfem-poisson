function inter = intersection( pos, pes )
% inter = intersection( pos, pes )
%
% Find the point where a linear function is zero. The linear function is
% defined by the points. 
%
% pes1 |
%      |   .\
%      |   . \
%      |___.__\__pos2____
%      |  pos1 \  .     
%      |        \ .
% pes2 |         \.
%      |
%
% INPUT
%   pos     [x1 x2]  
%   pes     [y1 y2]  
%
% OUTPUT
%   inter   the x where the line cross the x-axis
%

% when the two points are the same check that the weights are almost zeros
% and take the intersection as the point entered
if pos(1,1) == pos(2,1) & pos(1,2) == pos(2,2)
%     if 1e-12 > abs(pes(1,1))
%         inter = pos(1,:);
%         return
%     else
%         diff(pes)
%         error('for the same point the weight should be zeros')
%     end
        inter = pos(1,:);
        return
end

if prod( pes ) > 0
   error( 'interseccion: There is no zero in the interval' )
end
%
alpha = pes(2) / (pes(2) - pes(1));
beta = pes(1) / (pes(2) - pes(1));
inter = pos(1,:) * alpha - pos(2,:) * beta;
