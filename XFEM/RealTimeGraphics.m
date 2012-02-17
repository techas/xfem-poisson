function m = RealTimeGraphics( handl, graphiks, t, X, T, ...
                          levelSet, ph, h, hE, pospg, N, Nxi, Neta, opts );
% to draw Real-Time graphics

% R. Cottereau 03/2008
if (graphiks.Step~=Inf) & (graphiks.Movie)
    figure(handl);
%    set(gca,'PlotBoxAspectRatio',[1.4 1.2 1]);
    set(gca,'PlotBoxAspectRatio',[1.2 1.2 1]);
    set(gca,'nextplot','replacechildren');
    switch graphiks.Type
        % EVOLUTION OF THE WATER FRONT
        case 'Front'
            hold off;trimesh(T,X(:,1),X(:,2),'Color','k');
            hold on;plotLS( X, levelSet );
            if graphiks.GasNatural
                hold on;plotTubes;
                view(2);set(gca,'XLim',[-.2 1.2],'YLim',[-1.2 0]);
                hold on; fill3([-.2 -.2 1.2 1.2], [-1.2 ph ph -1.2], ...
                                      0.01*[1 1 1 1],'b','FaceAlpha',.5);
                box on;
            end
            if graphiks.Paper
                hold on;plotTubesPaper;
                view(2);set(gca,'XLim',[0 1.2],'YLim',[-1.2 0]);
                hold on; fill3([-.2 -.2 1.2 1.2], [-1.2 ph ph -1.2], ...
                                      0.01*[1 1 1 1],'b','FaceAlpha',.5);
                box on;
            end
        % EVOLUTION OF FLUJOS (time consuming and very unstable!!)
        case 'Flux'
            XbelowPhreatic = find( X(:,2) < ph );
            BelowPhreatic = ismember( T, XbelowPhreatic );
            Tb = T( ( sum(BelowPhreatic,2) == 3 ), : );
            [Xpg,Fpg] = Flujos(h,hE,X,Tb,pospg,N,Nxi,Neta,levelSet,opts);
%            Fe = repmat([-1000/(505+.495)],[size(Fpg,1) 1]);
%            max(abs(Fpg(:,2)-Fe))./max(abs(Fpg(:,2)+Fe))
            if graphiks.GasNatural
                fill([-.3 -.3 1.3 1.3], [-1.3 ph ph -1.3], ...
                                                       'b','FaceAlpha',.5);
            end
            hold on
            DibujaVec( Xpg, X, T, Fpg, levelSet );
            if graphiks.GasNatural
                hold on;plotTubes;
                view(2);set(gca,'XLim',[0 1],'YLim',[-1 0]);
            end

        % EVOLUTION OF THE LEVEL SET
        case 'Level Set'
            trisurf(T,X(:,1),X(:,2),levelSet);
            view(3);caxis([-1 1]*opts.cutLevelSet);colorbar;
            set(gca,'ZLim',[-1 1]*opts.cutLevelSet);
            if graphiks.GasNatural
                set(gca,'XLim',[0 1],'YLim',[-1 0]);
            end
        % EVOLUTION OF THE PIEZOMETRIC HEIGHT
        case 'h'
%            he = zeros(size(h));
%            ind = find(levelSet > 0);he(ind) = 1000/(505+.495)*X(ind,2);
%            ind = find(0 > levelSet);he(ind) = 1+1/(505+.495)*(X(ind,2)-1);
%            trisurf(T,X(:,1),X(:,2),h-he);
             trisurf(T,X(:,1),X(:,2),h);
             set(gca,'XLim',[0 1],'YLim',[0 1],'ZLim',[0 1]);
        case 'Gradient Norm q'
            [Xpg,Fpg] = Flujos(h,hE,X,T,pospg,N,Nxi,Neta,levelSet,opts);
    end
    % MAKE MOVIE
    if graphiks.Movie
        m = getframe;
    end
else
    m = {0};
end

%==========================================================================
function plotTubes

centerTube1 = [.25 -.5];
radiusTube1 = .025;
numberNodesTube1 = 15;
thetaTube1 = 2*pi*[numberNodesTube1:-1:1]'/numberNodesTube1;
nodeTube1 = repmat(centerTube1,numberNodesTube1,1) + ...
                            radiusTube1*[cos(thetaTube1) sin(thetaTube1)];
fill3(nodeTube1(:,1),nodeTube1(:,2),0.01*ones(size(nodeTube1,1),1), ...
          'w','FaceAlpha',1);

centerTube2 = [.7 -.84];
radiusTube2 = .05;
numberNodesTube2 = 28;
thetaTube2 = 2*pi*[numberNodesTube2:-1:1]'/numberNodesTube2;
nodeTube2 = repmat(centerTube2,numberNodesTube2,1) + ...
                            radiusTube2*[cos(thetaTube2) sin(thetaTube2)];
centerTube3 = [.7 -.76];
radiusTube3 = .015;
numberNodesTube3 = 10; % has to be even
thetaTube3 = 2*pi*[numberNodesTube3:-1:1]'/numberNodesTube3;
nodeTube3 = repmat(centerTube3,numberNodesTube3,1) + ...
                            radiusTube3*[cos(thetaTube3) sin(thetaTube3)];
X31 = nodeTube3(1,1);
Y31 = nodeTube3(1,2);
RefineBox23 = .008;
box23 = [X31 Y31
         X31 Y31-(centerTube3(2)-centerTube2(2))
         X31-2*radiusTube3 Y31-(centerTube3(2)-centerTube2(2))
         X31-2*radiusTube3 Y31];
tube2Inbox23 = inpolygon( nodeTube2(:,1), nodeTube2(:,2), ...
                                                  box23(:,1), box23(:,2) );
tube2Inbox23first = find( tube2Inbox23, 1, 'first' );
tube2Inbox23last = find( tube2Inbox23, 1, 'last' );
tube3Inbox23first = numberNodesTube3/2+1;
tube3Inbox23last = numberNodesTube3;
lbox23 = nodeTube3( tube3Inbox23first, 2)-nodeTube2( tube2Inbox23first, 2);
box23base = linspace( 0, lbox23, ceil(lbox23/RefineBox23) )';
nbox23 = length(box23base)-1;
box23up = repmat(nodeTube3( tube3Inbox23first, : ), [nbox23 1] ) - ...
          [ zeros( nbox23, 1 ) flipud( box23base(2:end,1) ) ];
box23down = repmat(nodeTube3( 1, : ), [nbox23 1] ) - ...
          [ zeros( nbox23, 1 ) box23base(2:end,1) ];
indTube3 = [tube3Inbox23first:tube3Inbox23last 1];

nodeTube2 = [nodeTube2( 1:tube2Inbox23first-1, : )
             box23up
             nodeTube3( indTube3, : )
             box23down
             nodeTube2( tube2Inbox23last+1:end, : )];

fill3(nodeTube2(:,1),nodeTube2(:,2),0.01*ones(size(nodeTube2,1),1), ...
          'w','FaceAlpha',1);
      
%==========================================================================
function plotTubesPaper

centerTube1 = [.9 -.7];
radiusTube1 = .05;
numberNodesTube1 = 20;
thetaTube1 = 2*pi*[numberNodesTube1:-1:1]'/numberNodesTube1;
nodeTube1 = repmat(centerTube1,numberNodesTube1,1) + ...
                            radiusTube1*[cos(thetaTube1) sin(thetaTube1)];
fill3(nodeTube1(:,1),nodeTube1(:,2),0.01*ones(size(nodeTube1,1),1), ...
          'w','FaceAlpha',1);

      