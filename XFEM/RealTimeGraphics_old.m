function m = RealTimeGraphics( handl, graphiks, t, X, T, ...
                          levelSet, ph, h, hE, pospg, N, Nxi, Neta, opts );
% to draw Real-Time graphics

% R. Cottereau 03/2008
if (graphiks.Step~=Inf) & (graphiks.Movie)
    figure(handl);
    set(gca,'PlotBoxAspectRatio',[1 1 1]);
    set(gca,'nextplot','replacechildren');
    switch graphiks.Type
        % EVOLUTION OF THE WATER FRONT
        case 'Front'
            fill([0 -0 1 1], [-1 ph ph -1], 'b','FaceAlpha',.5);
            hold on;trimesh(T,X(:,1),X(:,2),'Color','k');
            plotLS( X, levelSet );
            if graphiks.GasNatural
                hold on;plotTubes;
                view(2);set(gca,'XLim',[0 1],'YLim',[-1 0]);
            end
        % EVOLUTION OF FLUJOS (time consuming and very unstable!!)
        case 'Flux'
            XbelowPhreatic = find( X(:,2) < ph );
            BelowPhreatic = ismember( T, XbelowPhreatic );
            Tb = T( find( sum(BelowPhreatic,2) == 3 ), : );
            [Xpg,Fpg] = Flujos(h,hE,X,Tb,pospg,N,Nxi,Neta,levelSet,opts);
%            Fe = repmat([-1000/(505+.495)],[size(Fpg,1) 1]);
%            max(abs(Fpg(:,2)-Fe))./max(abs(Fpg(:,2)+Fe))
            if graphiks.GasNatural
                fill([0 -0 1 1], [-1 ph ph -1], 'b','FaceAlpha',.5);
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

centerTube1 = [.3 -.5];
radiusTube1 = .03;
numberNodesTube1 = 10;
thetaTube1 = 2*pi*[numberNodesTube1:-1:1]'/numberNodesTube1;
nodeTube1 = repmat(centerTube1,numberNodesTube1,1) + ...
                            radiusTube1*[cos(thetaTube1) sin(thetaTube1)];
fill3(nodeTube1(:,1),nodeTube1(:,2),0.01*ones(size(nodeTube1,1),1), ...
          'w','FaceAlpha',1);

centerTube2 = [.69 -.69];
radiusTube2 = .03;
numberNodesTube2 = 30;
thetaTube2 = 2*pi*[numberNodesTube2:-1:1]'/numberNodesTube2;
nodeTube2 = repmat(centerTube2,numberNodesTube2,1) + ...
                            radiusTube2*[cos(thetaTube2) sin(thetaTube2)];
fill3(nodeTube2(:,1),nodeTube2(:,2),0.01*ones(size(nodeTube2,1),1), ...
          'w','FaceAlpha',1);
