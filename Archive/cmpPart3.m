%
% cmpPart3
%
enrichedNodes = debugData.enrichedNodes;
enrichedElements = debugData.enrichedElements;

if isfield( debugData, 'polis' )
   polis = debugData.polis;
else
   [ Seg, SegsBnd ] = CrossedSegments( T, enrichedElements, levelSet, opts.tolerance );
   [ polis, dummy, Segi ] = MakePoligonalFromSegments( X, T, SegsBnd, Seg, levelSet );
end

if isfield( debugData, 'Ei' )
   Ei = debugData.Ei;
else
   [ Seg, SegsBnd ] = CrossedSegments( T, enrichedElements, levelSet, opts.tolerance );
   [ dummy, Ei, Segi ] = MakePoligonalFromSegments( X, T, SegsBnd, Seg, levelSet );
end

if isfield( debugData, 'theNormals' )
   theNormals = debugData.theNormals;
end


xmin = min(X(:,1));
xmax = max(X(:,1));
ymin = min(X(:,2));
ymax = max(X(:,2));

analytic = 0;
%% ---
figure(2), clf
figure(3), clf
% figure(5), clf
%

%% FEM
hh=figure(4); clf, hold on
set(hh,'name','fem');
load datafem.mat
plotFlux, axis equal tight, axis off
% plot(X(killedNodes,1),X(killedNodes,2),'or')
set(hhh,'linewidth', 1.0, 'Clipping', 'off', 'AutoScaleFactor', 3 )
%
figure(6), clf
subplot(311)
plotNormalFluxJump



%% XFEM
hh=figure(2); clf, hold on
set(hh,'name','xfem');
load dataxfem.mat
plotFlux, axis equal tight, axis off
set(hhh,'linewidth', 1.0, 'Clipping', 'off', 'AutoScaleFactor', 3 )
% plot(X(killedNodes,1),X(killedNodes,2),'or')
% figure(5), hold on
% col = 'b'; 
% cmpPlotFluxInterface
%
figure(6)
subplot(312)
plotNormalFluxJump


%% XFEMP
hh=figure(3); clf, hold on
set(hh,'name','xfem+');
load dataxfemp.mat
plotFlux, axis equal tight, axis off
set(hhh,'linewidth', 1.0, 'Clipping', 'off'  )
% plot(X(killedNodes,1),X(killedNodes,2),'or')
% figure(5), hold on
% col = 'r'; 
% cmpPlotFluxInterface
%
figure(6)
subplot(313)
plotNormalFluxJump



%
% figure(5), hold on
% analytic = 1;
% col = 'k'; 
% cmpPlotFluxInterface

