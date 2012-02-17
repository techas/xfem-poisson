% my figure
clc
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

if useEnrichment == 1
   txt = 'XFEM';
else
   txt = 'XFEM+';
end
%%
triangle = 2;
nPoints = 3;
%%
myfh=figure(1); clf; 
set(myfh, 'name', txt)
msp=2;
nsp=3;
isp=1;
xmin = min(X(:,1));
xmax = max(X(:,1));
ymin = min(X(:,2));
ymax = max(X(:,2));
% - 1 flux ---------------------------------------------------------------
% fh(isp)=subplot(msp,nsp,isp); isp = isp+1; hold on;
plotFlux
% title( 'flux' )
% %
% % energy computation
% u = [ h; hE ];
% Nu = length(u);
% eH = u' * Ktot(1:Nu,1:Nu) * u;
% xlabel(sprintf('||u||=%1.5g',eH))
% 
% % - 2 enr elem and nodes -------------------------------------------------
% fh(isp)=subplot(msp,nsp,isp); isp = isp+1; hold on;
% title( 'Enr nodes and elems' )
% plotEnrichedElements
% 
% % - interface segments & normals -----------------------------------------
% fh(isp)=subplot(msp,nsp,isp); isp = isp+1; hold on;
% title( 'Discr Interphase' )
% plotInterfaceDiscretization
% 
% % -normal flux acreoss the interface -------------------------------------
% subplot(msp,2,3); hold on;
% title('normal flux accros the interface')
% plotNormalFluxJump
% 
% % - solution along the interface -----------------------------------------
% subplot(msp,2,4); isp = isp+1; hold on;
% title('solution along the interface')
% plotSolutionAlongInterface
% 
% % ------------------------------------------------------------------------
% if ~isempty( fh ) 
%    if size(X,2) == 3
%       hlink = linkprop( fh(:), {'xlim','ylim','zlim','CameraPosition','CameraUpVector'} );
%    else
%       hlink = linkprop( fh(:), {'xlim','ylim'} );
%    end
%    setappdata( fh(1), 'graphics_linkprop', hlink ); 
% end
% 
% % - flux jump * N --------------------------------------------------------
% % num of points in each segment
% if useEnrichment == 2
%    % gauss quad
%    xi = [-1/sqrt(3); 1/sqrt(3)];
%    wi = [1; 1];
%    Ne = [(1-xi)/2 (1+xi)/2];
%    Ne_xi = [-0.5 0.5];
%    for p = 1:length(polis)
%       for s = 2:length(polis{p})-2
%          % two consecutive segments
%          p1 = polis{p}(s,:);
%          p2 = polis{p}(s+1,:);
%          p3 = polis{p}(s+2,:);
%          Xe1 = [p1;p2];
%          Xe2 = [p2;p3];
%          pos = [Ne*Xe1; Ne*Xe2];
%          % fluxes
%          qp = FluxosX( X, T, levelSet, h, hE, pos, 1, opts.tolerance );
%          qn = FluxosX( X, T, levelSet, h, hE, pos, 0, opts.tolerance );
%          % delta flux
%          dq = qp - qn;
%          % normal to each segment (linear levelset assumed!)
%          dp = p1-p2;
%          nn = [-dp(2) dp(1)];
%          nn1 = repmat(nn/norm(nn),2,1);
% %          nn1 = theNormals{1}(s,:);
% 
%          dp = p2-p3;
%          nn = [-dp(2) dp(1)];
%          nn2 = repmat(nn/norm(nn),2,1);
%          nn = [nn1;nn2];
%          % flux jump
%          normalq = sum(dq.*nn,2);
%          % distance elem 1
%          di1 = Ne*Xe1 - repmat(p1,2,1);
%          di1 = sqrt(sum(di1.^2,2));
%          % distance elem 2
%          di2 = Ne*Xe2 - repmat(p2,2,1);
%          di2 = sqrt(sum(di2.^2,2));
%          di = [di1; di2+max(di1(:))];
%          % normalFlux * shape functions 
%          qN = normalq .* [Ne(2,:)' ; Ne(1,:)'];
%          % integrate
%          I1 = sum(qN(1:2).*wi) * det(Ne_xi*di1);
%          I2 = sum(qN(3:4).*wi) * det(Ne_xi*di2);
%          Ie = I1+I2;
%          % output
%          fprintf('dof %i : \\int N*Aq*n = %+g \n',s,Ie)
%       end   
%    end
% end
