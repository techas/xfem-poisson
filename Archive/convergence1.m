%
%
addpath('../XFEM/');
addpath('../XFEM+/');
%
% convergence plot
%
%close all
clc
global cond useEnrichment debugData lambda
cond = [1 1000];
lambda= [1 0; .5 .5];

%%
caseBoundaryCondition = 10; % u=0 around a circle (radius 1)
caseLoad = 4;               % uniform source term f=-4
initialLevelSet = 18;       % circle of radius .5
%%
levelSetTol = 0.01; % %
%%
opts.cutLevelSet = 1;
opts.numberOfGaussPoints = 4;
opts.elementType = 2; % 2=triangle
opts.numberOfElementNodes = 3;
opts.BoundaryCondition = caseBoundaryCondition;
opts.LoadCondition = caseLoad;
[pospg,pespg] = quadrature( opts.elementType, opts.numberOfGaussPoints );
[N,Nxi,Neta] = shapeFunction( opts.elementType, opts.numberOfElementNodes, pospg );

%    0=FEM
%    1=XFEM 
%    2=XFE+hat
%    3=XFE+cordero
%    4=XFE+least squares (not tested)
%    5=XFE+ 2 semi hat 
%    6=XFE+ inverse hat
% methods = [1 2 3 5];
methods = [0 1 5];
methodsLabel ={'FE','XFE','XFE+hat','XFE+cordero','XFE+LS','XFE+semi hat','XFE+inv hat'};

meshes = ['MeshDisc1.mat'
          'MeshDisc2.mat'
          'MeshDisc3.mat'];
% meshes = ['MeshDisc1.mat'
%           'MeshDisc2.mat'
%           'MeshDisc3.mat'
%           'MeshDisc4.mat'
%           'MeshDisc5.mat'];

% init
hnorm = zeros(length(methods),size(meshes,1));
maxRelErr = hnorm;
avrRelErr = hnorm;
hmin = hnorm;

for met = 1:length(methods)
   useEnrichment = methods(met);
   for mesh = 1:size(meshes,1)
      fprintf('\n\n--- scheme: %s,   mesh: %s ---\n',methodsLabel{met},meshes(mesh,:))
      debugData = {};
      
      load( meshes(mesh,:) )
      hh = calcMinEdgeLength( X, T );
      opts.tolerance = hh * levelSetTol;
      levelSet = predefinedLevelSet( X, initialLevelSet, opts );
      %
      warning('levelSet(<tol) = 0')
      levelSet(abs(levelSet)<=opts.tolerance) = 0;
      fprintf('#levelset nodes set to zero: %g\n', ...
         sum(abs(levelSet)<=opts.tolerance) );
      %
      [h,hE,Ktot,R] = PoissonXFEM(X,T,pospg,pespg,N,Nxi,Neta,levelSet,opts);
      
      %% energy computation
      u = [ h; hE ];
      Nu = length(u);
      eH = u' * Ktot(1:Nu,1:Nu) * u / 2;
      hnorm(met,mesh) = eH;
      
      %% Analytic solution
      r0 = 0.5;
      [pos,qun,qan] = circleNumericalAnalyticalFlux(r0,X,T,h,hE,levelSet,opts);
      
      re = abs(qun-qan) ./ abs(qan);
      maxRelErr(met,mesh) = max(re(:));
      avrRelErr(met,mesh) = mean(re(:));
      hmin(met,mesh) = hh;
      
      if met==1 && 1
         figure(50+mesh), clf
         plotMesh(X,T);
         hold on
         plot(pos(:,1),pos(:,2),'.r')
      end
      
   end
end

save('todo13.mat')

figure(44), clf, hold on
for met = 1:length(methods)
   x = hmin(met,:);
   y = maxRelErr(met,:);
   plot(log10(x),log10(y),[giveMeLineSpec(met) 'o'],'linewidth',2)
end
legend(methodsLabel{methods+1})
xlabel('log(h)')
ylabel('max err rel')
title('flux normal to the interface')


figure(45), clf, hold on
for met = 1:length(methods)
   x = hmin(met,:);
   y = avrRelErr(met,:);
   plot(log10(x),log10(y),[giveMeLineSpec(met) 'o'],'linewidth',2)
end
legend(methodsLabel{methods+1})
xlabel('log(h)')
ylabel('(log) mean err rel')
title('flux normal to the interface')


