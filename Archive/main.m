addpath('../XFEM/');
addpath('../XFEM+/');
%
% MAIN FILE FOR COMPUTATION OF THE "GAS NATURAL" CASE
%
%close all
clear all
clc
% profile on;
% GLOBAL VARIABLES
% cond: permeability value [ wet_soil water ]
global cond useEnrichment debugData lambda
cond = [1 1000];
% useEnrichment:
%    0=FEM
%    1=XFEM
%    2=XFE+hat
%    3=XFE+cordero
%    4=XFE+least squares (not tested)
%    5=XFE+ 2 semi hat
%    6=XFE+ inverse hat
useEnrichment = 5;
% lambda only used in the cordero version!
%lambda = eye(2); % Variar los puntos sobre el Level Set: Otras opciones:
%lambda = [1 3; 3 1]/4;
lambda= [1 0; .5 .5];

% MESH: VARIABLES X AND T
%  meshI.mat (I=1:4): uniform rectangular mesh
%  meshDiscI.mat (I=1:4): uniform circular mesh
%  meshHole1.mat: rectangular mesh with a circular hole in the center
%  meshGasNatural1.mat. mesh for the "Gas Natural" case
% load 'meshDisc2.mat'
% load 'meshsmall.mat'
load 'meshSquare02.mat'
%filename = 'error_FE_Mesh4';

% BOUNDARY CONDITIONS
%  1: test cases with h=0 on all sides (a rectangular mesh)
%  2: test cases with h=0 below and q=0 on other sides (rectangular mesh)
%  3: test cases with h=0 on left, h==1 on right and q=0 on other sides
%  4: test cases with h=0 below, h==1 on top and q=0 on other sides
%     (rectangular mesh)
%  5: (2) and fixed h=1 at node closest to [.68 -.68]
%  6: test cases with h=0 on right, h==1 on left for points higher than
%     .55, h==2 on left for points lower than .55, and q=0 on other sides
%     (rectangular mesh)
%  7: (2) and fixed h=1 at node closest to [.68 -.68]
%  8: (2) and fixed h=1 at node closest to [.86 -.67]
%  9: q=0 up and down, h=0 on the left, h=1 on the right
% 10: u=0 around a circle (radius 1)
% 11: q=0 up and down, h=1 on the left and right
% caseBoundaryCondition = 10;
caseBoundaryCondition = 4;

% LOADING CONDITIONS
%  1: uniform source term f=1
%  2: uniform source term f=0
%  3: source term f=1 within a circle ([.5 .5],.1) and f=0 elsewhere.
%  4: uniform source term f=-4
% caseLoad = 4;
caseLoad = 2;

% INITIAL LEVEL SET FUNCTION
%  1: no level set
%  2: horizontal split at y=.55
%  3: "almost" horizontal split around y=.55
%  4: Ellipsoidal Level Set in direction [-1 1] from point [.75 -.75]
%  5: Corner split at [.48 .52]
%  6: Ellipsoidal Level Set in direction [-1 1] from point [.5 .5]
%  7: horizontal split at y=-.6
%  8: Ellipsoidal Level Set with angle 157/180*pi from point [.672 -.688]
%  9: same as 8 but for when the water table reaches the jet
%  10/11: Ellipsoid intersecting the exterior boundary (for tests)
%  12: Ellipsoidal Level Set with angle 157/180*pi from point [.86 -.67]
%  13: tube at y in [-.2 +.2]
%  14: same as 13 but with sine walls
%  15: horizontal line at y=.55 with an indent in the middle
%  16: horizontal line at y=.5 with an indent in the middle
%  110: line split around y=.55, with an angle of 10 degrees
%  17: circle of radius .25
%  18: circle of radius .5
%  19: circle of radius .75
% initialLevelSet = 3;
% initialLevelSet = 18;
initialLevelSet = 113;

% TIME DISCRETIZATION
time.initial = 0;
time.t = time.initial;
time.numberSteps = 1;

% 'REAL TIME' GRAPHICAL OUTPUT
% to draw outputs while the computation is performing
% Step: [integer] number of time steps between two drawings. Use Step=Inf
%       if you do not want real time graphics.
% Type: 'Front'. 'Flux', 'h', 'Gradient Norm q', or 'Level Set'. Watch out
%       that 'Flux' is not very stable
% Movie: [logical] to indicate whether to make a movie
% MovieFileName: to indicate the name of the file where the movie is stored
% GasNatural: logical equal to 1 to enforce special formats for the Gas
%       Natural case
graphiks.Step = 10;
graphiks.Type = 'Flux';
graphiks.Movie = ( graphiks.Step ~= Inf );
graphiks.MovieFileName = 'moving_levelset.avi';
graphiks.GasNatural = 0;
graphiks.Paper = 0;
graphiks.TimeStep = 10800;

% Update: number of time step between total updating of the level set
%         as the distance to phi=0. Choose Inf to desactivate this
%         function
Update = Inf;

% GENERAL OPTIONS
opts.tolerance = calcMinEdgeLength( X, T )* 0.1;
% Level set tolerance: to detect if a node is on the interface
opts.cutLevelSet = 1;
opts.numberOfGaussPoints = 4;
opts.elementType = 2; % 2=triangle
opts.numberOfElementNodes = 3;
opts.TimeDisplayStep = 1;
opts.BoundaryCondition = caseBoundaryCondition;
opts.LoadCondition = caseLoad;

% QUADRATURE
[pospg,pespg] = quadrature( opts.elementType, opts.numberOfGaussPoints );
[N,Nxi,Neta] = shapeFunction( opts.elementType, ...
   opts.numberOfElementNodes, pospg );
mass = MatrixM( X, T, pospg, pespg, N, Nxi, Neta );
CarLength = .03*ones(size(X,1),1);

% INITIALIZATION
levelSet = predefinedLevelSet( X, initialLevelSet, opts );
%%%%%%
warning('levelSet(<tol) = 0')
levelSet(abs(levelSet)<=opts.tolerance) = 0;
%%%%%%
dt = 1e-10;
time.step = dt;
h = [];
hE = [];
out = struct('t', zeros( 1, time.numberSteps ), ...
   'levelSet', zeros( size(X,1), time.numberSteps ) );
% LOOP ON TIME
tic


for i1=1:time.numberSteps
   
   if floor(i1/opts.TimeDisplayStep)==(i1/opts.TimeDisplayStep)
      disp(['Time step ' num2str(i1) ' of ' num2str(time.numberSteps)])
      toc,tic
   end
   
   % DARCY MODEL WHEN THE WATER TABLE IS HIGH
   [ h, hE, Ktot, R ] = PoissonXFEM( X, T, pospg, pespg, N, Nxi, Neta, ...
      levelSet, opts );
   
   %         test_movie_1 = graphiks.Step ~= Inf;
   %         test_movie_2 = floor(i1/graphiks.Step) == i1/graphiks.Step;
   %         if test_movie_1 && test_movie_2
   %             RealTimeGraphics( 100+i1, graphiks, i1, X, T, levelSet, ...
   %                             Inf, h, hE, pospg, N, Nxi, Neta, opts );
   %         end
   %         break
   Vn = NormalVelocityInFlow( X, T, pospg, N, Nxi, Neta, levelSet, ...
      h, hE,   opts);
   

   [ Xpg, qH ] = Flujos( h, hE, X, T, pospg, N, Nxi, Neta, levelSet, opts );   

%    plotStuff
%    pause
%    plotFrame
%    str = sprintf('xfem+-frame%3i.pdf',i1);
%    saveas(gcf,str)
%    pause
   
   
   % RESOLUTION OF HAMILTON-JACOBI + REDISTANCING
   if i1 ~= time.numberSteps
      [ levelSet, dt ] = EvolutionLevelSet( X, T, levelSet, Vn, ...
         N, Nxi, Neta, CarLength, opts );
   end
   
   % SAVING FOR POST-TREATMENT
   out.t( i1 ) = time.t;
   out.levelSet ( :, i1 ) = levelSet;
   
   % UPDATING TIME
   time.t = time.t + dt;
   time.step = dt;
   
end


% % compute energies for the error
% f = 4;
% r0 = .5;
% [ eL2, eH1, eH121, ecomp1, eH122, ecomp2, exL2,exH1 ] = ...
%                        ComputeEnergy( X, T, h, hE, levelSet, r0, f, opts );
% [ sum(eL2) sum(exL2) sum(eH1) sum(exH1)]
% [sum(eH121) sum(ecomp1)]
% [sum(eH122) sum(ecomp2)]
%
%eval(['save draft/files_test/' filename '.mat X T eL2 eH1 eH121 eH122 ' ...
%         'ecomp1 ecomp2 h hE levelSet opts exH1 exL2'])

% switch initialLevelSet
%     case 2
%         figure; fill( [0 0 1 1], [0 0.55 0.55 0], 0.7*[1 1 1] )
%         hold on; fill( [0 0 1 1], [0.55 1 1 0.55], [1 1 1] )
%     case 3
%         figure; fill( [0 0 1 1], [0 .65 .45 0], 0.7*[1 1 1] )
%         hold on; fill( [0 0 1 1], [.65 1 1 .45], [1 1 1] )
%end

plotStuff

return
% energy computation
u = [ h; hE ];
Nu = length(u);
eH = u' * Ktot(1:Nu,1:Nu) * u

return

% verify that it is the same as above when nu1=nu2=1
% [ eL2, eH1 ] = ComputeEnergy( X, T, h, hE, levelSet, r0, 0, opts );
% sum(eH1)

return

% MOVIE
% MakeMovieGasNatural( out, X, T, time, graphiks );
% profile viewer

% energy computation
u = [ h; hE ];
Nu = length(u);
eH = u' * Ktot(1:Nu,1:Nu) * u /2
% error estimation: Eex = pi*f^2/8*(r0^4/nu1+(R^4-r0^4)/nu2)
eex = pi*f^2/8*(r0^4/cond(2)+(1-r0^4)/cond(1))
x0 = 0.5625;
nu1 = cond(1); nu2 = cond(2);
a = nu2 / (1-(1-nu2/nu1)*x0);
x = 0:0.001:1;
hex = zeros(size(x));
hex(x<x0) = a*x(x<x0)/nu1;
hex(x>=x0) = 1-a/nu2*(1-x(x>=x0));




