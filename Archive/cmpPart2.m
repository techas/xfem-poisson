%
% run model
%
debugData=[];
debug=[];


% lambda only used in the cordero version!
%lambda = eye(2); % Variar los puntos sobre el Level Set: Otras opciones:
%lambda = [1 3; 3 1]/4;
lambda= [1 0; .5 .5];


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
killedNodes = abs(levelSet) < opts.tolerance;
levelSet(killedNodes) = 0;
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
   
%    if floor(i1/opts.TimeDisplayStep)==(i1/opts.TimeDisplayStep)
%       disp(['Time step ' num2str(i1) ' of ' num2str(time.numberSteps)])
%       toc,tic
%    end
   
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
%    Vn = NormalVelocityInFlow( X, T, pospg, N, Nxi, Neta, levelSet, ...
%       h, hE,   opts);
   

   [ Xpg, qH ] = Flujos( h, hE, X, T, pospg, N, Nxi, Neta, levelSet, opts );   

%    plotStuff
%    pause
%    plotFrame
%    str = sprintf('xfem+-frame%3i.pdf',i1);
%    saveas(gcf,str)
%    pause
   
   
   % RESOLUTION OF HAMILTON-JACOBI + REDISTANCING
%    if i1 ~= time.numberSteps
%       [ levelSet, dt ] = EvolutionLevelSet( X, T, levelSet, Vn, ...
%          N, Nxi, Neta, CarLength, opts );
%    end
   
   % SAVING FOR POST-TREATMENT
%    out.t( i1 ) = time.t;
%    out.levelSet ( :, i1 ) = levelSet;
   
   % UPDATING TIME
%    time.t = time.t + dt;
%    time.step = dt;
   
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


% energy computation
u = [ h; hE ];
Nu = length(u);
format long
eH = u' * Ktot(1:Nu,1:Nu) * u


% plotStuff
if useEnrichment==0
   save('datafem.mat')
end
if useEnrichment==5
   save('dataxfemp.mat')
end
if useEnrichment==1
   save('dataxfem.mat')
end

return

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




