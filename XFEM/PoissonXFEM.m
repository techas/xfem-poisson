function [h,hE,Ktot,R] = PoissonXFEM( X, T, pospg, pespg, N, Nxi, Neta, ...
   levelSet, opts )
% POISSONXFEM to solve the Poisson equation using the XFE method
%
%  syntax [h,levelSet,out] = PoissonXFEM( X, T, pospg, pespg, N, Nxi,
%                                                  Neta, levelSet, opts );
%
%  X:            nodal coordinates
%  T:            connectivity matrix
%  pospg, pespg: Gauss points and weights
%  N, Nxi, Neta: FE interpolation functions and derivatives
%  levelSet:     vector giving the value at each node of the distance to
%                the phase interface. Values on one side are positive and
%                on the other negative
%  opts: structured array containing options

% R. Cottereau 03/2008 (based on code by S. Zlotnik and P. Diez)

global useEnrichment
Ne = size( X, 1 );
R = [];


% create FE matrix and right-hand side
[ K, f, enrichedNodes ] = createMatrix( X, T, levelSet, pospg, pespg, ...
   N, Nxi, Neta, opts.LoadCondition, opts.tolerance );
Nr = length( enrichedNodes );

% Boundary conditions
ksize = size( K, 1 );
[Abc,bbc] = makeBoundaryConditions( X, ksize, opts.BoundaryCondition );

if useEnrichment > 0
   if useEnrichment==2
      R = makeEnrichedValuesConditions( X, T, Nxi, Neta, levelSet, opts);
   end
   if useEnrichment==3
      R = makeEnrichedValuesConditionsCordero( X, T, N, Nxi, Neta, ...
         levelSet, opts);
   end
   if useEnrichment==4
      R = makeEnrichedValuesConditionsTest( X, T, levelSet, opts);
   end
   if useEnrichment==5
      R = makeEnrichedValuesConditionsTest3( X, T, Nxi, Neta, levelSet, ...
         opts);
%       keyboard
%       R([ 1 2],:) = [];
   end
   if useEnrichment==6
      R = makeEnrichedValuesConditionsTest4( X, T, Nxi, Neta, levelSet, ...
         opts);
   end
   %
   Abc = [ R ; Abc ];
   bbc = [ zeros(size(R,1),1); bbc ];
   %
   
%    if size(R,1) < 201
%       rankR = rank(R);
%       if rankR < size(R,1)
%          warning('R IS RANK DEFICIENT!!!!')
%          fprintf('rank(R)=%g\n#rows(R)=%g\n',rankR,size(R,1))
%       end
%    else
%       disp('rank not checked')
%    end

   [AbcE,bbcE] = makeEnrichedBoundaryConditions( X, T, levelSet, ...
      ksize, enrichedNodes, opts.BoundaryCondition );
   
   Z0 = zeros(size(AbcE,1),size(Abc,2)-size(AbcE,2));
   Abc = [Abc ; Z0 AbcE];
   bbc = [bbc ; bbcE];
end

% Mount and solve the system
Z0 = zeros( size( Abc, 1 ), size( Abc, 1 ) );
Ktot = [K    Abc'
        Abc  Z0];
ftot = [f; bbc];
aux = Ktot \ ftot;

% output
h  = aux( 1:Ne );
if ~useEnrichment
   Nr = 0;
end
hE = aux( Ne + (1:Nr) );
% multip = aux( (Ne+Nr+1):end );

%max(abs(R(1:81,:)*[h;hE]))
%keyboard