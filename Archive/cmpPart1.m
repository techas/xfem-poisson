%
% MAIN FILE FOR COMPUTATION OF THE "GAS NATURAL" CASE
%
addpath('../XFEM/');
addpath('../XFEM+/');
%close all
clear all
clc

global cond useEnrichment debugData lambda debug
cond = [1 1000];

switch 13
   case 1
      %% anda OK!
      load 'meshDisc2.mat'
      caseBoundaryCondition = 10;
      caseLoad = 4;
      initialLevelSet = 18;
      opts.tolerance = calcMinEdgeLength( X, T ) * 0.05;

   case 2 
      %% caso para mostrar triangulos
      % los dos flujos se ven bien en las flechas
      % la comparativa muestra que xfem+ es mejor
      load 'meshDisc1.mat'
      caseBoundaryCondition = 10;
      caseLoad = 4;
      initialLevelSet = 18;
      opts.tolerance = calcMinEdgeLength( X, T ) * 0.1;
  
   case 3
      load 'meshDisc1.mat'
      caseBoundaryCondition = 10;
      caseLoad = 4;
      initialLevelSet = 18;
      opts.tolerance = calcMinEdgeLength( X, T ) * 0.345;

   case 4
      load 'meshDisc2.mat'
      caseBoundaryCondition = 10;
      caseLoad = 4;
      initialLevelSet = 18;
      opts.tolerance = calcMinEdgeLength( X, T ) * 0.1;
   case 5
      % caso que anda muy mal por Neumann BC
      load 'meshSquare03.mat'
      caseBoundaryCondition = 4;
      caseLoad = 2;
      initialLevelSet = 3;
      opts.tolerance = calcMinEdgeLength( X, T ) * 0.05;
   case 6
      % caso donde hay un elemento malo NO enriquecido 
      load 'meshSquare03.mat'
      caseBoundaryCondition = 4;
      caseLoad = 2;
      initialLevelSet = 112;
      opts.tolerance = calcMinEdgeLength( X, T ) * 0.0;
   case 7
      % caso que muestra el caso bueno y el malo con la configuracion de
      % los elementos partidos en dos.
      load 'meshSquare02.mat'
      caseBoundaryCondition = 4;
      caseLoad = 2;
      initialLevelSet = 110;
      opts.tolerance = calcMinEdgeLength( X, T ) * 0.05;
      
   case 10 
      % Caso malo
      % level set horizontal + malla estructurada
      load 'meshSquare02.mat'
      % se arregla con malla no estructurada. p.e.
%       load 'meshHet.mat'
      caseBoundaryCondition = 4;
      caseLoad = 2;
      % y se aregla tambien variando el levelSet. e.g. 114
      initialLevelSet = 115;
      opts.tolerance = calcMinEdgeLength( X, T ) * 0.05;
   
   %% FIGURAS PARA EL PAPER   
   case 11 
      % figure 3
      load 'meshSquare02.mat'
      caseBoundaryCondition = 4;
      caseLoad = 2;
      initialLevelSet = 113;
      opts.tolerance = calcMinEdgeLength( X, T ) * 0.05;
      
   case 12
      % figure 4 & 5 & 6
      load 'meshSquare02.mat'
      caseBoundaryCondition = 4;
      caseLoad = 2;
      initialLevelSet = 116;
      opts.tolerance = calcMinEdgeLength( X, T ) * 0.05;

   case 13
      % figure 7 part 1
      load 'meshSquareUnstructured-2.mat'
      caseBoundaryCondition = 4;
      caseLoad = 2;
      initialLevelSet = 117;
      opts.tolerance = calcMinEdgeLength( X, T ) * 0.01;

   case 14
      % figure 8 Disk
      load 'meshDisc2.mat'
      caseBoundaryCondition = 10;
      caseLoad = 4;
      initialLevelSet = 18;
      opts.tolerance = calcMinEdgeLength( X, T ) * 0.01;
end



%% -- Solve
% useEnrichment:
%    0=FEM
useEnrichment = 0;
cmpPart2

% useEnrichment:
%    1=XFEM
useEnrichment = 1;
cmpPart2

%    5=XFE+ 2 semi hat
useEnrichment = 5;
cmpPart2

%% -- Postproc
cmpPart3


