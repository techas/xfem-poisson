function ...
   [Abc,bbc] = makeEnrichedValuesConditionsCordero( X, T, N, Nxi, Neta, levelSet, opts)
% Condiciones de difusi?n
global cond lambda
nu1 = cond(1); nu2 = cond(2);
alfa = (nu1 + nu2)/(nu1 - nu2);
% Nodos enriquecidos
[type,enrichedNodes] = classifyElements( levelSet, T, opts.tolerance );
enrichedElements = find( type > 0 );
% N?mero de nodos
numberOfNodes = size(X,1);
% N?mero de nodos enriquecidos
numberOfenrichedNodes = length( enrichedNodes );
% N?mero de elementos enriquecidos
numberOfenrichedElements = length( enrichedElements );
% N?mero de ecuaciones del sistema general
numberOfEquations = numberOfNodes + numberOfenrichedNodes;
%
% Matriz de condiciones sobre los nodos enriquecidos (uno por elemento)
Abc = zeros(numberOfenrichedElements,numberOfEquations);
bbc = zeros(numberOfenrichedElements,1);
B = [ Nxi(1,:); Neta(1,:) ]; I = eye(3);
%
for ieleM = 1:numberOfenrichedElements
   ielem = enrichedElements(ieleM);
   Te = T(ielem,:); Xe = X(Te,:); LSe = levelSet(Te);
   
   RE = vectorFind( enrichedNodes, Te );
   FE = RE + numberOfNodes;
   Jacob = B*Xe ; %dJ = det(Jacob);
   G = Jacob\B;
   v = G'*G*LSe;
   [xi,eta] = LevelSetCero(LSe);
   Npc = shapeFunction( opts.elementType, 3, lambda*[xi,eta] );
   absLSe = abs(LSe);
   R = Npc*absLSe;
   g1 = absLSe - alfa*LSe;
   C = cell(1,2);
   for inpc = 1:size(Npc,1)
      C{inpc} = R(inpc)*I + g1*Npc(inpc,:);
   end
   d = [v'*C{1}; v'*C{2}];
   iv = 2*ieleM+(-1:0);
   Abc(iv,[Te FE]) = [repmat(v',2,1) d];
   bbc(iv,1) = 0;
end
%  [rank(Abc) size(Abc,1)]

