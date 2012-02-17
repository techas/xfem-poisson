function LS = mkLS( X, polis, h )
% LS = mkLS( X, polis, h )
%
% Esta funcion arma un campo level set poniendo en cada nodo 
% la distancia con signo a las rectas descritan en polis.
%
% INPUT
%	X           Posiciones nodales.
%	polis{i}    cell de poligonos
%              poligono = Solo los nodos:
%                         El primero debe ser igual al ultimo
%                         Debe ser cerrado
%	h           Distancia maxima
%
% OUTPUT
%  LS          level set nodal
%
for i = 1:length( X )
   pto = X(i,:);
   if h == 0
      LS(i) = dpls( pto, polis );
   else
      LS(i) = min( dpls( pto, polis ), h );
   end
   if Adentrols( pto, polis )
      LS(i) = -LS(i);
   end
end
%------------------------------------------------------------------------
function res = dpls( pto, polis )
% distancia de un punto a un conjunto de poligonos
%  pto = [x y]
%  polis{i} un poligono
%
cp = size( polis, 2 );
for i = 1:cp
   dist(i) = dppoli( pto, polis{i} );
end
res = min( dist );

%------------------------------------------------------------------------
function res = dppoli( pto, poli )
% distancia de un punto a un poligono
% pto = [x y]
% poli = solo los nodos. se cierra el primero y el ultimo

n1 = size( poli, 1 ) - 1;
dist1 = zeros(n1,1);
for i = 1:n1
   rec = [ poli(i,1) poli(i,2) ; poli(i+1,1) poli(i+1,2) ];   
   dist1(i) = dpr( pto, rec );
end
res = min( dist1 );


%------------------------------------------------------------------------
function dist = dpr( p, r )
% dist punto recta
% punto = [x y]
% recta = [x1 y1 ; 
%          x2 y2 ]
p1  = r(1,:);
p2  = r(2,:);
dir = p2 - p1;
modd = norma(dir);
if modd ~= 0
   dir = dir / modd;
   pto = p - p1;
   f   = sum( dir.*pto );
   pr  = p1 + dir * f;
   if f < 0 | f > modd
      % proy fuera
      dist = min( dpp( p1, p ), dpp( p2, p ) );
   else
      %proy dentro
      dist = dpp( pr, p );
   end
else
   error( 'MkLS: un vector tiene long cero' )
end

%------------------------------------------------------------------------
function [pr,f,modd] = proy( x, r )
% proyeccion punto sobre recta
% punto = [x y]
% recta = [x1 y1 ; 
%          x2 y2 ]
p1  = r(1,:);
p2  = r(2,:);
dir = p2 - p1;
modd = norma( dir );
dir = dir / modd;
pto = x - p1;
f   = sum( dir.*pto );
pr  = p1 + dir * f;

%------------------------------------------------------------------------
function d = dpp( p1, p2 )
% dist punto punto
% punto = [x y]
d = ((p1(1)-p2(1)).^2 + (p1(2)-p2(2)).^2 ).^0.5;

function n = norma( v )
% norma
n = (v(:,1).^2 + v(:,2).^2).^0.5;


%------------------------------------------------------------------------
function nada
figure( 1 );
plot( punto(1), punto(2), 'bp', 'LineWidth', 2 );
p1 = [0 0];
for rad = 0 : pi/3 : 2*pi
   p2 = 2*[sin( rad ) cos( rad )];
   recta = [p1 ; p2];
   [pr,f,modd] = pro( punto, recta );
   disp( sprintf( 'p2: %f\t%f\t f: %f', p2(1), p2(2), f ) )
   
   
   figure( 1 );
   hold on;
   plot( recta(:,1), recta(:,2), 'r', 'LineWidth', 1 );
   text( p1(1), p1(2), 'P1', 'FontSize', 12, 'Color', 'r', ...
      'FontWeight', 'bold' )
   text( p2(1), p2(2), 'P2', 'FontSize', 12, 'Color', 'r', ...
      'FontWeight', 'bold' )
   
   if f < 0 | f > modd
      plot( pr(1), pr(2), 'kx', 'LineWidth', 2 );
   else   
      plot( pr(1), pr(2), 'gp', 'LineWidth', 2 );
   end
   box on; 
   axis equal tight;
   hold off;
end


%------------------------------------------------------------------------
function res = Adentrols( pto, polis )
% dice si el punto esta adentro de alg�n poligono
% pto = [x y]
% polis(:,:,i) = un poligono
cp = size( polis, 2 );
for i = 1:cp
   if Adentro( pto, polis{i} )
      res = 1;
      return
   end
end
res = 0;


%------------------------------------------------------------------------
function res = Adentro( v, poli )
% pto v = [x y]
% poli = solo los nodos. se cierra el primero y el ultimo
res = 0;
cantptos = size( poli, 1 );
if cantptos < 4
   return
end
pold = poli(cantptos,:);
for i = 1:cantptos
   pnew = poli(i,:);
   if pnew(1) > pold(1)
      P1 = pold;
      P2 = pnew;
   else
      P1 = pnew;
      P2 = pold;
   end
   if ( (pnew(1)<v(1)) == (v(1)<=pold(1)) ) & ...
         ( ((v(2)-P1(2))*(P2(1)-P1(1))) < ((P2(2)-P1(2))*(v(1)-P1(1))) )
      if res
         res = 0;
      else
         res = 1;
      end
   end
   pold = pnew;
end
