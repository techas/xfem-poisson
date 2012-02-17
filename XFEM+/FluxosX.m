function q = FluxosX( X, T, levelSet, h, hE, pos, lls, tol )
% compute the value of the fluxes in points pos
%
% syntax: q = FluxosX( X, T, levelSet, h, hE, pos, lls )
%
% X,T: coordinates and connectivity matrices
% levelSet: value of the value set function at the nodes in X
% h,hE: value of the solution at nodes, for the linear FE functions (h) and
%       the enriched functions (hE)
% pos: positions at which the fluxes should be evaluated
% lls: logical indicating which side of the levelset should be considered
%      (when evaluating the fluxes on the interface, lls=1 : LS>0, lls=0 : LS<0)

% R. Cottereau and S. Zlotnik 04/2011

% constants
global cond
q = zeros( size(pos) );

% is it standard fem?
isfem = length(hE) == 0;


if ~isfem
   % enriched elements
   [ type, enrichedNodes ] = classifyElements( levelSet, T, tol );
   enrichedElements = find( type > 0 );
   Nee = length( enrichedElements );
   
   triangle = 2;
   nPoints = 3;
   
   % loop on enriched elements
   for i1=1:Nee
      
      % local element and nodes
      Te = T( enrichedElements(i1), : );
      Xe = X( Te, : );
      LSe = levelSet( Te );
      %%
      LSe(abs(LSe)<tol) = 0;
      
      % calculo de N (linear FE functions) at pos
      lpos = inpolygon( pos(:,1), pos(:,2), Xe(:,1), Xe(:,2) );
      if sum(lpos) > 1
         Np = sum(lpos);
         [xi,eta] = invMap( Xe, pos(lpos,1), pos(lpos,2) );
         [Ne,Nxi,Neta] = shapeFunction( triangle, nPoints, [xi,eta] );
         B0 = [ Nxi(1,:); Neta(1,:) ];
         Be = ( B0 * Xe ) \ B0;
         % compute the ridge function
         Re = Ne*abs(LSe)-abs(Ne*LSe);
         if lls
            ind = LSe<0;
            gRe = -2 * Be(:,ind)*LSe(ind);
         else
            ind = LSe>0;
            gRe = 2 * Be(:,ind)*LSe(ind);
         end
         
         % value of the local flux
         ue = h( Te );
         ae = hE( vectorFind( enrichedNodes, Te ) );
         conde = cond( ~lls + 1 );
         q(lpos,:) = conde*( repmat((Be*ue),[1 Np]) + ((Be*ae)*Re') + (gRe*(Ne*ae)') )' ;
      end
   end
   
   
   
else
   %% flux of fem solution
   triangle = 2;
   nGaussPoints = 4;
   nElementNodes = 3;
   [pospg,pespg] = quadrature( triangle, nGaussPoints );
   [N,Nxi,Neta] = shapeFunction( triangle, nElementNodes, pospg );
   
   
   for i1=1:size(T,1)
      % local element and nodes
      Te = T( i1, : );
      Xe = X( Te, : );
      %     LSe = levelSet( Te );
      %%
      %     LSe(abs(LSe)<tol) = 0;
      
      % calculo de N (linear FE functions) at pos
      lpos = inpolygon( pos(:,1), pos(:,2), Xe(:,1), Xe(:,2) );
      if sum(lpos) > 1
         Np = sum(lpos);
         %        [xi,eta] = invMap( Xe, pos(lpos,1), pos(lpos,2) );
         %        [Ne,Nxi,Neta] = shapeFunction( triangle, nPoints, [xi,eta] );
         B0 = [ Nxi(1,:); Neta(1,:) ];
         Be = ( B0 * Xe ) \ B0;
         % compute the ridge function
         %        Re = Ne*abs(LSe)-abs(Ne*LSe);
         %        if lls
         %            ind = LSe<0;
         %            gRe = -2 * Be(:,ind)*LSe(ind);
         %        else
         %            ind = LSe>0;
         %            gRe = 2 * Be(:,ind)*LSe(ind);
         %        end
         
         % value of the local flux
         ue = h( Te );
         %        ae = hE( vectorFind( enrichedNodes, Te ) );
         conde = cond( ~lls + 1 );
         q(lpos,:) = conde * (repmat(Be*ue, [1 Np] ))';
      end
   end
end
