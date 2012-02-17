function [outz,outw,outN,outNxi,outNeta] = ...
   shapeFunctionOnTheInterface(Xe,LSe,type,nPoints)
% [outz,outN,outNxi,outNeta] = ...
%    shapeFunctionOntheInterface(Xe,LSe,type,nPoints)
%
% returns
%
% linear elements assumed
%
% type = 1 -> level set crosses one node
% type = 2 -> level set does not cross any node
%
switch type
   case 1
   %% Case Ls crosses one node
   P1 = Xe(1,:);
   P2 = intersection( Xe([2 3],:), LSe([2 3]) );
   case 2
   %% Case Ls does not cross cross any node
   P1 = intersection( Xe([1 2],:), LSe([1 2]) );
   P2 = intersection( Xe([1 3],:), LSe([1 3]) );   
end

% 1D quadrature
line = 0;
nNodes = 2;
[zi,outw] = quadrature( line, nNodes );
% 1D shape functions
[Nl,Nxil] = shapeFunction( line, nPoints, zi );
% 1D integration points over the interface
outz = Nl*[P1; P2];
% map points to reference elements
[xi,eta] = invMap( Xe, outz(:,1), outz(:,2) );
triangle = 2;
nNodes = 3;
% evaluate shape function @ mapped points
[outN,outNxi,outNeta]= shapeFunction(triangle,nNodes,[xi,eta]);


if 0
   figure(88), clf, hold on
   ix = [1:3 1];
   plot(Xe(ix,1),Xe(ix,2),'-k')
   plot([P1(1) P2(1)],[P1(2) P2(2)],'-r')
   plot(outz(:,1),outz(:,2),'mx')   
   keyboard
end

