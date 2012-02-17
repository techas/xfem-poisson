function h = plotMesh( X, T, options )
% h = plotMesh( X, T, options )
%
% plot a mesh of triangular elements
%
% INPUT
%   X,T      a mesh
%   options 
%            .plotNodes
%            .lineSpecNodes
%            .labelNodes (0=off, 1=on)
%            .labelNodesFontSize
%            .labelNodesColor
%            .plotElements
%            .lineSpecElements
%            .labelElements (0=off, 1=on)
%            .labelElementsFontSize
%            .labelElementsColor
%
% OUTPUT
%   h        handler to created objects
%
if nargin == 2
   options = [];
end
if ~isfield( options, 'plotNodes' )
   options.plotNodes = 0;
end
if ~isfield( options, 'lineSpecNodes' )
   options.lineSpecNodes = '.k';
end
if ~isfield( options, 'labelNodes' )
   options.labelNodes = 0;
end
if ~isfield( options, 'labelNodesFontSize' )
   options.labelNodesFontSize = 8;
end
if ~isfield( options, 'labelNodesColor' )
   options.labelNodesColor = 'k';
end
if ~isfield( options, 'plotElements' )
   options.plotElements = 1;
end
if ~isfield( options, 'lineSpecElements' )
   options.lineSpecElements = '-b';
end
if ~isfield( options, 'labelElements' )
   options.labelElements = 0;
end
if ~isfield( options, 'labelElementsFontSize' )
   options.labelElementsFontSize = 8;
end
if ~isfield( options, 'labelElementsColor' )
   options.labelElementsColor = 'b';
end
if ~isfield( options, 'addLocalNumbering' )
   options.addLocalNumbering = 0;
end
hnodes = [];
helems  = [];
order = [1  2  3  1];          
hold on
if options.plotNodes && ~options.labelNodes
   hnodes = plot( X(:,1), X(:,2), options.lineSpecNodes, ...
      'Tag', 'nodeLineTag' );
end
if options.plotNodes && options.labelNodes
   hnodes = zeros( 1, length( X ) );
   for I = 1:length( X )
       hnodes(I) = text( X(I,1), X(I,2), int2str( I ), ...
          'FontSize', options.labelNodesFontSize, ...
          'Color', options.labelNodesColor, ...
          'Tag', 'nodeNumberTag' );
   end
end
if options.plotElements
   helems = zeros( 1, length( T ) );
   for j = 1:size( T, 1 )
      if options.labelElements
         prom = mean( X( T(j,1:3),:) );
         helems(j) = text( prom(1), prom(2), int2str( j ), ...
            'FontSize', options.labelElementsFontSize, ...
            'Color', options.labelElementsColor, ...
            'Tag', 'elementNumberTag' );
      else
         helems(j) = plot( X(T(j,order),1), X(T(j,order),2), ...
            options.lineSpecElements, ...
            'lineWidth', 0.5, ...
            'Tag', 'elementLineTag' );
      end
      if options.addLocalNumbering
         prom = mean( X( T(j,1:3),:) );
         pos =  X( T(j,1:3),:) + repmat(prom,3,1); 
         for I=1:3
             text( pos(:,1), pos(:,2), '%i', I, ...
               'FontSize', options.labelElementsFontSize-1, ...
               'Color', 'k' );
         end
      end
   end
end
hold off
h = [hnodes helems];
