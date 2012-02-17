%
enrichedNodes = debugData.enrichedNodes;
enrichedElements = debugData.enrichedElements;

if isfield( debugData, 'polis' )
   polis = debugData.polis;
else
   [ Seg, SegsBnd ] = CrossedSegments( T, enrichedElements, levelSet, opts.tolerance );
   [ polis, dummy, Segi ] = MakePoligonalFromSegments( X, T, SegsBnd, Seg, levelSet );
end

if isfield( debugData, 'Ei' )
   Ei = debugData.Ei;
else
   [ Seg, SegsBnd ] = CrossedSegments( T, enrichedElements, levelSet, opts.tolerance );
   [ dummy, Ei, Segi ] = MakePoligonalFromSegments( X, T, SegsBnd, Seg, levelSet );
end

if isfield( debugData, 'theNormals' )
   theNormals = debugData.theNormals;
end
