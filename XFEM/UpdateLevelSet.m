function levelSet = UpdateLevelSet( X, T, levelSet, it, phreaLevel, ...
                                                          Step, opts)

% R. Cottereau 03/2008

itt = it / Step;

if ( itt~=0 ) && ( floor(itt) == itt )
    
    type = classifyElements( levelSet, T, 0 );
    enrichedElements = find( type > 0 );

    % find the segments that are crossed by the level set
    [ Seg, SegsBnd ] = CrossedSegments( T, enrichedElements, levelSet );
                                                       
    % make poligonals from the segments crossed by the level set
    [ polis, Ei ] = MakePoligonalFromSegments( X, T, SegsBnd, Seg, ...
                                                              levelSet );

    % close the poligonals
    polis = ClosePoligonals( polis, Ei, X, levelSet );
    
    % create the level set
    LSsave = levelSet;
    levelSet = makeLS( X, polis, opts.cutLevelSet );

    % keep the levelSet over the phreatic level constant (= no update)
    overPhreatic = find( X(:,2) > phreaLevel );
    levelSet( overPhreatic ) = LSsave( overPhreatic );

    % possible reshape
    levelSet = reshape( levelSet, length( levelSet ), 1 );
    
end
