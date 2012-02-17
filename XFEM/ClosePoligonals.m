function [polis, Ei] = ClosePoligonals( polis, Ei, X, levelSet );
% CLOSEPOLIGONALS to close the poligonals that were created to identify the
% level set function.
%
%  syntax: polis = ClosePoligonals( polis, X, levelSet );
%
%  polis:    list of poligonals [Np*1 cell array of np*2 matrices]. Note
%            that the np may be different in each cell
%  X:        nodal coordinates [Nn*d matrix]
%  levelset: nodal values of the level set function [Nn*1 vector]
%
%  polis:    (in output) list of closed poligonals [NP*1 cell array of 
%            nP*2 matrices]. Note that the sum of np on Np should be equal
%            to the sum of NP on nP
%
%  SEE ALSO: MAKEPOLIGONALSFROMSEGMENTS

% R. Cottereau 04/2008

% close poligonals that do not touch the exterior boundary (by closing them
% with the closest point ... this is probably not very stable in
% complicated situations (!!!)
[polis, Ei] = ClosePoligonalInternal( X, polis, Ei );

% close poligonals that touch the exterior boundary
if length(polis) > 1
    polis = ClosePoligonalMultipleBoundary( X, polis, levelSet );
end

%==========================================================================
function [polis2, Ei2] = ClosePoligonalInternal( X, polis, Ei );
% this function searches the ends of the poligonals contained in polis,
% selects those that are not on the exterior boundary and connects each of
% them to the closest ends

% constants
n = length(polis);

% get the ends in all the poligonals
[ Ends, isBnd ]= getEnds( polis, X );

% SPECIAL CASE FOR GAS NATURAL !!!!!
if Ends(1,1)<.6
    Ends( [1 2], : ) = Ends( [2 1], : );
    isBnd( [1 2], 1 ) = isBnd( [2 1], 1 );
    polis{1} = polis{1}( end:-1:1, : );
    Ei{1} = Ei{1}( end:-1:1, 1 );
end

% for all internal nodes, find the closest internal node
internal = find( isBnd == 0 );
Ni = length( internal ) / 2;
pairs = zeros( Ni, 2 );
for i1 = 1:Ni
    P1 = Ends( internal( 1 ), : );
    P2 = Ends( internal( 2:end ), : );
    dist = ( P2( :, 1 ) - P1(1) ).^2 + ( P2( :, 2 ) - P1(2) ).^2;
    ind = find( dist == min( dist ) );
    pairs( i1, : ) = [ internal( 1 ) internal( 1+ind ) ];
    ind = setdiff( 2:2*(Ni-i1+1), 1+ind );
    internal = internal( ind );
end

% create new poligonal cell
n2 = n - ceil( Ni/2 );
polis2 = cell( n2, 1 );
Ei2 = cell( n2, 1 );

% do a loop on all poligonals
% NOTE THAT IN THE CURRENT IMPLEMENTATION IT IS NOT POSSIBLE TO HAVE
% MULTIPLY CONNECTED POLIGONALS, THERE MAY BE ONLY ONE 'HOLE' ...
WasTreated = zeros( n, 1 );
ind1 = 1;
for i1 = 1:n
    
    % check that the poligonal was not considered
    if WasTreated( i1 ) == 0
        
        % start the new poligonal
        polis2{ ind1 } = polis{ i1 };
        Ei2{ ind1 } = Ei{ i1 };
    
        % look for internal ends at the end of the poligonal and connect
        % them to the next poligonal (read in 'pairs')
        if isBnd( i1*2 ) == 0
            
            [ i2, j2 ] = find( pairs == i1*2 );
            i2 = pairs( i2, setdiff( [1 2], j2 ) );
            nextPoli = ceil( i2 / 2 );
            
            % when the connecting point is at the end of the next poligonal
            % change the ordering before connecting
            if nextPoli ~= i1
                if nextPoli == i2/2
                    polis{ nextPoli } = polis{ nextPoli }( end:-1:1, : );
                    Ei{ nextPoli } = Ei{ nextPoli }( end:-1:1 );
                end
                polis2{ ind1 } = [ polis2{ ind1 } ; polis{ nextPoli } ];
                Ei2{ ind1 } = [ Ei2{ ind1 } ; 0; Ei{ nextPoli } ]; 
            end
            
            % indicating that the next poligonal was already treated
            WasTreated( nextPoli ) = i1;
            
        end

        % look for internal ends at the beginning of the poligonal
        % change the ordering before connecting
        if isBnd( (i1-1)*2+1 ) == 0
            
            polis2{ ind1 } = polis2{ ind1 }( end:-1:1, : );
            Ei2{ ind1 } = Ei2{ ind1 }( end:-1:1, : );
            [ i2, j2 ] = find( pairs == (i1-1)*2+1 );
            i2 = pairs( i2, setdiff( [1 2], j2 ) );
            nextPoli = ceil( i2 / 2 );
            
            % when the connecting point is at the end of the next poligonal
            % change the ordering before connecting
            if ( nextPoli ~= i1 ) & ( WasTreated( nextPoli ) ~= i1 )
                if nextPoli == i2/2
                    polis{ nextPoli } = polis{ nextPoli }( end:-1:1, : );
                    Ei{ nextPoli } = Ei{ nextPoli }( end:-1:1 );
                end
                polis2{ ind1 } = [ polis2{ ind1 } ; polis{ nextPoli } ];
                Ei2{ ind1 } = [ Ei2{ ind1 } ; 0; Ei{ nextPoli } ]; 
            end
            
            % indicating that the next poligonal was already treated
            WasTreated( nextPoli ) = i1;
            
        end
        
        % increment the counter for polis2
        ind1 = ind1+1;
            
    end
    
end

%==========================================================================
function polis = ClosePoligonalBoundary( X, polis, levelSet );
% use to add the necessary points to close the poligonal when the first and
% last nodes are not the same
% if the poligonal reaches the exterior boundary, close by a large loop
% (check the current level set values to choose whether the loop goes up or
% down. (ONLY WORKS FOR RECTANGULAR EXTERIOR BOUNDARY !!)
% if the poligonal does not reach the exterior boundary, close it by
% simply repeating the first node at the end

x_lo = min( X(:,1) );
x_up = max( X(:,1) );
y_lo = min( X(:,2) );
y_up = max( X(:,2) );
margen = 2 * max(abs(x_lo) + abs(x_up) , abs(y_lo) + abs(y_up) );

for i1=1:length(polis)
    poli1 = polis{i1};
    if sum(poli1(1,:)==poli1(end,:),2)~=2
        if (abs(poli1(end,1)-x_lo) < 1e-8)
           % the last node of the poligonal is on the left side
           SuppPoli = ClosePoliFromOnePoint(X(:,1),X(:,2), ...
                x_lo,y_lo,x_up,y_up, ...
                poli1(end,1),poli1(end,2),poli1(1,1),poli1(1,2), ...
                levelSet,margen);
        else
            if (abs(poli1(end,2)-y_up) < 1e-8)
           % the last node of the poligonal is on the upper side
           % change x by -y and y by x
               SuppPoli = ClosePoliFromOnePoint(-X(:,2),X(:,1), ...
                    -y_up,x_lo,-y_lo,x_up, ...
                    -poli1(end,2),poli1(end,1),-poli1(1,2),poli1(1,1), ...
                    levelSet,margen);
               SuppPoli = [SuppPoli(:,2) -SuppPoli(:,1)];
            else
                if (abs(poli1(end,1)-x_up) < 1e-8)
           % the last node of the poligonal is on the right side
           % change x by -x
                  SuppPoli = ClosePoliFromOnePoint(-X(:,1),X(:,2), ...
                        -x_up,y_lo,-x_lo,y_up, ...
                        -poli1(end,1),poli1(end,2),-poli1(1,1),poli1(1,2), ...
                        levelSet,margen);
                    SuppPoli = [-SuppPoli(:,1) SuppPoli(:,2)];
                else
                    if (abs(poli1(end,2)-y_lo) < 1e-8)
          % the last node of the poligonal is on the lower side
          % change x by y and y by -x
                      SuppPoli = ClosePoliFromOnePoint(X(:,2),-X(:,1), ...
                            y_lo,-x_up,y_up,-x_lo, ...
                            poli1(end,2),-poli1(end,1),poli1(1,2),-poli1(1,1), ...
                            levelSet,margen);
                        SuppPoli = [-SuppPoli(:,2) SuppPoli(:,1)];
                    else
          % internal node: simply close the loop by repeating the
          % first node
                        SuppPoli = poli1(1,:);
                    end
                end

            end
        end
        poli1 = [poli1; SuppPoli];
    end
    polis{i1} = poli1;
end
%==========================================================================
function polis = ClosePoligonalMultipleBoundary( X, polis, levelSet );
% use to add the necessary points to close the poligonal when the first and
% last nodes are not the same
% if the poligonal reaches the exterior boundary, close by a large loop
% (check the current level set values to choose whether the loop goes up or
% down. (ONLY WORKS FOR RECTANGULAR EXTERIOR BOUNDARY !!)
% if the poligonal does not reach the exterior boundary, close it by
% simply repeating the first node at the end

x_lo = min( X(:,1) );
x_up = max( X(:,1) );
y_lo = min( X(:,2) );
y_up = max( X(:,2) );
margen = 2 * max(abs(x_lo) + abs(x_up) , abs(y_lo) + abs(y_up) );

nPolis = length(polis);
xEnds = zeros(2*nPolis,2);
for i1=1:length(polis)
    xEnds((i1-1)*2+[1:2],:) = polis{i1}([1 end],:);
end
IsBnd = IsOnBoundary( X, xEnds );

% Supposition that only one level set exists: I end up with only one
% poligonal and all the parts connect to one another
% I also suppose that only one cut can take places on each side
poliTot = [];
iPoli = 1;
poli1 = polis{ iPoli };
for i1=1:nPolis
    IndNext = find( IsBnd == IsBnd(2*iPoli) );
    IndNext = setdiff( IndNext, 2*iPoli+[-1 0] );
    jPoli = ceil(IndNext/2);
    poli2 = polis{jPoli};
    if floor(IndNext/2)==IndNext/2
        poli2 = poli2(end:-1:1,:);
    end
    xf = poli1(end,1);
    yf = poli1(end,2);
    x0 = poli2(1,1);
    y0 = poli2(1,2);
    Bndf = IsOnBoundary(X, [xf yf]);
    if (x0~=xf | y0~=yf)
        if Bndf==1
           % the last node of the poligonal is on the left side
           SuppPoli = ClosePoliFromOnePoint(X(:,1),X(:,2), ...
                x_lo,y_lo,x_up,y_up, xf, yf, x0, y0, ...
                levelSet,margen);
        else
            if Bndf==2
           % the last node of the poligonal is on the upper side
           % change x by -y and y by x
               SuppPoli = ClosePoliFromOnePoint(-X(:,2),X(:,1), ...
                    -y_up,x_lo,-y_lo,x_up, -yf, xf, -y0, x0, ...
                    levelSet,margen);
               SuppPoli = [SuppPoli(:,2) -SuppPoli(:,1)];
            else
                if Bndf==3
           % the last node of the poligonal is on the right side
           % change x by -x
                  SuppPoli = ClosePoliFromOnePoint(-X(:,1),X(:,2), ...
                        -x_up,y_lo,-x_lo,y_up,  -xf, yf, -x0, y0, ...
                        levelSet,margen);
                    SuppPoli = [-SuppPoli(:,1) SuppPoli(:,2)];
                else
                    if Bndf==4
          % the last node of the poligonal is on the lower side
          % change x by y and y by -x
                      SuppPoli = ClosePoliFromOnePoint(X(:,2),-X(:,1), ...
                            y_lo,-x_up,y_up,-x_lo, yf, -xf, y0, -x0, ...
                            levelSet,margen);
                        SuppPoli = [-SuppPoli(:,2) SuppPoli(:,1)];
                    else
          % internal node: simply close the loop by repeating the
          % first node
                        SuppPoli = [];
                    end
                end

            end
        end
        poliTot = [ poliTot; poli1; SuppPoli ];
%        poliTot = [poli1; SuppPoli];
        iPoli = jPoli;
        poli1 = poli2;
    end
end
polis = {};
polis{1} = [ poliTot ; poliTot(1,:) ];

%==========================================================================
function SuppPoli = ClosePoliFromOnePoint(  x, y, x_lo, y_lo, x_up, y_up, ...
                                       x0, y0, xf, yf, levelSet, margen )

SuppPoli = [x_lo-margen y0];
BoundaryNodes = find(x==x_lo);
BoundaryNodes = BoundaryNodes(find(y(BoundaryNodes)> y0));
BoundaryNode = BoundaryNodes(find(y(BoundaryNodes)==min(y(BoundaryNodes))));
% if negative go up, else go down
UpDown = -sign(levelSet(BoundaryNode));

% the first node of the poligonal is on the same side as the last
if ((abs(xf-x_lo)<1e-8) & ...
                    ((UpDown==+1 & (yf > y0)) | (UpDown==-1 & (yf < y0))))
    SuppPoli = [SuppPoli
                SuppPoli(end,1) yf];
%    SuppPoli = [SuppPoli
%                SuppPoli(end,1) yf
%                xf yf];
else
    % go in the direction of the positive level set
    SuppPoli = [SuppPoli;
                SuppPoli(end,1) SuppPoli(end,2)+UpDown*margen];
    if (UpDown==1 & (abs(yf-y_up)<1e-8)) | ...
                                        (UpDown==-1 & (abs(yf-y_lo)<1e-8))
        % the node is on the following face
        SuppPoli = [SuppPoli
                    xf SuppPoli(end,2)];
%        SuppPoli = [SuppPoli
%                    xf SuppPoli(end,2)
%                    xf yf];
    else
        SuppPoli = [SuppPoli
                    SuppPoli(end,1)+2*margen SuppPoli(end,2)];
        if (abs(xf-x_up)<1e-8)
            % the node is on the right face
            SuppPoli = [SuppPoli
                        SuppPoli(end,1) yf];
%            SuppPoli = [SuppPoli
%                        SuppPoli(end,1) yf
%                        xf yf];
        else
            SuppPoli = [SuppPoli
                        SuppPoli(end,1) SuppPoli(end,2)-UpDown*2*margen];
            if (UpDown==1 & (abs(yf-y_lo)<1e-8)) | ...
                                        (UpDown==-1 & (abs(yf-y_up)<1e-8))
                % the node is on the following face
                SuppPoli = [SuppPoli
                            xf SuppPoli(end,2)];
%                SuppPoli = [SuppPoli
%                            xf SuppPoli(end,2)
%                            xf yf];
            else
                % the node is on the left face (where we started)
                SuppPoli = [SuppPoli
                            SuppPoli(end,1)-2*margen SuppPoli(end,2)
                            SuppPoli(end,1)-2*margen yf];
%                SuppPoli = [SuppPoli
%                            SuppPoli(end,1)-2*margen SuppPoli(end,2)
%                            SuppPoli(end,1)-2*margen yf
%                            xf yf];
            end
        end
    end
end
%==========================================================================
function [ Ends, isBnd ] = getEnds( polis, X );
% gets the ending elements in each of poligonals in polis and checks
% whether these ends are on the boundary of the domain
%
%  syntax: [ Ends, isBnd ] = getEnds( polis, X );
%
%  polis:    list of poligonals [Np*1 cell array of np*2 matrices]. Note
%            that the np may be different in each cell
%  X:        nodal coordinates [Nn*d matrix]
%  
%  Ends:  ending nodes in all the poligonals [(Np*2)*2 matrix]
%  isBnd: for each node in Ends, indicates whether it is on the exterior
%         boundary [(Np*2)*1 vector of integers] (see IsOnBoundary)

% R. Cottereau 04/2008

% constants
n = length( polis );

% initializations
Ends = zeros( n*2, 2 );
isBnd = zeros( n*2, 1 );

% loop on all poligonals to get the ends
for i1 = 1:n
    Ends( 2*i1-1 , : ) = polis{i1}( 1, : );
    Ends( 2*i1 , : ) = polis{i1}( end, : );
end

% check whether the ends are on the exterior boundary of the domain
isBnd = IsOnBoundary( X, Ends );
