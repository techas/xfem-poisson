function levelSet = predefinedLevelSet( X, numls, opts )
% levelSet = predefinedLevelSet( X, h, lim, PARAM, ADV, numls )
%
% Create an initial level set field
%
% INPUT    
%   X          nodal coords
%   opts.cutLevelSet          min edge size
%   numls      number of leves set to create
%
% OUTPUT
%   levelSet   (column vector)
%
% numls
%    case 1      All ones
%    case 2      Horizontal split
%
x_lo = min( X(:,1) );
x_up = max( X(:,1) );
y_lo = min( X(:,2) );
y_up = max( X(:,2) );
%
switch numls
    case 1
        levelSet = ones( size( X, 1 ), 1 );

    case 2
        % Horizontal Split
        midy = 0.55;
        margen = 2*max( y_up, x_up );
        % Points
        poli1 = [x_lo-margen    midy ; ...
                 x_up+margen    midy ;...
                 x_up+margen    y_up+margen ;...
                 x_lo-margen    y_up+margen ;...
                 x_lo-margen    midy ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );

    case 3
        %
%         midy = 0.5;
%         margen = 2*max( y_up, x_up );
%         Points
%         poli1 = [x_lo-margen    midy+1.25 ; ...
%                  x_up+margen    midy-1.25 ;...
%                  x_up+margen    y_up+margen ;...
%                  x_lo-margen    y_up+margen ;...
%                  x_lo-margen    midy+1.25 ];
        midy = 0.41;
        margen = 2*max( y_up, x_up );
        % Points
        cc = 1.3;
        poli1 = [x_lo-margen    midy+cc ; ...
                 x_up+margen    midy-cc ; ...
                 x_up+margen    y_up+margen ;...
                 x_lo-margen    y_up+margen ;...
                 x_lo-margen    midy+cc ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );

    case 4
        fixedPointLS = [-5 0];
        axisEllipsLS = [1 0];
        numberNodesLS = 200;
        radiusLongLS = 10;
        radiusShortLS = .5;
%        fixedPointLS = [.69 -.69];
%        axisEllipsLS = [-1 1];
%        numberNodesLS = 200;
%        radiusLongLS = .03;
%        radiusShortLS = .01;
        
        axisEllipsLS = axisEllipsLS/norm(axisEllipsLS);
        thetaLS = 2*pi*[numberNodesLS-1:-1:0]'/numberNodesLS;
        nodeLS = [radiusLongLS*cos(thetaLS) radiusShortLS*sin(thetaLS)];
        nodeLS(:,1) = nodeLS(:,1) + radiusLongLS;
        if axisEllipsLS(1)==0;
            thetaEllipsLS = sign(axisEllipsLS(2))*pi/2;
        else
            thetaEllipsLS = atan(axisEllipsLS(2)/axisEllipsLS(1));
        end
        if axisEllipsLS(1)<0;
            thetaEllipsLS = thetaEllipsLS +pi;
        end
        thetaNodeLS = zeros(size(nodeLS,1),1);
        ind1 = find(nodeLS(:,1)==0);
        thetaNodeLS(ind1) = sign(nodeLS(ind1,1))*pi/2;
        ind1 = find(nodeLS(:,1)~=0);
        thetaNodeLS(ind1) = atan(nodeLS(ind1,2)./nodeLS(ind1,1));
        thetaNodeLS = thetaNodeLS + thetaEllipsLS;
        radiusNodeLS = sqrt(nodeLS(:,1).^2+nodeLS(:,2).^2);
        nodeLS = repmat(radiusNodeLS,1,2) .* ...
                                      [cos(thetaNodeLS) sin(thetaNodeLS)];
        nodeLS = nodeLS + repmat(fixedPointLS,size(nodeLS,1),1);
        
        polis{1} = nodeLS;
        levelSet = makeLS( X, polis, opts.cutLevelSet );

    case 5

        midy = 0.52;
        midx = 0.42;
        margen = 2*max( y_up, x_up );
        % Points
        poli1 = [x_lo-margen    midy ; ...
                 midx           midy ;...
                 midx           y_up+margen ;...
                 x_lo-margen    y_up+margen ;...
                 x_lo-margen    midy ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );

    case 6
        fixedPointLS = [.5 .5];
        axisEllipsLS = [-1 1];
        numberNodesLS = 20;
        radiusLongLS = .15;
        radiusShortLS = .1;
        
        axisEllipsLS = axisEllipsLS/norm(axisEllipsLS);
        thetaLS = 2*pi*[numberNodesLS-1:-1:0]'/numberNodesLS;
        nodeLS = [radiusLongLS*cos(thetaLS) radiusShortLS*sin(thetaLS)];
        nodeLS(:,1) = nodeLS(:,1) + radiusLongLS;
        if axisEllipsLS(1)==0;
            thetaEllipsLS = sign(axisEllipsLS(2))*pi/2;
        else
            thetaEllipsLS = atan(axisEllipsLS(2)/axisEllipsLS(1));
        end
        if axisEllipsLS(1)<0;
            thetaEllipsLS = thetaEllipsLS +pi;
        end
        thetaNodeLS = zeros(size(nodeLS,1),1);
        ind1 = find(nodeLS(:,1)==0);
        thetaNodeLS(ind1) = sign(nodeLS(ind1,1))*pi/2;
        ind1 = find(nodeLS(:,1)~=0);
        thetaNodeLS(ind1) = atan(nodeLS(ind1,2)./nodeLS(ind1,1));
        thetaNodeLS = thetaNodeLS + thetaEllipsLS;
        radiusNodeLS = sqrt(nodeLS(:,1).^2+nodeLS(:,2).^2);
        nodeLS = repmat(radiusNodeLS,1,2) .* ...
                                      [cos(thetaNodeLS) sin(thetaNodeLS)];
        nodeLS = nodeLS + repmat(fixedPointLS,size(nodeLS,1),1);
        
        polis{1} = nodeLS;
        levelSet = mkLS( X, polis, opts.cutLevelSet );
    case 7
        % Vertical Split
        midx = 0.5625;
        margen = 2*max( y_up, x_up );
        % Points
        poli1 = [midx           y_lo-margen ; ...
                 midx           y_up+margen ;...
                 x_up+margen    y_up+margen ;...
                 x_up+margen    y_lo-margen ;...
                 midx           y_lo-margen ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );
    case 8
        fixedPointLS = [.698 -.7525];
        anglLS = 150/180*pi;
%        fixedPointLS = [.7 -.7];
%        anglLS = 153/180*pi;
        axisEllipsLS = [cos(anglLS) sin(anglLS)];
        numberNodesLS = 100;
        radiusLongLS = .03;
        radiusShortLS = .02;
        
        axisEllipsLS = axisEllipsLS/norm(axisEllipsLS);
        thetaLS = 2*pi*[numberNodesLS-1:-1:0]'/numberNodesLS;
        nodeLS = [radiusLongLS*cos(thetaLS) radiusShortLS*sin(thetaLS)];
        nodeLS(:,1) = nodeLS(:,1) + radiusLongLS;
        if axisEllipsLS(1)==0;
            thetaEllipsLS = sign(axisEllipsLS(2))*pi/2;
        else
            thetaEllipsLS = atan(axisEllipsLS(2)/axisEllipsLS(1));
        end
        if axisEllipsLS(1)<0;
            thetaEllipsLS = thetaEllipsLS +pi;
        end
        thetaNodeLS = zeros(size(nodeLS,1),1);
        ind1 = find(nodeLS(:,1)==0);
        thetaNodeLS(ind1) = sign(nodeLS(ind1,1))*pi/2;
        ind1 = find(nodeLS(:,1)~=0);
        thetaNodeLS(ind1) = atan(nodeLS(ind1,2)./nodeLS(ind1,1));
        thetaNodeLS = thetaNodeLS + thetaEllipsLS;
        radiusNodeLS = sqrt(nodeLS(:,1).^2+nodeLS(:,2).^2);
        nodeLS = repmat(radiusNodeLS,1,2) .* ...
                                      [cos(thetaNodeLS) sin(thetaNodeLS)];
        nodeLS = nodeLS + repmat(fixedPointLS,size(nodeLS,1),1);
        
        polis{1} = nodeLS;
        levelSet = makeLS( X, polis, opts.cutLevelSet );
    case 9
        fixedPointLS = [.698 -.7525];
        anglLS = 150/180*pi;
%        fixedPointLS = [.7 -.7];
%        anglLS = 153/180*pi;
        axisEllipsLS = [cos(anglLS) sin(anglLS)];
        numberNodesLS = 200;
        radiusLongLS = .30;
        radiusShortLS = .1;
        
        axisEllipsLS = axisEllipsLS/norm(axisEllipsLS);
        thetaLS = 2*pi*[numberNodesLS-1:-1:0]'/numberNodesLS;
        nodeLS = [radiusLongLS*cos(thetaLS) radiusShortLS*sin(thetaLS)];
        nodeLS(:,1) = nodeLS(:,1) + radiusLongLS;
        if axisEllipsLS(1)==0;
            thetaEllipsLS = sign(axisEllipsLS(2))*pi/2;
        else
            thetaEllipsLS = atan(axisEllipsLS(2)/axisEllipsLS(1));
        end
        if axisEllipsLS(1)<0;
            thetaEllipsLS = thetaEllipsLS +pi;
        end
        thetaNodeLS = zeros(size(nodeLS,1),1);
        ind1 = find(nodeLS(:,1)==0);
        thetaNodeLS(ind1) = sign(nodeLS(ind1,1))*pi/2;
        ind1 = find(nodeLS(:,1)~=0);
        thetaNodeLS(ind1) = atan(nodeLS(ind1,2)./nodeLS(ind1,1));
        thetaNodeLS = thetaNodeLS + thetaEllipsLS;
        radiusNodeLS = sqrt(nodeLS(:,1).^2+nodeLS(:,2).^2);
        nodeLS = repmat(radiusNodeLS,1,2) .* ...
                                      [cos(thetaNodeLS) sin(thetaNodeLS)];
        nodeLS = nodeLS + repmat(fixedPointLS,size(nodeLS,1),1);
        
        polis{1} = nodeLS;
        levelSet = makeLS( X, polis, opts.cutLevelSet );
    case 10
        fixedPointLS = [1.2 -0.9];
        anglLS = 153/180*pi;
        axisEllipsLS = [cos(anglLS) sin(anglLS)];
        numberNodesLS = 200;
        radiusLongLS = .90;
        radiusShortLS = .2;
        
        axisEllipsLS = axisEllipsLS/norm(axisEllipsLS);
        thetaLS = 2*pi*[numberNodesLS-1:-1:0]'/numberNodesLS;
        nodeLS = [radiusLongLS*cos(thetaLS) radiusShortLS*sin(thetaLS)];
        nodeLS(:,1) = nodeLS(:,1) + radiusLongLS;
        if axisEllipsLS(1)==0;
            thetaEllipsLS = sign(axisEllipsLS(2))*pi/2;
        else
            thetaEllipsLS = atan(axisEllipsLS(2)/axisEllipsLS(1));
        end
        if axisEllipsLS(1)<0;
            thetaEllipsLS = thetaEllipsLS +pi;
        end
        thetaNodeLS = zeros(size(nodeLS,1),1);
        ind1 = find(nodeLS(:,1)==0);
        thetaNodeLS(ind1) = sign(nodeLS(ind1,1))*pi/2;
        ind1 = find(nodeLS(:,1)~=0);
        thetaNodeLS(ind1) = atan(nodeLS(ind1,2)./nodeLS(ind1,1));
        thetaNodeLS = thetaNodeLS + thetaEllipsLS;
        radiusNodeLS = sqrt(nodeLS(:,1).^2+nodeLS(:,2).^2);
        nodeLS = repmat(radiusNodeLS,1,2) .* ...
                                      [cos(thetaNodeLS) sin(thetaNodeLS)];
        nodeLS = nodeLS + repmat(fixedPointLS,size(nodeLS,1),1);
        
        polis{1} = nodeLS;
        levelSet = makeLS( X, polis, opts.cutLevelSet );
    case 11
        fixedPointLS = [.7 -.7];
        anglLS = 153/180*pi;
        axisEllipsLS = [cos(anglLS) sin(anglLS)];
        numberNodesLS = 200;
        radiusLongLS = .60;
        radiusShortLS = .1;
        
        axisEllipsLS = axisEllipsLS/norm(axisEllipsLS);
        thetaLS = 2*pi*[numberNodesLS-1:-1:0]'/numberNodesLS;
        nodeLS = [radiusLongLS*cos(thetaLS) radiusShortLS*sin(thetaLS)];
        nodeLS(:,1) = nodeLS(:,1) + radiusLongLS;
        if axisEllipsLS(1)==0;
            thetaEllipsLS = sign(axisEllipsLS(2))*pi/2;
        else
            thetaEllipsLS = atan(axisEllipsLS(2)/axisEllipsLS(1));
        end
        if axisEllipsLS(1)<0;
            thetaEllipsLS = thetaEllipsLS +pi;
        end
        thetaNodeLS = zeros(size(nodeLS,1),1);
        ind1 = find(nodeLS(:,1)==0);
        thetaNodeLS(ind1) = sign(nodeLS(ind1,1))*pi/2;
        ind1 = find(nodeLS(:,1)~=0);
        thetaNodeLS(ind1) = atan(nodeLS(ind1,2)./nodeLS(ind1,1));
        thetaNodeLS = thetaNodeLS + thetaEllipsLS;
        radiusNodeLS = sqrt(nodeLS(:,1).^2+nodeLS(:,2).^2);
        nodeLS = repmat(radiusNodeLS,1,2) .* ...
                                      [cos(thetaNodeLS) sin(thetaNodeLS)];
        nodeLS = nodeLS + repmat(fixedPointLS,size(nodeLS,1),1);
        
        polis{1} = nodeLS;
        levelSet = makeLS( X, polis, opts.cutLevelSet );
    case 12
        fixedPointLS = [.9 -.7];
        anglLS = 135/180*pi;
        axisEllipsLS = [cos(anglLS) sin(anglLS)];
        numberNodesLS = 100;
        radiusLongLS = .04;
        radiusShortLS = .025;
        
        axisEllipsLS = axisEllipsLS/norm(axisEllipsLS);
        thetaLS = 2*pi*[numberNodesLS-1:-1:0]'/numberNodesLS;
        nodeLS = [radiusLongLS*cos(thetaLS) radiusShortLS*sin(thetaLS)];
        nodeLS(:,1) = nodeLS(:,1) + radiusLongLS;
        if axisEllipsLS(1)==0;
            thetaEllipsLS = sign(axisEllipsLS(2))*pi/2;
        else
            thetaEllipsLS = atan(axisEllipsLS(2)/axisEllipsLS(1));
        end
        if axisEllipsLS(1)<0;
            thetaEllipsLS = thetaEllipsLS +pi;
        end
        thetaNodeLS = zeros(size(nodeLS,1),1);
        ind1 = find(nodeLS(:,1)==0);
        thetaNodeLS(ind1) = sign(nodeLS(ind1,1))*pi/2;
        ind1 = find(nodeLS(:,1)~=0);
        thetaNodeLS(ind1) = atan(nodeLS(ind1,2)./nodeLS(ind1,1));
        thetaNodeLS = thetaNodeLS + thetaEllipsLS;
        radiusNodeLS = sqrt(nodeLS(:,1).^2+nodeLS(:,2).^2);
        nodeLS = repmat(radiusNodeLS,1,2) .* ...
                                      [cos(thetaNodeLS) sin(thetaNodeLS)];
        nodeLS = nodeLS + repmat(fixedPointLS,size(nodeLS,1),1);
        
        polis{1} = nodeLS;
        levelSet = makeLS( X, polis, opts.cutLevelSet );
    case 13
        % Horizontal Tube at y in [-.2 +.2]
        tubey = 0.2;
        margen = 2*max( y_up, x_up );
        % Points
        poli1 = [x_lo-margen -tubey ; ...
                 x_up+margen -tubey ;...
                 x_up+margen  tubey ;...
                 x_lo-margen  tubey ;...
                 x_lo-margen -tubey ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );

    case 14
        % Sine Tube at y in [-.2 +.2]
        xx = (x_lo:0.03:x_up)';
        tubey = 0.2+0.05*cos(4*pi*xx);
        tubeY = max(tubey);
        margen = 2*max( y_up, x_up );
        % Points
        poli1 = [x_lo-margen  -tubeY 
                 xx           -tubey
                 x_up+margen  -tubeY
                 x_up+margen   tubeY 
                 xx(end:-1:1)  tubey(end:-1:1)
                 x_lo-margen   tubeY
                 x_lo-margen  -tubeY          ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );

    case 15
        % Horizontal Split
        midy = 0.55;
        margen = 2*max( y_up, x_up );
        r = 0.2;
        theta = linspace(-pi/2,pi/2,40);
        xt = 0.5 + r * sin(theta);
        yt = 0.6 - r * cos(theta);
        ind = (yt <= midy);
        % Points
        poli1 = [x_lo-margen    midy ; ...
                 xt(ind)'       yt(ind)'; ...
                 x_up+margen    midy ;...
                 x_up+margen    y_up+margen ;...
                 x_lo-margen    y_up+margen ;...
                 x_lo-margen    midy ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );

    case 16
        % Horizontal Split
        midy = 0.5;
        midy2 = .4375;
        margen = 2*max( y_up, x_up );
        % Points
        poli1 = [x_lo-margen    midy ; ...
                 .3125          midy ; ...
                 .3750          midy2; ...
                 .6250          midy2; ...
                 .6875          midy ; ...
                 x_up+margen    midy ;...
                 x_up+margen    y_up+margen ;...
                 x_lo-margen    y_up+margen ;...
                 x_lo-margen    midy ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );

    case 17
        % circle of radius .25
        R = .25;
        NR = 360;
        theta = linspace(0,2*pi,NR+1)';
        poli1 = R*[ cos(theta) sin(theta) ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );

    case 18
        % circle of radius .5
        R = .5;
        NR = 360;
        theta = linspace(0,2*pi,NR+1)';
        poli1 = R*[ cos(theta) sin(theta) ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );

    case 19
        % circle of radius .75
        R = .75;
        NR = 360;
        theta = linspace(0,2*pi,NR+1)';
        poli1 = R*[ cos(theta) sin(theta) ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );

    case 110
        angle = 14;
        midy = 0.45*y_up-0.32;
        margen = 2*max( y_up, x_up );
        % Points
        delta = (margen + 0.5*(x_up - x_lo) )*tan(angle*pi/180);
        poli1 = [x_lo-margen    midy+delta ; ...
                 x_up+margen    midy-delta ;...
                 x_up+margen    y_up+margen ;...
                 x_lo-margen    y_up+margen ;...
                 x_lo-margen    midy ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );

   case 111
        poli1 = [-0.8 1
                 1 -0.8
                 10 -0.8
                 10 10
                 -0.8 10
                 -0.8 1];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );
        
   case 112
        poli1 = [-1 0.3
           -0.4 0.3
           0.4 -0.35
           1 -0.35
           10 -0.35
           10 10
           -10 10
           -10 0.3
           -1 0.3];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );
        
    case 113
        % Horizontal Split
        midy = 0.1;
        margen = 2*max( y_up, x_up );
        % Points
        poli1 = [x_lo-margen    midy ; ...
                 x_up+margen    midy ;...
                 x_up+margen    y_up+margen ;...
                 x_lo-margen    y_up+margen ;...
                 x_lo-margen    midy ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );

   case 114
        midy = 0.25/2;
        margen = 2*max( y_up, x_up );
        % Points
        cc = 0.1;
        poli1 = [x_lo-margen    midy+cc ; ...
                 x_up+margen    midy-cc ; ...
                 x_up+margen    y_up+margen ;...
                 x_lo-margen    y_up+margen ;...
                 x_lo-margen    midy+cc ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );
        
    case 115
        % Horizontal Split
        midy = 0.25/2;
        margen = 2*max( y_up, x_up );
        % Points
        poli1 = [x_lo-margen    midy ; ...
                 x_up+margen    midy ;...
                 x_up+margen    y_up+margen ;...
                 x_lo-margen    y_up+margen ;...
                 x_lo-margen    midy ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );
        
   case 116
        midy = 0.1;
        margen = 2*max( y_up, x_up );
        % Points
        cc = 2/8 * 3;
        poli1 = [x_lo-margen    midy+cc ; ...
                 x_up+margen    midy-cc ; ...
                 x_up+margen    y_up+margen ;...
                 x_lo-margen    y_up+margen ;...
                 x_lo-margen    midy+cc ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );
        
   case 117
        midy = 0.1;
        margen = 2*max( y_up, x_up );
        % Points
        cc = 2/8 * 1;
        poli1 = [x_lo-margen    midy+cc ; ...
                 x_up+margen    midy-cc ; ...
                 x_up+margen    y_up+margen ;...
                 x_lo-margen    y_up+margen ;...
                 x_lo-margen    midy+cc ];
        polis{1} = poli1;
        levelSet = mkLS( X, polis, opts.cutLevelSet );
        
    otherwise
        error( 'predefinedLevelSet: unknow case number' )
end
% Make it a column vector, just in case
levelSet = reshape( levelSet, length( levelSet ), 1 );
