function [normGradLS,gradLS] = GradientLevelSetNodal( X, T, levelSet, ...
                                                     N, Nxi, Neta, opts )
% Computation of the nodal gradient of the level set function. A least
% square fit is performed over a patch of elements, and only the UpWind
% Gauss points are considered for the fit
%
%  syntax: [normGradLS,gradLS] = GradientLevelSetNodal( X, T, levelSet, ...
%                                                     N, Nxi, Neta, opts )
%
%  X, T: nodal coordinates and connectivity matrix
%  levelSet: level set function given as a list of nodal values
%  N, Nxi, Neta: shape function and derivatives on the reference element
%  opts: options structured array
%
%  gradLS: gradient of the level Set at nodes
%  normGradLS: norm of the previous

% R. Cottereau 03/2008

Xt = unique(T(:));
nt = size(Xt,1);
npoin = size(X,1);
ngaus = opts.numberOfGaussPoints;
gradLS = zeros(npoin,2);
normGradLS = zeros(npoin,1);

for i1 = 1:nt
    inode = Xt(i1);
    [EltsInContact,IndexNodeI] = find(T==inode);
    TContact = T(EltsInContact,:);
    NumberEltsInContact = length(EltsInContact);
    dPhi = zeros(ngaus*NumberEltsInContact,2);
    UpWind = zeros(ngaus*NumberEltsInContact,1);
    Xpg = zeros(size(dPhi));
    for ielt=1:NumberEltsInContact
        Te = TContact(ielt,:);
        Xe = X(Te,:);
        LSe = levelSet(Te);
        for igaus = 1:ngaus
            dN = [ Nxi(igaus,:) ; Neta(igaus,:)];
            Jacob = dN * Xe;
            dLS = dN * LSe;
            dPhi((ielt-1)*ngaus+igaus,:) = Jacob \ dLS;
            if ( dN( :, IndexNodeI(ielt) )' * dLS >= 0 )
                UpWind((ielt-1)*ngaus+igaus) = 1;
            end
%            Nxi_igaus = Nxi(igaus,:);
%            Neta_igaus = Neta(igaus,:);
%            Jacob = [Nxi_igaus*(Xe(:,1))	Nxi_igaus*(Xe(:,2))
%                     Neta_igaus*(Xe(:,1))	Neta_igaus*(Xe(:,2))];
%            LS_xi = Nxi_igaus * LSe;
%            LS_eta = Neta_igaus * LSe;
%            dPhi((ielt-1)*ngaus+igaus,:) = (Jacob\[LS_xi;LS_eta])';
%            if (( Nxi_igaus(IndexNodeI(ielt))*LS_xi + ...
%                              Neta_igaus(IndexNodeI(ielt))*LS_eta ) >= 0 )
%                UpWind((ielt-1)*ngaus+igaus) = 1;
%            end
        end
        Xpg( (ielt-1)*ngaus + [1:ngaus],:) = Isopar(Xe,N);
    end
    if sum(UpWind)~=0
        ind = find(UpWind==1);
    else
        ind= 1 : ngaus*NumberEltsInContact ;
    end
    dPhi = dPhi(ind,:);
    Xpg = Xpg(ind,:);
    MatCoord = [ones(size(Xpg,1),1) Xpg];
    gradLS(inode,:) = [1 X(inode,:)]*((MatCoord'*MatCoord)\(MatCoord'*dPhi));
    normGradLS(inode) = norm(gradLS(inode,:));
end
