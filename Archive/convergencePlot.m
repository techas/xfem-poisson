function convergencePlotForPaper()
%
%
%
clc
datafem = [];
dataxfem = [];
dataxfemp = [];
datah = 2./[16 32 64 128];

%% MESH 2
% fem
load datafem-um2.mat
someStuffForTheConvergencePlot
mm = maxFluxJump(polis,X,T,levelSet,h,hE,opts);
datafem = [datafem mm];
% xfem
load dataxfem-um2.mat
mm = maxFluxJump(polis,X,T,levelSet,h,hE,opts);
dataxfem = [dataxfem mm];
% xfemp
load dataxfemp-um2.mat
mm = maxFluxJump(polis,X,T,levelSet,h,hE,opts);
dataxfemp = [dataxfemp mm];

%% MESH 4
% fem
load datafem-um4.mat
someStuffForTheConvergencePlot
mm = maxFluxJump(polis,X,T,levelSet,h,hE,opts);
datafem = [datafem mm];
% xfem
load dataxfem-um4.mat
mm = maxFluxJump(polis,X,T,levelSet,h,hE,opts);
dataxfem = [dataxfem mm];
% xfemp
load dataxfemp-um4.mat
mm = maxFluxJump(polis,X,T,levelSet,h,hE,opts);
dataxfemp = [dataxfemp mm];

%% MESH 8
% fem
load datafem-um8.mat
someStuffForTheConvergencePlot
mm = maxFluxJump(polis,X,T,levelSet,h,hE,opts);
datafem = [datafem mm];
% xfem
load dataxfem-um8.mat
mm = maxFluxJump(polis,X,T,levelSet,h,hE,opts);
dataxfem = [dataxfem mm];
% xfemp
load dataxfemp-um8.mat
mm = maxFluxJump(polis,X,T,levelSet,h,hE,opts);
dataxfemp = [dataxfemp mm];

%% MESH 16
% fem
load datafem-um16.mat
someStuffForTheConvergencePlot
mm = maxFluxJump(polis,X,T,levelSet,h,hE,opts);
datafem = [datafem mm];
% xfem
load dataxfem-um16.mat
mm = maxFluxJump(polis,X,T,levelSet,h,hE,opts);
dataxfem = [dataxfem mm];
% xfemp
load dataxfemp-um16.mat
mm = maxFluxJump(polis,X,T,levelSet,h,hE,opts);
dataxfemp = [dataxfemp mm];

save('tmpData.mat', 'datafem', 'dataxfem', 'dataxfemp', 'datah')

figure(88), clf
subplot(311), plot(datah,datafem,'o')
set(gca,'xscale', 'log', 'yscale', 'log')
subplot(312), plot(datah,dataxfem,'o')
set(gca,'xscale', 'log', 'yscale', 'log')
subplot(313), plot(datah,dataxfemp,'o')
set(gca,'xscale', 'log', 'yscale', 'log')


function out = maxFluxJump(polis,X,T,levelSet,h,hE,opts)
%
%
%
np = 10;
mm = [];
for p = 1:length(polis)
   for s = 1:length(polis{p})-1
      p1 = polis{p}(s,:);
      p2 = polis{p}(s+1,:);
      dp = p1-p2;
      pos = [linspace(p1(1),p2(1),np)' linspace(p1(2),p2(2),np)'];
      % flux
      qp = FluxosX( X, T, levelSet, h, hE, pos, 1, opts.tolerance );
      qn = FluxosX( X, T, levelSet, h, hE, pos, 0, opts.tolerance );
      % delta flux
      dq = qp - qn;
      % normal to the segment (linear levelset assumed!)
      nn = [-dp(2) dp(1)];
      nn = nn/norm(nn);
      % flux jump
      normalq = sum(dq.*repmat(nn,size(dq,1),1),2);
      mm = [mm; normalq];
   end
end

[out,ix] = max(mm);
fprintf('max in %f of %i (p=%i)\n',ix/np,length(polis{1}),length(polis))

out = norm(mm);