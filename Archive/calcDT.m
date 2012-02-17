function [dt,maxVel] = calcDT( vel, h )
%
%
%

maxVel = max( vel );

theta = 0.1;
dt = sqrt(theta * 3/4) * h / maxVel;
