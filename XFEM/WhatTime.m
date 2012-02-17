function t = WhatTime( i, time)
% WHATTIME to get the value of time t at step i
%
%  syntax: t = WhatTime( i, time)
%
%  i:    time step index [integer]
%  time: time structured array
%
%  t:    value of time in seconds [scalar]

% R. Cottereau 04/2008

t = time.initial + i * time.step;