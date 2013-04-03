function fh = physio_plot_cropped_phys_to_acqwindow(ons_secs, sqpar)
% plot parts of the time series to be processed into regressors
%
% USAGE
%   fh = physio_plot_cropped_phys_to_acqwindow(ons_secs, sqpar, y)
%
% INPUT
%   ons_secs    - output of physio_crop_scanphysevents_to_acq_window
%   sqpar       - output of physio_crop_scanphysevents_to_acq_window
%
% OUTPUT
%   fh          figure handle of output figure
% 
% Author: Lars Kasper
%
% Copyright (C) 2013, Institute for Biomedical Engineering, ETH/Uni Zurich.
%
% This file is part of the PhysIO toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$
[fh, MyColors] = physio_get_default_fig_params();
set(fh,'Name','Cutout actual scans - all events and gradients');

Ndummies    = sqpar.Ndummies;
Nslices     = sqpar.Nslices;
sampling    = 1;

%% 1. Plot uncropped data
t        = ons_secs.raw.t;
cpulse      = ons_secs.raw.cpulse;
r           = ons_secs.raw.r;
c           = ons_secs.raw.c;
spulse      = ons_secs.raw.spulse;
svolpulse   = ons_secs.raw.svolpulse;

maxr = max(abs(r));
ampsv = maxr/2.25;
amps = maxr / 3;
ampc = maxr / 2;

y = [c, r];
x = y (1:sampling:end, :);
stem(spulse(1:Ndummies*Nslices),amps*ones(Ndummies*Nslices,1),'r--');
hold on;
stem(svolpulse(Ndummies+1:end),ampsv*ones(length(svolpulse)-Ndummies,1),'--g', 'LineWidth',2);
stem(spulse((Ndummies*Nslices+1):end), amps*ones(length(spulse)-Ndummies*Nslices,1), 'c--') ;
stem(cpulse, ampc*ones(length(cpulse),1), 'm--') ;
plot(t(1:sampling:end), x, '--');


%% 2. Plot cropped data

t        = ons_secs.t;
cpulse      = ons_secs.cpulse;
r           = ons_secs.r;
c           = ons_secs.c;
spulse      = ons_secs.spulse;
svolpulse   = ons_secs.svolpulse;



%plot physiological regressors and scan events

y = [c, r];
x = y (1:sampling:end, :);
hs(1) = stem(spulse(1:Ndummies*Nslices),amps*ones(Ndummies*Nslices,1),'r');
hold on;
hs(end+1) = stem(svolpulse(Ndummies+1:end),ampsv*ones(length(svolpulse)-Ndummies,1),'g', 'LineWidth',2);
hs(end+1) = stem(spulse((Ndummies*Nslices+1):end), amps*ones(length(spulse)-Ndummies*Nslices,1), 'c') ;
hs(end+1) = stem(cpulse, ampc*ones(length(cpulse),1), 'm') ;
hs2 = plot(t(1:sampling:end), x, '-')';
hs = [hs, hs2];
hs(end+1) = plot(t,r,'ko');

xlabel('t (s)'); ylabel('Amplitude (a. u.)');
title('Cutout region for physiological regressors');
legend( hs, {['dummy scan event marker (N = ' int2str(Ndummies*Nslices) ')'], ...
    ['volume event marker (N = ' int2str(length(svolpulse)-Ndummies) '), without dummies'], ...
    ['scan event marker (N = ' int2str(length(spulse)-Ndummies*Nslices) ')'], ...
    ['ECG event marker (N = ' int2str(length(cpulse)) ')'], ...
    'filtered ECG', 'resp signal', 'used resp signal'});
ylim(1.2*maxr*[-1 1]);
end