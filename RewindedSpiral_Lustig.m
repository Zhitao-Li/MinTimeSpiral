clear;
clc;
close all;

addpath('Lib/');
addpath('Lib/socp/bin/');
addpath('Lib/minTimeGradientMex/bin/');

if ~exist('Results/', 'dir')
    mkdir('Results');
end

%%
useInterpolate = false;
secondMomentNull = false;

headerFilePath = 'Results/SpiralTraj.h';

%%
T = 0.000001; % Second
slewR = 15000; % Gauss/cm/s
maxAmp = 3.7; % Gauss/cm

alpha1 = 0.9;
alpha2 = 0.6;

NumArms = 48;
FOV = [40, 40, 20, 20]; % cm
kRad = [0, 0.2, 0.3, 1]; % Normalized
inplaneRes = 1.17; % mm
baseRes = ceil(mean(FOV)*10/inplaneRes);

%%
[~, gsp, ~, ~, ~] = ...
    vdSpiralDesign(0, NumArms, inplaneRes, ...
                   FOV, kRad, ...
                   maxAmp, alpha1*slewR/1000, T*1e3, 'pchip');

nstep = 10;

%%
gsp = gsp(:, 1:2);
gsp1 = gsp(1:nstep:end, :);
[k, ~, s, m1, m2] = calcgradinfo(gsp1, nstep*T);
sz = size(k);
t0 = nstep*T*sz(1);

if secondMomentNull
    g1 = mintimegrad(60, ...
                     gsp1(sz(1), :), [0, 0], ...
                     [-k(sz(1), :); -m1(sz(1), :); -m2(sz(1), :)], ...
                     nstep*T, maxAmp, alpha2*slewR, t0,...
                     3);
else
    g1 = mintimegrad(60, ...
                     gsp1(sz(1), :), [0, 0], ...
                     [-k(sz(1), :); -m1(sz(1), :)], ...
                     nstep*T, maxAmp, alpha2*slewR, t0,...
                     3);
end


if useInterpolate
    g = interp1(linspace(1, length(g1), length(g1)), g1, ...
                linspace(1, length(g1), nstep*length(g1)), ...
                'pchip');
    plotgradinfo([gsp; g], T);
    createCPPHeader2(gsp, g, alpha1*slewR, alpha2*slewR,...
                     NumArms, ...
                     FOV, kRad, baseRes, inplaneRes, ...
                     T, maxAmp, nstep, headerFilePath);
else % If not using interpolation, then the gradient dwell time will increase
    plotgradinfo([gsp1; g1], nstep*T);
    createCPPHeader2(gsp, g1, alpha1*slewR, alpha2*slewR,...
                     NumArms, ...
                     FOV, kRad, baseRes, inplaneRes, ...
                     T, maxAmp, nstep, headerFilePath);
end
