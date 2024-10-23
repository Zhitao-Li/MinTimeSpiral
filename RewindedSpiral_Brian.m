clear;
clc;
close all;

addpath('Lib/');
addpath('Lib/socp/bin/');

if ~exist('Result', 'dir')
    mkdir('Result');
end

%%
secondMomentNull = false;

headerFilePath = 'Results/SpiralTraj.h';

%%
T = 0.000001; % Second
slewR = 15000; % Gauss/cm/s
maxAmp = 3.7; % Gauss/cm

alpha1 = 0.9;
alpha2 = 0.6;

NumArms = 48;
FOV0 = 22; % cm
FOV1 = 0;
FOV2 = 0;
baseRes = 256;
inplaneRes = 10.0.*FOV0./baseRes; % mm

%%
[~, g] = vds(alpha1*slewR, maxAmp, T, ...
             NumArms, ...
             FOV0, FOV1, FOV2, ...
             1/(0.2*inplaneRes));

nstep = 10;

%%
gsp = [real(g(:)), imag(g(:))];
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


plotgradinfo([gsp1; g1], nstep*T);
createCPPHeader(gsp1, g1, alpha1*slewR, alpha2*slewR, ...
                NumArms, ...
                FOV0, FOV1, FOV2, baseRes, inplaneRes, ...
                T, nstep, maxAmp, headerFilePath);





