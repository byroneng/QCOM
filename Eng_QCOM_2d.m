%% This script plots the output from the 2d QCOM model
%% ATMOS 6150
%% Byron Eng
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

theory = false; %Theoretical solutions

v = dlmread('v.dat');
w = dlmread('w.dat');
theta = dlmread('theta.dat');
Pi = dlmread('pi.dat');

%Theoretical solution
if theory
    z = 0:.01:1;
    y = 0:.01:(2*sqrt(2));
    kc = pi/sqrt(2);
    wo = sin(pi*z);
    kw = 100;
    
    w = kron(wo',cos(kc*y));
    v = -(1/kc)*kron((diff(wo)/.01)',sin(kc*y));
    theta = kron(((kw/(kc*kc))*((diff(diff(diff(diff(wo)/.01)/.01)/.01)/.01)+((kc^4)*wo(3:length(wo)-2))))',cos(kc*y));
    Pi = kron((kw/(kc*kc))*((diff(diff(diff(wo)/.01)/.01)/.01)-(kc*kc*(diff(wo(2:length(wo)-1))/.01)))',cos(kc*y));
    
end

%Normalize
v = v/max(max(v));
w = w/max(max(w));
theta = theta/max(max(theta));
Pi = Pi/max(max(Pi));

figure('OuterPosition',[0 0 900 800])

subplot(2,2,1)
contour(v)
colorbar
title('v')

subplot(2,2,2)
contour(w)
colorbar
title('w')

subplot(2,2,3)
contour(theta)
colorbar
title('theta')

subplot(2,2,4)
contour(Pi)
colorbar
title('\pi')
