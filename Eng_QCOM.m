%% This script plots the output from the 2d QCOM model
%% ATMOS 6150
%% Byron Eng
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

animate = false;
plotKE = false;
plotPROFILES = false;
normalize = true; %also makes the clouds look better

%% Read in the data
v = dlmread('v.dat');
w = dlmread('w.dat');
theta = dlmread('theta.dat');
Pi = dlmread('pi.dat');
qc = dlmread('qc.dat');
qctest = (1:12)*.0001;
qctest = repmat(qctest',1,22);
qctest = flipud(qctest);

if plotKE
    tv = dlmread('tv.dat');
    tw = dlmread('tw.dat');
    ttheta = dlmread('ttheta.dat');
    tPi = dlmread('tpi.dat');
end

%Normalize
if normalize
v = v/max(max(v));
w = w/max(max(w));
theta = theta/max(max(theta));
Pi = Pi/max(max(Pi));
end %if normalize

%%Plots
figure('OuterPosition',[0 0 900 800])

subplot(2,2,1)
contourf(v);
if normalize
    set(gca, 'Clim', [-1 1])
end% if normalize
colorbar
ch = colormap;
ch(64,1:3) = 1;
colormap(ch)
title('v')
hold on
h = pcolor(ones(12,22));
alpha(h,(1-qctest))
shading flat
hold off


subplot(2,2,2)
contourf(w)
if normalize
    set(gca, 'Clim', [-1 1])
end% if normalize
colorbar
title('w')
hold on
h = pcolor(ones(12,22));
alpha(h,(1-qctest))
shading flat
hold off

subplot(2,2,3)
contourf(theta)
if normalize
    set(gca, 'Clim', [-1 1])
end% if normalize
colorbar
title('theta')
hold on
h = pcolor(ones(12,22));
alpha(h,(1-qctest))
shading flat
hold off

subplot(2,2,4)
contourf(Pi)
if normalize
    set(gca, 'Clim', [-1 1])
end% if normalize
colorbar
title('\pi')
hold on
h = pcolor(ones(12,22));
alpha(h,(1-qctest))
shading flat
hold off

if plotKE
%Kinetic Energy
i=1;
ke = 0*tv;
for i=1:length(tv)
    ke(i) = (tv(i)^2 + tw(i)^2)/numel(theta);
end

figure
plot(100*(1:length(ke)),ke)
xlabel('Time [s]')
ylabel('Kinetic Energy')
title('Mean Kinetic Energy Profile')
end %if plotKE

if plotPROFILES
figure('OuterPosition',[0 0 900 800])

%theta
kth = 50;
subplot(2,2,1)
mtheta = mean(theta,2);
plot(mtheta,kth*((1:length(mtheta))-1))
title('Mean \theta Perturbation [K]')
ylabel('Height [m]')
xlabel('\theta [K]')
ylim([0 500])


%conductive vertical heat flux
cdhf = -diff(mtheta) * kth;
tohf = repmat(cdhf(1),1,length(cdhf));
cvhf = tohf-cdhf';

subplot(2,2,2)
plot(cdhf,kth*((1:length(cdhf))-1))
title('Conductive Vertical Heat Flux')
ylabel('Height')
xlabel('Heat Flux')
xlim([-20 max(cvhf)+5])
ylim([0 500])

%Convective vertical heat flux

% cvhf = mean(w.*theta,2);

subplot(2,2,3)
plot(cvhf,kth*((1:length(cvhf))-1))
title('Convective Vertical Heat Flux')
ylabel('Height [m]')
xlabel('Heat Flux')
xlim([-20 max(cvhf)+5])
ylim([0 500])

%Total heat flux
subplot(2,2,4)
% plot(cdhf+cvhf(1:length(cdhf)),50*(1:length(cdhf)))
plot(tohf,kth*((1:length(tohf))-1))
title('Total Heat Flux')
ylabel('Height [m]')
xlabel('Heat Flux')
xlim([-20 max(cvhf)+5])
ylim([0 500])

end %if plotPROFILES


if animate
    %use av.dat, aw.dat, etc.
    %use entire domain + boundaries

    av = dlmread('av.dat');
    aw = dlmread('aw.dat');
    atheta = dlmread('atheta.dat');
    aPi = dlmread('api.dat');
    
    gridht = 12;
    nframes = size(av,1)/gridht;
    filename = 'animated.gif';
    
    figure(4)
    
    for i=1:nframes
        subplot(2,2,1)
            contour(av((i-1)*gridht + (1:gridht),:))
%            caxis([-1.5 1.5])
            colorbar
            title('v')

        subplot(2,2,2)
            contour(aw((i-1)*gridht + (1:gridht),:))
%            caxis([-4 4])
            colorbar
            title('w')

        subplot(2,2,3)
            contour(atheta((i-1)*gridht + (1:gridht),:))
%            caxis([-10 10])
            colorbar
            title('theta')

        subplot(2,2,4)
            contour(aPi((i-1)*gridht + (1:gridht),:))
%            caxis([0 .05])
            colorbar
            title('\pi')
            
      drawnow
      frame = getframe(4);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if i == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end %if-else
      
    end %for loop    
end %if animate
