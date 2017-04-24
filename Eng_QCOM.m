%% This script plots the output from the 2d QCOM model
%% ATMOS 6150
%% Byron Eng
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

animate = true;
plotKE = false;
plotPROFILES = false;
normalize = true; %also makes the clouds look better
cloud = true; %currently only applies to the static plot

%% Read in the data
v = dlmread('v.dat');
w = dlmread('w.dat');
theta = dlmread('theta.dat');
Pi = dlmread('pi.dat');
qc = dlmread('qc.dat');
qc(1,:) = 0;
qc(end,:) = max(qc);

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
if cloud
    ch(64,1:3) = 1;
    hold on
    h = pcolor(ones(12,22));
    alpha(h,(qc))
    hold off
end% if cloud
colormap(ch)
shading flat
if normalize
    title('v','FontSize',16)
else
    title('v [m/s]','FontSize',16)
end


subplot(2,2,2)
contourf(w)
if normalize
    set(gca, 'Clim', [-1 1])
end% if normalize
if cloud
    hold on
    h = pcolor(ones(12,22));
    alpha(h,(qc))
    hold off
end %if cloud
colorbar
if normalize
    title('w','FontSize',16)
else
    title('w [m/s]','FontSize',16)
end
shading flat

subplot(2,2,3)
contourf(theta)
if normalize
    set(gca, 'Clim', [-1 1])
end% if normalize
if cloud
    hold on
    h = pcolor(ones(12,22));
    alpha(h,(qc))
    hold off
end %if cloud
colorbar
if normalize
    title('\theta_v - \theta_0','FontSize',16)
else
    title('\theta_v - \theta_0 [K]','FontSize',16)
end
shading flat

subplot(2,2,4)
contourf(Pi)
if normalize
    set(gca, 'Clim', [-1 1])
end% if normalize
if cloud
    hold on
    h = pcolor(ones(12,22));
    alpha(h,(qc))
    hold off
end %if cloud
colorbar
if normalize
    title('\pi\prime','FontSize',16)
else
    title('\pi\prime [Pa]','FontSize',16)
end
shading flat

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
    aqc = dlmread('aqc.dat');
    
    %Normalize
    if normalize
    av = av/max(max(v));
    aw = aw/max(max(w));
    end %if normalize
    
    gridht = 102;
    nframes = size(av,1)/gridht;
    filename = 'animated.gif';
    
    figure('OuterPosition',[0 0 1414 1000]) 
    
    for i=1:nframes
        aqc((((i-1)*gridht) + 1),:) = 0;
        aqc(i*gridht,:) = max(max(aqc));
        subplot(2,2,1)
            avnow = av((i-1)*gridht + (1:gridht),:);
            if normalize
                avnow = avnow/max(max(avnow));
            end %if normalize
            contourf(avnow)
            if normalize
                set(gca, 'Clim', [-1 1])
            end% if normalize
            colorbar
            ch = colormap;
            if cloud
                ch(64,1:3) = 1;
                hold on
                h = pcolor(ones(12,22));
                alpha(h,(aqc((i-1)*gridht + (1:gridht),:)))
                hold off  
            end %if cloud
            colormap(ch)
            title('v')  
            shading flat

        subplot(2,2,2)
            awnow = aw((i-1)*gridht + (1:gridht),:);
            if normalize
                awnow = awnow/max(max(awnow));
            end %if normalize
            contourf(awnow)
            if normalize
                set(gca, 'Clim', [-1 1])
            end% if normalize
            colorbar
            ch = colormap;
            if cloud
                ch(64,1:3) = 1;
                hold on
                h = pcolor(ones(12,22));
                alpha(h,(aqc((i-1)*gridht + (1:gridht),:)))
                hold off
            end %if cloud
            colormap(ch)
            title('w')
            shading flat

        subplot(2,2,3)
            athetanow = atheta((i-1)*gridht + (1:gridht),:);
            if normalize
                athetanow = athetanow/max(max(athetanow));
            end %if normalize
            contourf(athetanow)
            if normalize
                set(gca, 'Clim', [-1 1])
            end% if normalize
            colorbar
            ch = colormap;
            if cloud
                ch(64,1:3) = 1;
                hold on
                h = pcolor(ones(12,22));
                alpha(h,(aqc((i-1)*gridht + (1:gridht),:)))
                hold off
            end %if cloud
            colormap(ch)
            title('\theta_v - \theta_0')
            shading flat

        subplot(2,2,4)
            if normalize
                aPi((i-1)*gridht + (1:gridht),:) = ...
                    aPi((i-1)*gridht + (1:gridht),:)/...
                    max(max(aPi((i-1)*gridht + (1:gridht),:)));
            end %if normalize
            contourf(aPi((i-1)*gridht + (1:gridht),:))
            if normalize
                set(gca, 'Clim', [-1 1])
            end% if normalize
            colorbar
            ch = colormap;
            if cloud
            ch(64,1:3) = 1;
            hold on
            h = pcolor(ones(12,22));
            alpha(h,(aqc((i-1)*gridht + (1:gridht),:)))
            hold off
            end %if cloud
            colormap(ch)
            title('\pi')
            shading flat
            
      drawnow
      frame = getframe(gcf);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if i == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end %if-else
      
      close
        figure('OuterPosition',[0 0 1414 1000])
      
    end %for loop    
end %if animate
