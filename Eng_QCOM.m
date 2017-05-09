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
cloud = true; %cloud layer


%% Read in the data
params = dlmread('params.dat');
domht = params(1);
gridhtm = params(2);
gridht = params(3)+2;
domwtm = params(4);
gridwt = params(5);
gridwtm = domwtm/gridwt;

v = dlmread('v.dat');
w = dlmread('w.dat');
theta = dlmread('theta.dat');
Pi = dlmread('pi.dat');
qc = dlmread('qc.dat');
thetao = dlmread('thetao.dat');
pio = dlmread('pio.dat');
% qc(1,:) = 0;
% qc(end,:) = max(qc);

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
qc = qc/max(max(qc(3:100,:))); %excludes boundaries
end %if normalize

figure('OuterPosition',[100 100 1100 800])
subplot(1,2,1)
plot(thetao(2:end,1),0:gridhtm:domht)
set(gca,'FontSize',16)
xlabel('\theta_0 (K)','fontsize',25)
ylabel('Height (m)','fontsize',25)
title('Initial \theta Profile','fontsize',25)
xlim([250 300])

subplot(1,2,2)
plot(pio(2:end,1),0:gridhtm:domht)
set(gca,'FontSize',16)
xlabel('\pi_0 (K)','fontsize',25)
ylabel('Height (m)','fontsize',25)
title('Initial \pi Profile','fontsize',25)
xlim([40000 100000])

%%Plots
figure('OuterPosition',[0 0 900 800])

subplot(2,2,1)
contourf(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm, v);
set(gca,'FontSize',12)
xlabel('Domain Width (m)', 'FontSize', 16)
ylabel('Domain Height (m)', 'FontSize', 16)
if normalize
    set(gca, 'Clim', [-1 1])
end% if normalize
colorbar
colormap jet
ch = colormap;
if cloud
    ch(64,1:3) = 1;
    hold on
    h = pcolor(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm,ones(size(qc)));
    alpha(h,(qc))
%    al = get(gcf,'Alphamap');
%    al = 0:max(max(qc));
    al = 0:(.75/63):.75;
    alphamap(al)
    hold off
end% if cloud
colormap(ch)
shading flat
if normalize
    title('v','fontsize',25)
else
    title('v [m/s]','fontsize',25)
end

subplot(2,2,2)
contourf(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm,w)
set(gca,'FontSize',12)
xlabel('Domain Width (m)', 'FontSize', 16)
ylabel('Domain Height (m)', 'FontSize', 16)
if normalize
    set(gca, 'Clim', [-1 1])
end% if normalize
if cloud
    hold on
    h = pcolor(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm,ones(size(qc)));
    alpha(h,(qc))
    hold off
end %if cloud
colorbar
if normalize
    title('w','fontsize',25)
else
    title('w [m/s]','fontsize',25)
end
colormap(ch)
shading flat

subplot(2,2,3)
contourf(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm,theta)
set(gca,'FontSize',12)
xlabel('Domain Width (m)', 'FontSize', 16)
ylabel('Domain Height (m)', 'FontSize', 16)
if normalize
    set(gca, 'Clim', [-1 1])
end% if normalize
if cloud
    hold on
    h = pcolor(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm,ones(size(qc)));
    alpha(h,(qc))
    hold off
end %if cloud
colorbar
colormap(ch)
if normalize
    title('\theta_v^\prime','fontsize',25)
else
    title('\theta_v^\prime [K]','fontsize',25)
end
shading flat

subplot(2,2,4)
contourf(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm,Pi)
set(gca,'FontSize',12)
xlabel('Domain Width (m)', 'FontSize', 16)
ylabel('Domain Height (m)', 'FontSize', 16)
if normalize
    set(gca, 'Clim', [-1 1])
end% if normalize
if cloud
    hold on
    h = pcolor(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm,ones(size(qc)));
    alpha(h,(qc))
    hold off
end %if cloud
colorbar
if normalize
    title('\pi^\prime','fontsize',25)
else
    title('\pi^\prime [Pa]','fontsize',25)
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
xlabel('Time [s]', 'interpreter','latex','fontsize',25)
ylabel('Kinetic Energy', 'interpreter','latex','fontsize',25)
title('Mean Kinetic Energy Profile', 'interpreter','latex','fontsize',25)
end %if plotKE

if plotPROFILES
figure('OuterPosition',[0 0 900 800])

%theta
kth = 50;
subplot(2,2,1)
mtheta = mean(theta,2);
plot(mtheta,kth*((1:length(mtheta))-1))
title('Mean \theta Perturbation [K]', 'interpreter','latex','fontsize',25)
ylabel('Height [m]', 'interpreter','latex','fontsize',25)
xlabel('\theta [K]', 'interpreter','latex','fontsize',25)
ylim([0 500])


%conductive vertical heat flux
cdhf = -diff(mtheta) * kth;
tohf = repmat(cdhf(1),1,length(cdhf));
cvhf = tohf-cdhf';

subplot(2,2,2)
plot(cdhf,kth*((1:length(cdhf))-1))
title('Conductive Vertical Heat Flux', 'interpreter','latex','fontsize',25)
ylabel('Height', 'interpreter','latex','fontsize',25)
xlabel('Heat Flux', 'interpreter','latex','fontsize',25)
xlim([-20 max(cvhf)+5])
ylim([0 500])

%Convective vertical heat flux

% cvhf = mean(w.*theta,2);

subplot(2,2,3)
plot(cvhf,kth*((1:length(cvhf))-1))
title('Convective Vertical Heat Flux', 'interpreter','latex','fontsize',25)
ylabel('Height [m]', 'interpreter','latex','fontsize',25)
xlabel('Heat Flux', 'interpreter','latex','fontsize',25)
xlim([-20 max(cvhf)+5])
ylim([0 500])

%Total heat flux
subplot(2,2,4)
% plot(cdhf+cvhf(1:length(cdhf)),50*(1:length(cdhf)))
plot(tohf,kth*((1:length(tohf))-1))
title('Total Heat Flux', 'interpreter','latex','fontsize',25)
ylabel('Height [m]', 'interpreter','latex','fontsize',25)
xlabel('Heat Flux', 'interpreter','latex','fontsize',25)
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
    av = av/max(max(av));
    aw = aw/max(max(aw));
    aqc = aqc/(max(max(aqc(3:100,:)))); %excludes boundaries
    end %if normalize
    
%     gridht = 12;
    nframes = size(av,1)/gridht;
    filename = 'animated.gif';
    
    figure(4)
    set(4,'position',[200,200,800,700])

    for i=1:nframes
%        aqc((((i-1)*gridht) + 1),:) = 0;
%        aqc(i*gridht,:) = max(max(aqc));
        subplot(2,2,1)
            contourf(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm,av((i-1)*gridht + (1:gridht),:))
            set(gca,'FontSize',12)
            xlabel('Domain Width (m)', 'FontSize', 10)
            ylabel('Domain Height (m)', 'FontSize', 10)
            if normalize
                set(gca, 'Clim', [-1 1])
            end% if normalize
            colorbar
            colormap jet
            ch = colormap;
            if cloud
                aqcnow = aqc((i-1)*gridht + (1:gridht),:);
                ch(64,1:3) = 1;
                hold on
                h = pcolor(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm,ones(size(aqcnow)));
                alpha(h,(aqcnow))
%                al = get(gcf,'Alphamap');
%                al = 0:max(max(aqc))/63:max(max(aqc));
                al = 0:(.75/63):.75;
                alphamap(al)
                hold off  
            end %if cloud
            colormap(ch)
            title('v','fontsize',25)
%             hold on
%             h = pcolor(ones(12,22));
%             alpha(h,(aqc((i-1)*gridht + (1:gridht),:)))
             shading flat
%             hold off    

        subplot(2,2,2)
            contourf(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm,aw((i-1)*gridht + (1:gridht),:))
            set(gca,'FontSize',12)
            xlabel('Domain Width (m)', 'FontSize', 10)
            ylabel('Domain Height (m)', 'FontSize', 10)
            if normalize
                set(gca, 'Clim', [-1 1])
            end% if normalize
            colorbar
            ch = colormap;
            if cloud
                ch(64,1:3) = 1;
                hold on
                h = pcolor(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm,ones(size(aqcnow)));
                alpha(h,(aqcnow))
%                al = get(gcf,'Alphamap');
%                al = 0:max(max(aqc))/63:max(max(aqc));
%                alphamap(al)
                hold off
            end %if cloud
            colormap(ch)
            title('w','fontsize',25)
%             hold on
%             h = pcolor(ones(12,22));
%             alpha(h,(aqc((i-1)*gridht + (1:gridht),:)))
             shading flat
%             hold off

        subplot(2,2,3)
            athetanow = atheta((i-1)*gridht + (1:gridht),:);
            if normalize
                athetanow = athetanow/max(max(athetanow));
            end %if normalize
            contourf(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm,athetanow)
            set(gca,'FontSize',12)
            xlabel('Domain Width (m)', 'FontSize', 10)
            ylabel('Domain Height (m)', 'FontSize', 10)
            if normalize
                set(gca, 'Clim', [-1 1])
            end% if normalize
            colorbar
            ch = colormap;
            if cloud
                ch(64,1:3) = 1;
                hold on
                h = pcolor(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm,ones(size(aqcnow)));
                alpha(h,(aqcnow))
 %               al = get(gcf,'Alphamap');
 %               al = 0:max(max(aqc))/63:max(max(aqc));
 %               alphamap(al)
                hold off
            end %if cloud
            colormap(ch)
            title('\theta_v^\prime','fontsize',25)
%             hold on
%             h = pcolor(ones(12,22));
%             alpha(h,(aqc((i-1)*gridht + (1:gridht),:)))
             shading flat
%             hold off

        subplot(2,2,4)
            if normalize
                aPi((i-1)*gridht + (1:gridht),:) = ...
                    aPi((i-1)*gridht + (1:gridht),:)/...
                    max(max(aPi((i-1)*gridht + (1:gridht),:)));
            end %if normalize
            contourf(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm,aPi((i-1)*gridht + (1:gridht),:))
            set(gca,'FontSize',12)
            xlabel('Domain Width (m)', 'FontSize', 10)
            ylabel('Domain Height (m)', 'FontSize', 10)
            if normalize
                set(gca, 'Clim', [-1 1])
            end% if normalize
            colorbar
            ch = colormap;
            title('\pi^\prime','fontsize',25)
            if cloud
                ch(64,1:3) = 1;
                hold on
                h = pcolor(0:gridwtm:domwtm+gridwtm,0:gridhtm:domht+gridhtm,ones(size(aqcnow)));
                alpha(h,(aqcnow))
%                al = get(gcf,'Alphamap');
%                al = 0:max(max(aqc))/63:max(max(aqc));
%                alphamap(al)
                hold off
            end %if cloud
            colormap(ch)
            shading flat
            
      drawnow
      frame = getframe(4);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if i == 1
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end %if-else
      
      close
      figure(4)
      set(4,'position',[200,200,800,700])
      
    end %for loop    
end %if animate
