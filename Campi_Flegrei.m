%% Dimensions
close all;
clear all;
%%
M=table2array(readtable('./modvPS.txt'));
M2=reshape(M,[61,96,116,6]);
M3=permute(M2,[3,2,1,4]);
tit={'WE','SN','Al','vp','vs','vp/vs'};
figure
for l2=1:6
subplot(3,2,l2)
pcolor3(M3(:,:,:,l2));
xlabel('x');
ylabel('y');
zlabel('z');
colorbar;
title(tit(l2));
set(gca,'zdir','reverse');
end

% https://www.nhc.noaa.gov/gccalc.shtml
dx=790;
dy=111;
dz=100;
%% dimensions
dh=5; % Spacial grid step
dx=dh;
dy=dh;
dz=dh;
dt=10^-3; % [s]

nt=800; % Amount of time steps
ns=nt;

nx=71;
ny=51;
nz=41;
%% Define viscoelastic parameters
theta=0;
tau=5/8*theta;

lambda=10^9;
mu=10^9*.5;
C0=zeros(6,6);
C0(1,1)=lambda+2*mu;
C0(2,2)=C0(1,1);
C0(3,3)=C0(1,1);
C0(1,2)=lambda;
C0(1,3)=lambda;
C0(2,3)=lambda;
C0(4,4)=mu;
C0(5,5)=C0(4,4);
C0(6,6)=C0(4,4);

F0=zeros(6,6);
lambda2=lambda*theta+2/3*mu*(theta-tau);
mu2=mu*tau;
F0(1,1)=lambda2+2*mu2;
F0(2,2)=F0(1,1);
F0(3,3)=F0(1,1);
F0(1,2)=lambda2;
F0(1,3)=F0(1,2);
F0(2,3)=F0(1,2);
F0(4,4)=mu2;
F0(5,5)=F0(4,4);
F0(6,6)=F0(4,4);
%%
C=zeros(6,6,nx,ny,nz);
Eta=zeros(6,6,nx,ny,nz);
for i=1:6
    for j=1:6
        C(i,j,:,:,:)=C0(i,j);
        Eta(i,j,:,:,:)=F0(i,j);
    end
end
rho=ones(nx,ny,nz)*1000;
%% Source and source signals
M=2.7;
sx=23;
sy=25;
sz=3;
sn=length(sx);
freq=10;
singles=rickerWave(freq,dt,ns,M);
srcx=zeros(nt,1);
srcy=srcx;
srcz=srcx;
srcx=1*singles;
srcy=1*singles;
srcz=1*singles;
%% boundary condition input
lp=20;
lpn=2;
Rc=10^-3;
%%
[ux,uy,uz]=solver(dt,dx,dy,dz,nt,nx,ny,nz,sx,sy,sz,srcx,srcy,srcz,lp,C,Eta,rho,lpn,Rc);
%%
lim2=.1*[min(ux(:)),max(ux(:))];
for l2=1:nt
    tt=uy(sx,:,:,l2);
    tt2=reshape(tt,[ny,nz]);
    figure(3)
    
    imagesc(tt2,lim2);
    title(num2str(l2));
    colorbar;
    shg;
end
%%
tt=permute(ux,[2,1,3,4]);
lim2=.1*[min(ux(:)),max(ux(:))];
for l2=1:10:nt
    figure(3)
    pcolor3(tt(:,:,:,l2));
    set(gca,'zdir','reverse');
    caxis(lim2);
    xlabel(['x*' num2str(dx) '[m]']);
    ylabel(['y*' num2str(dx) '[m]']);
    zlabel(['z*' num2str(dx) '[m]']);
    title(num2str(l2));
    colorbar;
    shg;
end
%% ux and uz field
%{
col=[1,0,0;
    0,1,0;
    0,0,1;
    0,1,1;
    1,0,1;
    1,1,0;
    1,1,1;
    0,0,0];
rx=[50,100,150,200,250,280];
rz=[1,1,1,1,1,1];
rt=dt:dt:dt*nt;
t3=zeros(length(rx),nt);
t4=t3;
for i=1:size(t3,1)
    t3(i,:)=real(reshape(ux(rx(i),rz(i),:),[1,nt]));
    t4(i,:)=real(reshape(uz(rx(i),rz(i),:),[1,nt]));
end
lim2=.1*[min(min(real(ux(:)))),max(max(real(ux(:))))];
lim3=.1*[min(min(real(uz(:)))),max(max(real(uz(:))))];
for l=nt
    figure(4)
    set(gcf,'Visible','on');
    set(gcf,'position',[80,80,1300,600]);
    
    subplot(2,2,1)
    imshow(real(ux(:,:,l))',lim2);
    colorbar;
    hold on;
    contour(reshape(C(3,3,:,:),[nx,nz])','--','color','cyan');
    hold on;
    ax2=scatter(sx,sz,30,[1,0,0],'o');
    hold on;
    for i=1:length(rx)
        ax4=scatter(rx(i),rz(i),30,col(i,:),'filled');
        hold on;
    end
    ax3=plot([lp+1,lp+1],[1,nz-lp-1],'color','blue');
    hold on;
    ax3=plot([nx-lp-1,nx-lp-1],[1,nz-lp-1],'color','blue');
    hold on;
    hold on;
    ax3=plot([lp+1,nx-lp-1],[nz-lp-1,nz-lp-1],'color','blue');
    axis on;
    xlabel({['x*' num2str(dh) '[m]']});
    ylabel({['z*' num2str(dh) '[m]']});
    title({['t=' num2str(l*dt) 's'],['ux']});
    hold on;
    
    subplot(2,2,3)
    imshow(real(uz(:,:,l))',lim3);
    colorbar;
    hold on;
    contour(reshape(C(3,3,:,:),[nx,nz])','--','color','cyan');
    hold on;
    ax2=scatter(sx,sz,30,[1,0,0],'o');
    hold on;
    for i=1:length(rx)
        ax4=scatter(rx(i),rz(i),30,col(i,:),'filled');
        hold on;
    end
    ax3=plot([lp+1,lp+1],[1,nz-lp-1],'color','blue');
    hold on;
    ax3=plot([nx-lp-1,nx-lp-1],[1,nz-lp-1],'color','blue');
    hold on;
    hold on;
    ax3=plot([lp+1,nx-lp-1],[nz-lp-1,nz-lp-1],'color','blue');
    axis on;
    xlabel({['x*' num2str(dh) '[m]']});
    ylabel({['z*' num2str(dh) '[m]']});
    title({['uz']});
    hold on;
    
    subplot(2,2,2)
    for i=1:length(rx)
        plot(rt(1:l),t3(i,1:l),'color',col(i,:));
        hold on;
    end
    xlabel('t [s]');
    ylabel('ux [m]');
    xlim([0,rt(end)]);
    ylim([min(t3(:)),max(t3(:))]);
    
    subplot(2,2,4)
    for i=1:length(rx)
        plot(rt(1:l),t4(i,1:l),'color',col(i,:));
        hold on;
    end
    xlabel('t [s]');
    ylabel('uz [m]');
    xlim([0,rt(end)]);
    ylim([min(t4(:)),max(t4(:))]);
    
    legend([ax2,ax3,ax4],'source','PML boundary','receiver','Location',[0.5,0.03,0.005,0.005],'orientation','horizontal');
    shg;
    % print(gcf,['.\marmousi2\' num2str(l) '.png'],'-dpng','-r100');
end
hold off;
%% display

rx=50;
rz=1;
close;
t3=real(reshape(ux(rx,rz,:),[nt,1]));
t4=real(reshape(uz(rx,rz,:),[nt,1]));
rt=dt:dt:dt*nt;
ux2=ux/(max(real(uz(:))))*60;
uz2=uz/(max(real(uz(:))))*60;
%ux2(ux2<.5*max(ux2))=0;
%uz2(uz2<.5*max(uz2))=0;
for l=1:2:nt
    figure(20)
    set(gcf,'Visible','on');
    set(gcf,'position',[80,80,1000,600]);
    
    subplot(2,2,1)
    quiver(real(ux2(:,:,l)).',real(uz2(:,:,l)).',0);
    set(gca,'ydir','reverse');
    colorbar;
    xlim([0,nx]);
    ylim([-20,nz]);
    xlabel({['x*' num2str(dx) '[m]']});
    ylabel({['z*' num2str(dz) '[m]']});
    title({['t=' num2str(l*dt) 's'],['displacement vector [m]']});
    hold on;
    ax2=scatter(sx,sz,30,[1,0,0],'o');
    hold on;
    ax4=scatter(rx,rz,30,[0,1,0],'filled');
    hold on;
    ax3=plot([lp+1,lp+1],[1,nz-lp-1],'color','blue');
    hold on;
    ax3=plot([nx-lp-1,nx-lp-1],[1,nz-lp-1],'color','blue');
    hold on;
    ax3=plot([lp+1,nx-lp-1],[nz-lp-1,nz-lp-1],'color','blue');
    axis on;
    hold off;
    
    legend([ax2,ax4,ax3],'source','reiceiver','PML boundary','Location',[0.5,0.03,0.005,0.005],'orientation','horizontal');
    
    subplot(2,2,2)
    plot(rt(1:l),t3(1:l));
    xlabel('t [s]');
    ylabel('ux [m]');
    xlim([0,dt*nt]);
    ylim([min(t3),max(t3)]);
    
    subplot(2,2,4)
    plot(rt(1:l),t4(1:l));
    xlabel('t [s]');
    ylabel('uz [m]');
    xlim([0,dt*nt]);
    ylim([min(t4),max(t4)]);
    
    subplot(2,2,3)
    imagesc(reshape(C(3,3,:,:),[nx,nz])');
    colorbar;
    xlabel({['x*' num2str(dx) '[m]']});
    ylabel({['z*' num2str(dz) '[m]']});
    title('C33 [Pa]');
    
    shg;
    %print(gcf,['.\iso\' num2str(l) '.png'],'-dpng','-r100');
end
%}