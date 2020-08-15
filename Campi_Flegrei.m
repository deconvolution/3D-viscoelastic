%% Dimensions
close all;
clear all;
%%
M=table2array(readtable('./modvPS.txt'));
M2=reshape(M,[61,96,116,6]);
M3=permute(M2,[3,2,1,4]);
M3=M3(1:50,1:30,1:20,:);
tit={'WE','SN','Al','vp','vs','vp/vs'};
figure('name','model');
for l2=1:6
    subplot(3,2,l2)
    pcolor3(M3(:,:,:,l2));
    xlabel('x');
    ylabel('y');
    zlabel('z');
    colorbar;
    title(tit(l2));
    set(gca,'zdir','reverse');
    shg;
end
%%
% https://www.nhc.noaa.gov/gccalc.shtml
dx0=790;
dy0=111;
dz0=100;
[nx0,ny0,nz0,~,~]=size(M3);
%% PML input
lp=20;
lpn=2;
Rc=10^-3;
%% compute C0
lambda=reshape(M3(:,:,:,4).^2,[1,1,nx0,ny0,nz0]);
mu=reshape(M3(:,:,:,5).^2,[1,1,nx0,ny0,nz0]);
C0=zeros(6,6,nx0,ny0,nz0);
C0(1,1,:,:,:)=lambda;
C0(2,2,:,:,:)=lambda;
C0(3,3,:,:,:)=lambda;
C0(4,4,:,:,:)=mu;
C0(5,5,:,:,:)=mu;
C0(6,6,:,:,:)=mu;
C0(1,2,:,:,:)=lambda-2*mu;
C0(1,3,:,:,:)=lambda-2*mu;
C0(2,3,:,:,:)=lambda-2*mu;
%% dimensions
dx=10;
dy=10;
dz=10;

nx2=round(nx0/dx*dx0);
ny2=round(ny0/dy*dy0);
nz2=round(nz0/dz*dz0);

C=zeros(6,6,nx2+2*lp+2,ny2+2*lp+2,nz2+2*lp+2);
[~,~,nx,ny,nz]=size(C);

for i=1:6
    for j=1:6
        C(i,j,lp+2:nx-lp-1,lp+2:ny-lp-1,lp+2:nz-lp-1)=reshape(imresize3(reshape(C0(i,j,:,:,:),[nx0,ny0,nz0]),[nx2,ny2,nz2]),[1,1,nx2,ny2,nz2]);
    end
end

C(:,:,1:lp+1,:,:)=repmat(C(:,:,lp+2,:,:),[1,1,lp+1,1,1]);
C(:,:,nx-lp:nx,:,:)=repmat(C(:,:,nx-lp-1,:,:),[1,1,lp+1,1,1]);

C(:,:,:,1:lp+1,:)=repmat(C(:,:,:,lp+2,:),[1,1,1,lp+1,1]);
C(:,:,:,ny-lp:ny,:)=repmat(C(:,:,:,ny-lp-1,:),[1,1,1,lp+1,1]);

C(:,:,:,:,1:lp+1)=repmat(C(:,:,:,:,lp+2),[1,1,1,1,lp+1]);
C(:,:,:,:,nz-lp:nz)=repmat(C(:,:,:,:,nz-lp-1),[1,1,1,1,lp+1]);
%%
figure
pcolor3(reshape(C(3,3,:,:,:),[nx,ny,nz]));
xlabel('x');
ylabel('y');
zlabel('z');
colorbar;
title('C33');
set(gca,'zdir','reverse');
shg;
%% time step
dt=10^-3; % [s]
nt=800; % Amount of time steps
ns=nt;
%% Define viscoelastic parameters
theta=0;
tau=5/8*theta;

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
Eta=zeros(6,6,nx,ny,nz);
for i=1:6
    for j=1:6
        Eta(i,j,:,:,:)=F0(i,j);
    end
end
rho=ones(nx,ny,nz)*1000;
%% Source and source signals
M=2.7;
sx=25;
sy=25;
sz=200;
sn=length(sx);
freq=10;
singles=rickerWave(freq,dt,ns,M);
srcx=zeros(nt,1);
srcy=srcx;
srcz=srcx;
srcx=1*singles;
srcy=1*singles;
srcz=1*singles;
%%
[ux,uy,uz]=solver(dt,dx,dy,dz,nt,nx,ny,nz,sx,sy,sz,srcx,srcy,srcz,lp,C,Eta,rho,lpn,Rc);
%%
lim2=.01*[min(ux(:)),max(ux(:))];
for l2=1:nt
    tt=ux(:,sy(1),:,l2);
    tt2=reshape(tt,[nx,nz]);
    figure(3)
    
    imagesc(tt2,lim2);
    title(num2str(l2));
    colorbar;
    shg;
end
%%
tt=permute(ux,[2,1,3,4]);
lim2=.01*[min(ux(:)),max(ux(:))];

col=[1,0,0;
    0,1,0;
    0,0,1;
    0,1,1;
    1,0,1;
    1,1,0;
    1,1,1;
    0,0,0];
rx=[10,20,30,40,51];
ry=[10,20,30,40,51];
rz=[1,1,1,1,1,1];
rt=dt:dt:dt*nt;
t3=zeros(length(rx),nt);
t4=t3;
for i=1:size(t3,1)
    t3(i,:)=real(reshape(ux(rx(i),ry(i),rz(i),:),[1,nt]));
    t4(i,:)=real(reshape(uz(rx(i),ry(i),rz(i),:),[1,nt]));
end

for l2=1:10:nt
    figure(3)
    set(gcf,'Visible','on');
    set(gcf,'position',[0,0,1500,600]);
    
    subplot(1,2,1)
    pcolor3(tt(:,:,:,l2));
    caxis(lim2);
    colorbar;
    set(gca,'zdir','reverse');
    xlabel(['x*' num2str(dx) '[m]']);
    ylabel(['y*' num2str(dx) '[m]']);
    zlabel(['z*' num2str(dx) '[m]']);
    title(['t=' num2str(dt*l2) 's']);
    
    for i=1:size(sx,2)
        hold on;
        ax2=plot3(sx(i),sy(i),sz(i),'o','color','red');
    end
    for i=1:size(rx,2)
        hold on;
        ax3=plot3(rx(i),ry(i),rz(i),'.','color',col(i,:),'markersize',20);
    end
    
    subplot(2,2,2)
    for i=1:length(rx)
        plot(rt(1:l2),t3(i,1:l2),'color',col(i,:));
        hold on;
    end
    xlabel('t [s]');
    ylabel('ux [m]');
    xlim([0,rt(end)]);
    ylim([min(t3(:)),max(t3(:))]);
    
    legend([ax2,ax3],'source','receiver','Location',[0.5,0.03,0.005,0.005],'orientation','horizontal');
    hold off;
    shg;
    % print(gcf,['.\marmousi2\' num2str(l) '.png'],'-dpng','-r100');
end
%% slice
slice_x=sx(1);
slice_y=sy(1);
slice_z=1;
tt=ux;
lim2=.01*[min(tt(:)),max(tt(:))];

for l2=[1:5:nt,nt]
    figure(4)
    subplot(2,2,1)
    imagesc(reshape(tt(slice_x,:,:,l2),[ny,nz]),lim2);
    set(gca,'zdir','reverse');
    xlabel(['z*' num2str(dz) '[m]']);
    ylabel(['y*' num2str(dy) '[m]']);
    title({['t=' num2str(dt*l2) 's'],['x=' num2str(dx*slice_x) 'm']});
    colorbar;
    
    subplot(2,2,2)
    imagesc(reshape(tt(:,slice_y,:,l2),[nx,nz]),lim2);
    set(gca,'zdir','reverse');
    xlabel(['z*' num2str(dz) '[m]']);
    ylabel(['x*' num2str(dx) '[m]']);
    title({['y=' num2str(dy*slice_y) 'm']});
    colorbar;
    
    subplot(2,2,3)
    imagesc(reshape(tt(:,:,slice_z,l2),[nx,ny]),lim2);
    set(gca,'zdir','reverse');
    xlabel(['y*' num2str(dy) '[m]']);
    ylabel(['x*' num2str(dx) '[m]']);
    title({['z=' num2str(dz*slice_z) 'm']});
    colorbar;
    hold off;
    
    shg;
end
%% slice
for l2=1:5:nt
    figure(99)
    slice(ux(:,:,:,l2),slice_x,slice_y,slice_z);
    set(gca,'zdir','reverse');
end