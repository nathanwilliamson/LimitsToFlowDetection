%% Initialization
clc
clear
%close all hidden
%% Path to subroutines

addpath('propagator_subroutines');

default = 1.0e+03 * [0.6170    0.3483    1.0187    0.7653];

set(0, 'DefaultFigurePosition', default);


%% Subtract no flow  y/n = 1/0?
subNoFlow=1;
%% load experimental data and parameters

%pathname='C:\Users\aawilliamsonnh\Documents\MATLAB\MATLAB_mac-@xi\Glymphatics\glymp0117.Ft1';
%pathname='E:\glymp_c';
%pathname='C:\Users\aawilliamsonnh\Documents\MATLAB\Glymphatics\glymph040119\glymph040119.SA1';
pathname='/Users/williamsonnh/Documents/MATLAB/Glymphatics/glymph040119/glymph040119.SA1';
%%%

%delta=3e-3; %% I was using this value previously

%delta=2.793e-3;  %% params file shows 3 ms but if taking into account the 0.8692 difference between the b that I calculate and the b that paravision calculates I would get this delta
delta=2.6893e-3;  %% params file shows 3 ms but if taking into account the 0.8692 difference between the b that I calculate and the b that paravision calculates I would get this delta
%delta=2.745e-3; %%
j=1;

%%% Choose your experiment by deleting the % at the begining of the line
%%%EXPERIMENT:NewData
im_ind{j}=[70]; Delta(j)=25e-3;  actualVdot(j)=0;  % flow rate in mul/min. +q
gradDir(j)=1; j=j+1; 
 im_ind{j}=[69]; Delta(j)=25e-3;  actualVdot(j)=0;  %-q
 gradDir(j)=-1; j=j+1;
im_ind{j}=[72]; Delta(j)=25e-3;  actualVdot(j)=3;  %flow rate in mul/min
gradDir(j)=1; j=j+1; 
im_ind{j}=[71]; Delta(j)=25e-3;  actualVdot(j)=3; %-q
gradDir(j)=-1; j=j+1;
im_ind{j}=[74]; Delta(j)=25e-3;  actualVdot(j)=10;  %flow rate in mul/min
gradDir(j)=1; j=j+1; 
im_ind{j}=[73]; Delta(j)=25e-3;  actualVdot(j)=10; %-q
gradDir(j)=-1; j=j+1;
im_ind{j}=[76]; Delta(j)=25e-3;  actualVdot(j)=30;  %flow rate in mul/min
gradDir(j)=1; j=j+1; 
im_ind{j}=[75]; Delta(j)=25e-3;  actualVdot(j)=30; %-q
gradDir(j)=-1; j=j+1;
im_ind{j}=[78]; Delta(j)=25e-3;  actualVdot(j)=100;  %flow rate in mul/min
gradDir(j)=1; j=j+1; 
im_ind{j}=[77]; Delta(j)=25e-3;  actualVdot(j)=100; %-q
gradDir(j)=-1; j=j+1;
im_ind{j}=[80]; Delta(j)=25e-3;  actualVdot(j)=20;  %flow rate in mul/min
gradDir(j)=1; j=j+1; 
im_ind{j}=[79]; Delta(j)=25e-3;  actualVdot(j)=20; %-q
gradDir(j)=-1; j=j+1;
im_ind{j}=[82]; Delta(j)=25e-3;  actualVdot(j)=6;  %flow rate in mul/min
gradDir(j)=1; j=j+1; 
im_ind{j}=[81]; Delta(j)=25e-3;  actualVdot(j)=6; %-q
gradDir(j)=-1; j=j+1;
im_ind{j}=[84]; Delta(j)=25e-3;  actualVdot(j)=1;  %flow rate in mul/min
gradDir(j)=1; j=j+1; 
im_ind{j}=[83]; Delta(j)=25e-3;  actualVdot(j)=1; %-q
gradDir(j)=-1; j=j+1;
im_ind{j}=[86]; Delta(j)=25e-3;  actualVdot(j)=60;  %flow rate in mul/min
gradDir(j)=1; j=j+1; 
im_ind{j}=[85]; Delta(j)=25e-3;  actualVdot(j)=60; %-q
gradDir(j)=-1; j=j+1;
im_ind{j}=[88]; Delta(j)=25e-3;  actualVdot(j)=2;  %flow rate in mul/min
gradDir(j)=1; j=j+1; 
im_ind{j}=[87]; Delta(j)=25e-3;  actualVdot(j)=2; %-q
gradDir(j)=-1; j=j+1;
im_ind{j}=[90]; Delta(j)=25e-3;  actualVdot(j)=0.6;  %flow rate in mul/min
gradDir(j)=1; j=j+1; 
im_ind{j}=[89]; Delta(j)=25e-3;  actualVdot(j)=0.6; %-q
gradDir(j)=-1; j=j+1;
im_ind{j}=[95]; Delta(j)=25e-3;  actualVdot(j)=300;  %flow rate in mul/min
gradDir(j)=1; j=j+1; 
im_ind{j}=[94]; Delta(j)=25e-3;  actualVdot(j)=300; %-q
gradDir(j)=-1; j=j+1;
im_ind{j}=[118]; Delta(j)=25e-3;  actualVdot(j)=200;  %flow rate in mul/min
gradDir(j)=1; j=j+1; 
im_ind{j}=[117]; Delta(j)=25e-3;  actualVdot(j)=200; %-q
gradDir(j)=-1; j=j+1;
tic
clear j

col=get(groot,'DefaultAxesColorOrder');
%col=distinguishable_colors(length(im_ind));
symbol={'o','v','^','d','>','s','d'};

for j=1:length(im_ind)
[Image{j},dims{j},fov{j},grad{j}]=load_data_3D_MagPhase(pathname,im_ind{j});
end


for j=1:length(im_ind)
params{j}=parameters(delta,Delta(j),gradDir(j)*grad{j});
end

for k=1:length(im_ind)/2
velocity.params{k}=params{k*2-1};
velocity.actualVdot(k)=actualVdot(k*2-1);
velocity.im_ind(k)=im_ind(k*2-1);
end

%%
%velocity.porosity=0.44; %'very lose random packing' (wikipedia)

%velocity.porosity=0.367; % measured for beadpack
velocity.porosity=0.4; % roughly found from conservation of mass based on sum of velocity through beadpack and velocity through bulk.

velocity.Rtube=2.5E-3; %m

for k=1:length(im_ind)/2
       if subNoFlow==1
            Image_plusq(:,:,:,:,1)= squeeze(Image{k*2-1}(:,:,:,:,1));
            Image_minusq(:,:,:,:,1)= squeeze(Image{k*2}(:,:,:,:,1));
            Image_plusq(:,:,:,:,2)= squeeze(Image{k*2-1}(:,:,:,:,2)-Image{1}(:,:,:,:,2));   %phasing each E(q) in each pixel by its no-flow phase.;
            Image_minusq(:,:,:,:,2)= squeeze(Image{k*2}(:,:,:,:,2)-Image{2}(:,:,:,:,2));
        elseif subNoFlow==0
            Image_plusq= Image{k*2-1};
            Image_minusq= Image{k*2};
       end

    [velocity.Image3D{k},velocity.ImageMag3D{k},velocity.ImageMag3D_minusq{k},velocity.ImagePhase3D_plusq{k},velocity.ImagePhase3D_minusq{k},velocity.ImageMag3D_diff{k},velocity.ImagePhase3D_diff{k}]=analyzeVelocimetry3D_MagPhase(Image_plusq,Image_minusq,params{1});
   
    velocity.beads.actualv{k}=actualVdot(k*2-1)*1E-9/60 /(pi*(velocity.Rtube)^2) /velocity.porosity; %m/s;
    velocity.bulk.actualv{k}=actualVdot(k*2-1)*1E-9/60 /(pi*(velocity.Rtube)^2) ;  
end

for k=1:length(velocity.Image3D)

    %%% saggital image
    indices=[7,18, 31,31, 20,43];
   [velocity.Image{k},velocity.ImageMag{k}]=get2Dfrom3D(velocity.Image3D{k},velocity.ImageMag3D{k},indices);
    %%% coronal image
%     velocity.Image{k}=squeeze(velocity.Image3D{k}(7:18,19:42,30,:));
%     velocity.ImageMag{k}=squeeze(velocity.ImageMag3D{k}(7:18,19:42,30,:));
    %%% axial images
    indices=[9,9, 19,42, 20,43];
   [velocity.beads.ImageSlice{k},velocity.beads.ImageSliceMag{k}]=get2Dfrom3D(velocity.Image3D{k},velocity.ImageMag3D{k},indices);

   indices=[15,15, 19,42, 20,43];
   [velocity.bulk.ImageSlice{k},velocity.bulk.ImageSliceMag{k}]=get2Dfrom3D(velocity.Image3D{k},velocity.ImageMag3D{k},indices);

   %%% velocity profiles
            %%%it looks like the fastest velocity is found on line 32, not
            %%%31, however this still doesn't explain the decreased maximum
            %%%velocity
    indices=[8,10, 31,31, 20,43];  % avoiding the voxels on either side of the interface between beads and bulk
   [velocity.beads.v{k},unused]=get2Dfrom3D(velocity.Image3D{k},velocity.ImageMag3D{k},indices);   
    velocity.beads.vprofileMean{k}=squeeze(nanmean(velocity.beads.v{k},1));
    velocity.beads.vprofileSD{k}=squeeze(nanstd(velocity.beads.v{k},0,1));
    
    indices=[14,16, 31,31, 20,43];  % avoiding the voxels on either side of the interface between beads and bulk    
    [velocity.bulk.v{k},unused]=get2Dfrom3D(velocity.Image3D{k},velocity.ImageMag3D{k},indices);   
    velocity.bulk.vprofileMean{k}=squeeze(nanmean(velocity.bulk.v{k},1));
    velocity.bulk.vprofileSD{k}=squeeze(nanstd(velocity.bulk.v{k},0,1));
end

for k=1:length(velocity.Image3D)

 
    for j=1:length(velocity.Image3D{k}(:,1,1,1))
        clear varray
        for l=1:length(velocity.Image3D{k}(1,1,1,:))
               clear vq;
            vq=squeeze(velocity.Image3D{k}(j,:,:,l)); %7:10
            varray(:,l)=vq(:);
        end
        velocity.varray{k}(j,:,:)=varray;
    end

    clear varray vq;
    for l=1:length(velocity.Image3D{k}(1,1,1,:))
        clear vq;
        vq=squeeze(velocity.Image3D{k}(8:10,:,:,l)); 
        varray(:,l)=vq(:);
    end
    velocity.beads.varray{k}=varray;
    velocity.beads.vMean{k}=nanmean(varray,1);
    velocity.beads.vSD{k}=nanstd(varray,0,1); 
    
%slice indices
%6=start of beadpack
%11=beadpack-bulk interface
%15=last good bulk flow slice before RF roundoff (but I used 18 in the
%saggital image

%%% two options for averaging and solving for Vdot.
%%% option 1:
pixelSize=15/64*1E-3; %(FOV/#pts)

        velocity.Vdot{k}(:,:)=NaN(length(velocity.Image3D{1}(:,1,1,1)),length(velocity.Image3D{1}(1,1,1,:)));
    for j=6:10
        velocity.Vdot{k}(j,:)=sum(sum(velocity.Image3D{k}(j,:,:,:),3,'omitnan'),2,'omitnan')*(pixelSize^2)/1E-9*60*velocity.porosity;
    end  
    velocity.beads.VdotMean{k}=mean(velocity.Vdot{k}(7:10,:),1);
    velocity.beads.VdotSD{k}=std(velocity.Vdot{k}(7:10,:),0,1);
    for j=12:18
        velocity.Vdot{k}(j,:)=sum(sum(velocity.Image3D{k}(j,:,:,:),3,'omitnan'),2,'omitnan')*(pixelSize^2)/1E-9*60;
    end
    velocity.bulk.VdotMean{k}=mean(velocity.Vdot{k}(13:16,:),1);
    velocity.bulk.VdotSD{k}=std(velocity.Vdot{k}(13:16,:),0,1);

%%% option 2 for averaging:
% velocity.Vdot{k}=squeeze(nanmean(velocity.varray{k},2))*pi*(velocity.Rtube)^2/1E-9*60;
% velocity.beads.VdotMean{k}=nanmean(velocity.Vdot{k}(8:10,:),1)*velocity.porosity;
% velocity.beads.VdotSD{k}=nanstd(velocity.Vdot{k}(8:10,:),0,1)*velocity.porosity;
% velocity.bulk.VdotMean{k}=nanmean(velocity.Vdot{k}(14:16,:),1);
% velocity.bulk.VdotSD{k}=nanstd(velocity.Vdot{k}(14:16,:),0,1);
    
%%% calculate SNR for beadpack and bulk regions
clear nq1 nq2 sq1 sq2 sq3 sq4 noiseq signalqbeads signalqbulk;
    for l=1:5
        %nq1=Image_plusq(:,1:15,:,l,1);
        nq1=Image_plusq(:,1:15,:,l,1);
        noiseq(1,l)=std(nq1(:));
        nq2=Image_minusq(:,1:15,:,l,1);
        noiseq(2,l)=std(nq2(:));
        sq1=squeeze(velocity.ImageMag3D{k}(8:10,:,:,l));
        signalqbeads(1,l)=mean(sq1(:),'omitnan');
        sq2=squeeze(velocity.ImageMag3D_minusq{k}(8:10,:,:,l));
        signalqbeads(2,l)=mean(sq2(:),'omitnan');
        sq3=squeeze(velocity.ImageMag3D{k}(14:16,:,:,l));
        signalqbulk(1,l)=mean(sq3(:),'omitnan');
        sq4=squeeze(velocity.ImageMag3D_minusq{k}(14:16,:,:,l));
        signalqbulk(2,l)=mean(sq4(:),'omitnan');
    end
    velocity.beads.signal{k}=signalqbeads;
    velocity.bulk.signal{k}=signalqbulk;
    velocity.beads.noise{k}=noiseq;
    velocity.bulk.noise{k}=noiseq;
    velocity.beads.SNR{k}=signalqbeads./noiseq;
    velocity.bulk.SNR{k}=signalqbulk./noiseq;
    
end

velocity.beads.dims=size(velocity.beads.v{1});
velocity.bulk.dims=size(velocity.bulk.v{1});


toc

%% save data
if subNoFlow==1
    velocity_subNoFlow=velocity;
    save('velocityData_subNoFlow','velocity_subNoFlow');
elseif subNoFlow==0
    save('velocityData','velocity');
end
