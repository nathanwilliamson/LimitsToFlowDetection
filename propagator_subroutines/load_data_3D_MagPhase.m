function [Image,dims,fov,grad]=load_data_3D_MagPhase(pathname,im_ind)
%% Load PGSE image data, 2 slices, complex
%  close all
%  clear all
%  clc
%% define directory

 %im_ind=[9];

%reco_ind=[2,3];%!!!!!!assumes RECO files are 2-Real, 3-Imag!!!!!!!
reco_ind=[1,4];%trying to reconstruct the  image from magnitude and phase.


%pathname=uigetdir('/MATLAB/Glymphatics/glymp0117.Ft1','Pick Data Directory');
%pathname='C:\Users\aawilliamsonnh\Documents\MATLAB\MATLAB_mac-@xi\Glymphatics\glymp0117.Ft1';
%pathname='C:\Users\aawilliamsonnh\Documents\MATLAB\';
count=0;
for i=reco_ind
    count=count+1;
    %% read data 
    curdir=strcat(pathname,'/',num2str(im_ind),'/pdata/',num2str(i),'/2dseq'); %2dseq contains all the data (but in little endian) then we have to un-pack it (# grad pts, # slice, ect.)
    % curdir=strcat(pathname,'/',num2str(im_ind(1)),'/fid');
    
    version_file = [pathname, '/',num2str(im_ind(1)), '/acqp'];
    version_hand=fopen(version_file);
    version=fscanf(version_hand,'%c');
    version_place = strfind(version, '##$BYTORDA=');                                %finds if little or big endian
    tmp=textscan(version(version_place+11:version_place+17), '%s');
    version=tmp{1}{1};
    if strfind(version,'little')
        reco=fopen(curdir,'r','ieee-le');        % little endian
    else
        reco=fopen(curdir,'r');                  % big endian
    end                                                                             %cell refers to the entry bing within a (matlab) cell
    sliceNum_cell=find_multi_method(pathname, im_ind(1), 'PVM_SPackArrNSlices');    %Number of slices
    [dims_cell,num_dim]=find_multi_reco(pathname,im_ind(1),'RECO_size',i);          %dims_cell=number of pixels, eg 64x64
    [fov_cell,num_fov]=find_multi_reco(pathname,im_ind(1),'RECO_fov',i);          %fov_cell = field of view
    [slope_cell, num_grad]=find_multi_reco(pathname,im_ind(1),'RECO_map_slope',i);  %Slope = values which were used to normalize the signal intensity at each gradient point. We have to divide by this to un-normalize 
    %converting from strings within cells to numeric values
    dims=cellfun(@str2num,[dims_cell{:}]); %dimensions of image i.e. 64 x 64
    fov=cellfun(@str2num,[fov_cell{:}]); %field of view 
    slope=cellfun(@str2num,[slope_cell{:}]);
    sliceNum=cellfun(@str2num,[sliceNum_cell{:}]);
    
    num_grad=num_grad/sliceNum; 

  %  slope=slope2(1:2:end); %for 2 slices, it appears each entry in slope is repeated 2x, so this extracts only the each unique entry.
    
    
    A=fread(reco,'int16');
    tmp=length(A);
%% unpack data    
    i0=0;
    for m=1:num_grad
        for l=1:dims(3)
            for j=1:dims(2)
                A1(count,l,:,j,m) = A(i0+1:i0+dims(1)); %Unpacking the data. A1 will be [real/imag, #slice, #read dimension, #phase dimension, #grads]
                i0=i0+dims(1);
            end
        end
    end
    
    for j=1:num_grad
        Image_split(count,:,:,:,j)=A1(count,:,:,:,j)/slope(j); %un-normalizing the data

    end
    fclose(reco);
    fclose(version_hand);
end
Image_mag=squeeze(Image_split(1,:,:,:,:));
Image_phase=squeeze(Image_split(2,:,:,:,:));
Image(:,:,:,:,1)=Image_mag;
Image(:,:,:,:,2)=Image_phase;

%Image(=Image_mag.*(cos(Image_phase)+1i*(sin(Image_phase)));
%Image=complex(Image_real,Image_img); %storing the complex 3D image data


%% gradient values

[grad_cell, ~]=find_multi_method(pathname,im_ind(1),'PVM_DwGradAmp');
grad=cellfun(@str2num,[grad_cell{:}]);

