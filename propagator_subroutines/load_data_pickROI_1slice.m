%function [Image_comp_s1,Image_roi_s1,bw_mat,dims,fov,mean_roi_s1,grad]=load_data_pickROI_1slice(pathname,im_ind)
%% Load PGSE image data, 2 slices, complex
%  close all
%  clear all
%  clc
%% define directory

 im_ind=[9];

reco_ind=[2,3];%!!!!!!assumes RECO files are 2-Real, 3-Imag!!!!!!!

%pathname=uigetdir('/MATLAB/Glymphatics/glymp0117.Ft1','Pick Data Directory');
%pathname='C:\Users\aawilliamsonnh\Documents\MATLAB\MATLAB_mac-@xi\Glymphatics\glymp0117.Ft1';
pathname='C:\Users\aawilliamsonnh\Documents\MATLAB\';

for i=reco_ind
    %% read data 
    curdir=strcat(pathname,'/',num2str(im_ind(1)),'/pdata/',num2str(i),'/2dseq'); %2dseq contains all the data (but in little endian) then we have to un-pack it (# grad pts, # slice, ect.)
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
        for l=1:sliceNum
            for j=1:dims(2)
                A1(i-reco_ind(1)+1,l,:,j,m) = A(i0+1:i0+dims(1)); %Unpacking the data. A1 will be [real/imag, #slice, #read dimension, #phase dimension, #grads]
                i0=i0+dims(1);
            end
        end
    end
%         A_s1=A(1:tmp/2);
%         A_s2=A(tmp/2+1:end);
%         
% %         A_s1=A(1:2:end);
% %         A_s2=A(2:2:end);
%         
%         A1_s1(i-reco_ind(1)+1,:,:,:)=reshape(A_s1,[dims(1),dims(2),num_grad]);
%         A1_s2(i-reco_ind(1)+1,:,:,:)=reshape(A_s2,[dims(1),dims(2),num_grad]);
        
        A1_s1(i-reco_ind(1)+1,:,:,:)=A1(i-reco_ind(1)+1,1,:,:,:); %s1 and s2 are slice 1 and 2


    
    for j=1:num_grad
        Image_s1(i-reco_ind(1)+1,:,:,j)=A1_s1(i-reco_ind(1)+1,:,:,j)/slope(j); %un-normalizing the data

    end
    
    fclose(reco);
    fclose(version_hand);
end

% Image_comp=complex(A1(1,:,:,:,:),A1(2,:,:,:,:));
% Image_mag=abs(squeeze(Image_comp));
% 
% figure
% imagesc(squeeze(Image_mag(2,:,:,1)))
% colormap(gray)

%% return relevant matrices

Image_real_s1=squeeze((Image_s1(1,:,:,:)));  %now, separating the data into real and imaginary individual matrices and then combining again
Image_img_s1=squeeze((Image_s1(2,:,:,:)));
Image_comp_s1(1,:,:,:)=complex(Image_real_s1,Image_img_s1); %combining into complex matrices. Note the first dimension of the matrix is used as the experiment dimension. Below we will (potentially) add to this if there are more experiments
Image_mag_s1=abs(squeeze(Image_comp_s1));




clear A     %delete intermediaties
clear A_s1
clear A_s2
clear A1
clear A1_s1
clear A1_s2
clear reco
clear Image_s1
clear Image_s2


%% select ROI 

figure
imagesc(Image_mag_s1(:,:,1)) 
colormap(gray)

h = imellipse;  %draw an ellipse, for the ROI
position = wait(h);
setColor(h,'g')

bw_mat = createMask(h); %the ROI is used to creat a mask, which will be used to pull out only the data in the ROI in the following steps

% 
% [bw_mat , xi, yi] =roipoly;
% p = patch(xi,yi, 'g','FaceColor','none');
% set(p,'EdgeColor','g','LineWidth',2)

%% read and unpack data from other experiments in im_ind (skips the one which was previously read above

for k=2:length(im_ind)
    for i=reco_ind
        
        curdir=strcat(pathname,'/',num2str(im_ind(k)),'/pdata/',num2str(i),'/2dseq');
        % curdir=strcat(pathname,'/',num2str(im_ind(1)),'/fid');
        
        version_file = [pathname, '/',num2str(im_ind(k)), '/acqp'];
        version_hand=fopen(version_file);
        version=fscanf(version_hand,'%c');
        version_place = strfind(version, '##$BYTORDA=');
        tmp=textscan(version(version_place+11:version_place+17), '%s');
        version=tmp{1}{1};
        if contains(version,'little')
            reco=fopen(curdir,'r','ieee-le');        % little endian
        else
            reco=fopen(curdir,'r');                  % big endian
        end
        
    sliceNum_cell=find_multi_method(pathname, im_ind(k), 'PVM_SPackArrNSlices');
    [dims_cell,num_dim]=find_multi_reco(pathname,im_ind(k),'RECO_size',i);
    [slope_cell, num_grad]=find_multi_reco(pathname,im_ind(k),'RECO_map_slope',i);
    
    dims=cellfun(@str2num,[dims_cell{:}]);
    slope=cellfun(@str2num,[slope_cell{:}]);
    sliceNum=cellfun(@str2num,[sliceNum_cell{:}]);
    
    %slope=slope2(1:2:end);
    num_grad=num_grad/sliceNum;
    
    A=fread(reco,'int16');
    
    i0=0;
    for m=1:num_grad
        for l=1:sliceNum
            for j=1:dims(2)
                A1(i-reco_ind(1)+1,l,:,j,m) = A(i0+1:i0+dims(1));
                i0=i0+dims(1);
            end
        end
    end
    
      
        A1_s1(i-reco_ind(1)+1,:,:,:)=A1(i-reco_ind(1)+1,1,:,:,:);
% 
%         tmp=length(A);
%         A_s1=A(1:tmp/2);
%         A_s2=A(tmp/2+1:end);
%         
%         A1_s1(i-reco_ind(1)+1,:,:,:)=reshape(A_s1,[dims(1),dims(2),num_grad]);
%         A1_s2(i-reco_ind(1)+1,:,:,:)=reshape(A_s2,[dims(1),dims(2),num_grad]);

    
    for j=1:num_grad
        Image_s1(i-reco_ind(1)+1,:,:,j)=A1_s1(i-reco_ind(1)+1,:,:,j)/slope(j);
    end
    
        
        fclose(reco);
        fclose(version_hand);
    end
    

    
Image_real_s1=squeeze((Image_s1(1,:,:,:)));
Image_img_s1=squeeze((Image_s1(2,:,:,:)));
Image_comp_s1(k,:,:,:)=complex(Image_real_s1,Image_img_s1);



    
    
clear A
clear A_s1
clear A_s2
clear A1
clear A1_s1
clear A1_s2
clear reco
clear Image_s1
clear Image_s2
    
end

%% store complex ROI data for every experiment 

for k=1:length(im_ind)
    for i=1:num_grad
        Image_roi_s1(k,:,:,i)=squeeze(Image_comp_s1(k,:,:,i)).*bw_mat;
        tmp1=squeeze(Image_roi_s1(k,:,:,i));
        roi_s1{k,i}=tmp1(tmp1~=0);
        

    end
end

% Image_roi=squeeze(Image_roi);


for k=1:length(im_ind)
    for i=1:num_grad
        mean_roi_s1(i)=mean(mean(roi_s1{k,i}));
        Image_cut_s1(k,i,:)=roi_s1{k,i};
        

    end
end


%% gradient values

[grad_cell, ~]=find_multi_method(pathname,im_ind(1),'PVM_DwGradAmp');
grad=cellfun(@str2num,[grad_cell{:}]);

