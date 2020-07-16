function  [velocityImage,magImage_plusq,magImage_minusq,phaseImage_plusq,phaseImage_minusq,magDiffImage,phaseDiffImage]=analyzeVelocimetry3D_MagPhase(Image_plusq,Image_minusq,params)
    thresh= 5E5;
if length(size(Image_plusq))==4
    disp('image is the wrong dimension');

elseif length(size(Image_plusq))==5

    magImage_plusq=squeeze(Image_plusq(:,:,:,:,1));
    magImage_minusq=squeeze(Image_minusq(:,:,:,:,1));
    magDiffImage=magImage_plusq-magImage_minusq;
    %phaseDiffImage=unwrap(squeeze(Image_plusq(:,:,:,2)),[],2)-unwrap(squeeze(Image_minusq(:,:,:,2)),[],2); 
    phaseImage_plusq=squeeze(Image_plusq(:,:,:,:,2));
    phaseImage_minusq=squeeze(Image_minusq(:,:,:,:,2));    
   % phaseDiffImage=phaseImage_plusq-phaseImage_minusq;
    phaseDiffImage(:,:,:,1)=phaseImage_plusq(:,:,:,2)-phaseImage_plusq(:,:,:,1);
    phaseDiffImage(:,:,:,2)=phaseImage_plusq(:,:,:,3)-phaseImage_plusq(:,:,:,2);
    phaseDiffImage(:,:,:,3)=phaseImage_plusq(:,:,:,4)-phaseImage_plusq(:,:,:,3);
    phaseDiffImage(:,:,:,4)=phaseImage_plusq(:,:,:,5)-phaseImage_plusq(:,:,:,4);
    phaseDiffImage(:,:,:,5)=phaseImage_plusq(:,:,:,5)-phaseImage_plusq(:,:,:,3);

   
    %%% unwraping based on the median of nearby cells  
    p=1;
    for j=1:length(phaseDiffImage(1,1,1,:))
        %pad the image so that there are cells bordering all image cells
        %use NaN so that image doesn't display those values
        phaseImgPadded(:,:,:,j) = NaN(size(squeeze(phaseDiffImage(:,:,:,1)))+2*p); phaseImgPadded(p+1:end-p,p+1:end-p,p+1:end-p,j) = squeeze(phaseDiffImage(:,:,:,j));
        magImgPadded(:,:,:,j) = NaN(size(squeeze(magImage_plusq(:,:,:,1)))+2*p); magImgPadded(p+1:end-p,p+1:end-p,p+1:end-p,j) = squeeze(magImage_plusq(:,:,:,j));
        
        nanmask(:,:,:,j)=ones(size(squeeze(magImage_plusq(:,:,:,1))));
    end
    for h=1:length(phaseDiffImage(:,1,1,1))
        for i=1:length(phaseDiffImage(1,:,1,1))
            for j=1:length(phaseDiffImage(1,1,:,1))
                 h_=h+1;i_=i+1; j_=j+1;
                 if magImgPadded(h_,i_,j_,1)<thresh
                    phaseImgPadded(h_,i_,j_,:)=NaN;
                    magImgPadded(h_,i_,j_,:)=NaN;
                    
                    nanmask(h,i,j,:)=NaN;
                 end
            end
        end
    end
    for h=1:length(phaseDiffImage(:,1,1,1))
        for i=1:length(phaseDiffImage(1,:,1,1))
            for j=1:length(phaseDiffImage(1,1,:,1))
                 h_=h+1;i_=i+1; j_=j+1;
                 A=magImgPadded(h_-1:h_+1,i_-1:i_+1,j_-1:j_+1,1);
                 num_no_nans = sum(  ~isnan(A(:)));
                 if num_no_nans<3
                    phaseImgPadded(h_,i_,j_,:)=NaN;
                    magImgPadded(h_,i_,j_,:)=NaN;
                    
                    nanmask(h,i,j,:)=NaN;
                 end
            end
        end
    end
    
    
%  ci=0.7;
%     for h=1:length(phaseDiffImage(:,1,1,1))
%         for i=1:length(phaseDiffImage(1,:,1,1))
%             for j=1:length(phaseDiffImage(1,1,:,1))
%                  h_=h+1;i_=i+1; j_=j+1;
%                  if MagImage(h,i,j,1)>thresh
%                     for k=1:length(phaseDiffImage(1,1,1,:))
%                         %resolving wrapping by phase differeces between
%                         %neiboring cells, looking upwind.
%                         if imgPadded(h_,i_,j_,k) < actualv*(params.Delta*params.q_rad(k)*2)-ci*2*pi 
%                             imgPadded(h_,i_,j_,k)=imgPadded(h_,i_,j_,k)+2*pi;
%                         end                        
%                     end
%                 end
%             end
%         end
%     end

 %%% each unwrap step eats in on the wrapped areas a little more
 ci=0.65;
 for unwrap=1:10    
    for h=1:length(phaseDiffImage(:,1,1,1))
        for i=1:length(phaseDiffImage(1,:,1,1))
            for j=1:length(phaseDiffImage(1,1,:,1))
                 h_=h+1;i_=i+1; j_=j+1;
                 if magImage_plusq(h,i,j,1)>thresh
                    for k=1:length(phaseDiffImage(1,1,1,:))
                        
%                        finding cells that are more that pi off and unwrapping
%                        them.                
                            A=phaseImgPadded(h_-1:h_+1,i_-1:i_+1,j_-1:j_+1,k);
                            med=median(A(:),'omitnan');

                            if phaseImgPadded(h_,i_,j_,k)<med-ci*pi 
                                phaseImgPadded(h_,i_,j_,k)=phaseImgPadded(h_,i_,j_,k)+pi;
                            end
                        
                    end
                end
            end
        end
    end

 end
 %%% data is mostly wrapped backwards, not forwards, so doing forward
 %%% unwrapping secons and with only one step
 for unwrap=1:2  
    for h=1:length(phaseDiffImage(:,1,1,1))
        for i=1:length(phaseDiffImage(1,:,1,1))
            for j=1:length(phaseDiffImage(1,1,:,1))
                 h_=h+1;i_=i+1; j_=j+1;
                 if magImage_plusq(h,i,j,1)>thresh
                    for k=1:length(phaseDiffImage(1,1,1,:))
                        
%                        finding cells that are more that pi off and unwrapping
%                        them.                
                            A=phaseImgPadded(h_-1:h_+1,i_-1:i_+1,j_-1:j_+1,k);
                            med=median(A(:),'omitnan');

                    %        if imgPadded(h_,i_,j_,k)<med-ci*pi 
                    %            imgPadded(h_,i_,j_,k)=imgPadded(h_,i_,j_,k)+pi;
                            if phaseImgPadded(h_,i_,j_,k)>med+ci*pi
                                phaseImgPadded(h_,i_,j_,k)=phaseImgPadded(h_,i_,j_,k)-pi;
                            end
                        
                    end
                end
            end
        end
    end

 end 
 

    clear phaseDiffImage MagImage;
    phaseDiffImage=phaseImgPadded(1+p:end-p, 1+p:end-p, 1+p:end-p, :);
    magImage_plusq=magImgPadded(1+p:end-p, 1+p:end-p, 1+p:end-p, :);    
    
    magImage_minusq=magImage_minusq.*nanmask;
    phaseImage_plusq=phaseImage_plusq.*nanmask;
    phaseImage_minusq=phaseImage_minusq.*nanmask;
    magDiffImage=magDiffImage.*nanmask;
    
     for i=1:length(phaseDiffImage(1,1,1,:))-1
            velocityImage(:,:,:,i)=phaseDiffImage(:,:,:,i)/(params.Delta*(params.q_rad(i+1)-params.q_rad(i)));
  %         velocityImage=phaseDiffImage/(2*pi);
     end
            velocityImage(:,:,:,5)=phaseDiffImage(:,:,:,5)/(params.Delta*(params.q_rad(5)-params.q_rad(3)));


end
end

