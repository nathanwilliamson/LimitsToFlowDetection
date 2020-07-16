function  [velocityImage,MagImage]=analyzeVelocimetry3D(Image_plusq,Image_minusq,params)
if length(size(Image_plusq))==3
    thresh=1E5;
    MagImage=abs(Image_plusq);
    phaseDiffImage=unwrap(angle(Image_plusq),[],2)-unwrap(angle(Image_minusq),[],2);
    %phaseDiffImage_=angle(Image_plusq_)-angle(Image_minusq_); 

    for h=1:length(phaseDiffImage(:,1,1))
        for i=1:length(phaseDiffImage(1,:,1))
            for j=1:length(phaseDiffImage(1,1,:))
                if abs(Image_plusq(h,i,j))<thresh
                    phaseDiffImage(h,i,j)=0;
                    MagImage(h,i,j)=0;
                elseif abs(Image_minusq(h,i,j))<thresh
                    phaseDiffImage(h,i,j)=0;
               % elseif phaseDiffImage(h,i,j)>0
                %    phaseDiffImage(h,i,j)=0;
                end
            end
        end
    end
    for i=1:length(phaseDiffImage(1,1,:))
        velocityImage(:,:,i)=phaseDiffImage(:,:,i)/(params.Delta*params.q_rad(i)*2);
    end
elseif length(size(Image_plusq))==4
    thresh=1E5;
    MagImage=abs(Image_plusq);
    phaseDiffImage=unwrap(angle(Image_plusq),[],2)-unwrap(angle(Image_minusq),[],2);
    %phaseDiffImage_=angle(Image_plusq_)-angle(Image_minusq_); 

    for i=1:length(phaseDiffImage(1,1,1,:))
        velocityImage(:,:,:,i)=phaseDiffImage(:,:,:,i)/(params.Delta*params.q_rad(i)*2);
    end

end
end

