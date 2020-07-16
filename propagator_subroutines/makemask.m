function [mask_nan,mask_zero]=makemask(Image)
%thresh=max(abs(image))/10;
thresh=500000
dims=size(Image);
mask_nan=NaN(dims(1),dims(2),dims(3));
mask_zero=zeros(dims(1),dims(2),dims(3));
for j=1:dims(1)
    for k=1:dims(2)
        for l=1:dims(3)
            if abs(Image(j,k,l))>thresh
                mask_nan(j,k,l)=1;
                mask_zero(j,k,l)=1;
            end
        end
    end
end
end