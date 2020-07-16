function Enoflow=subtractNoFlow(image_roi)
[row_pix ,col_pix]=find(squeeze(image_roi(1,:,:,1))~=0);
k=1;
for i=1:length(row_pix)
    E=squeeze(image_roi(1,row_pix(k),col_pix(k),:));
    Enoflow(1,k,:)=E;
    k=k+1;
end
end