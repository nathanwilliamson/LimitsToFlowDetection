function [Image,ImageMag]=get2Dfrom3D(Image3D,ImageMag3D,ind)

   Image=squeeze(Image3D(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6),:));
   ImageMag=squeeze(ImageMag3D(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6),:));
end