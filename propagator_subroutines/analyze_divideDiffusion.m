function trans=analyze(image_roi,params)
E_noflow=subtractNoFlow(image_roi(1,:,:,:)); %the no flow E(q), for subtraction

N=params.N; qinv=params.qinv; Delta=params.Delta; q=params.q; dqinv=params.dqinv; b_=params.b;

q_rad=q*2*pi; % needs to be in radians/m for the velocity from phase shift calculation

[row_pix ,col_pix]=find(squeeze(image_roi(1,:,:,1))~=0);

trans.image.P=NaN(size(image_roi,2),size(image_roi,3),N);
trans.image.v=NaN(size(image_roi,2),size(image_roi,3));
trans.image.vavg_T=NaN(size(image_roi,2),size(image_roi,3));
trans.image.D_fit=NaN(size(image_roi,2),size(image_roi,3));
trans.image.MSD=NaN(size(image_roi,2),size(image_roi,3));

trans.Esum=0;

k=1;
for i=1:length(row_pix)
        E = squeeze(image_roi(2,row_pix(i),col_pix(i),:));
        E = E.*exp(-1i*angle(squeeze(E_noflow(1,k,:))));   %phasing each E(q) in each pixel by its no-flow phase.
        E = E.*exp(-1i*angle(E(1)));
        %E=E./abs(E); %%new line, trying to divide out diffusion
        
        E_ft = [conj(flipud(E(2:end)));E(1:end-1)];   % we do the conjugate of E for negative q because the imaginary signal has opposite sign after reflection.
        P = fftshift(fft(fftshift(E_ft),N));
        P = P/sum(dqinv*P); %note that this is now normalized.
        
        iTmax = find(angle(E) > (pi/2),1); % the index for the cutoff of finding velocity from a linear fit of phase vs Theta. Can't use all the data or you'll hit Nyquist's thm limit.
        T = angle(E(1:iTmax));
        v = sum((qinv-dqinv).*real(P))/sum(real(P))/Delta;  %simply the first moment of the propagator
        v_T = T'*q_rad(1:iTmax)/(q_rad(1:iTmax)'*q_rad(1:iTmax)) / Delta; %fit a line, since phase at q=0 is zero, we can use matrix form
        if isnan(v_T)==1
            v_T=0;
        end
        logE=log(abs(E)/abs(E(1)))';% by normalizing by E(1) we assume that E(1) is equal to E(q=0).
        Nfit=find(logE < -0.4,1); % the index for the cuttof for single exp. fit of D.
        b=b_(1:Nfit);
        vector_E=logE(1:Nfit)';
        D_fit=-vector_E'*b/(b'*b); % Vector calculus fit of the slope. Equation can be seen on https://en.wikibooks.org/wiki/Linear_Algebra/Topic:_Line_of_Best_Fit
        MSD = 2*D_fit*Delta;
        
        trans.image.I_0(row_pix(i),col_pix(i),:)=abs(E(1));
        trans.image.P(row_pix(i),col_pix(i),:)=P;
        trans.image.v(row_pix(i),col_pix(i))= v;
        trans.image.v_T(row_pix(i),col_pix(i)) = v_T;
        trans.image.D_fit(row_pix(i),col_pix(i)) = D_fit;
        trans.image.MSD(row_pix(i),col_pix(i)) = MSD;
        
        trans.vec.E(k,:)=E;       % note we are taking the phased E
        trans.vec.P(k,:)=P;       %note that these are normalized so taking an average here will me weighted by roi, not by spin density.
        
        k=k+1;
end
E=sum(trans.vec.E,1).'; %note the nonconjugate transpose .'
E_temp=E./abs(E).*sign(E);;
E_ft = [conj(flipud(E_temp(2:end)));E_temp(1:end-1)];

P = fftshift(fft(fftshift(E_ft),N)); 
P = P./sum(dqinv*P); 

iTmax = find(angle(E) > (pi/2),1); 
T = angle(E(1:iTmax));
v = qinv'*real(P)/sum(dqinv*real(P))/Delta;
v_T = T'*q_rad(1:iTmax)/(q_rad(1:iTmax)'*q_rad(1:iTmax)) / Delta;
if isnan(v_T)==1
    v_T=0;
end
logE=log(abs(E)/abs(E(1)))';% by normalizing by E(1) we assume that E(1) is equal to E(q=0).
Nfit=find(logE < -0.4,1);
b=b_(1:Nfit);
vector_E=logE(1:Nfit)';
D_fit=-vector_E'*b/(b'*b);
MSD=2*D_fit*Delta;

ROI.E=E;
ROI.vmax=max(max(trans.image.v));
ROI.vmax_T=max(max(trans.image.v_T));
ROI.P=P;
ROI.vavg=v;
ROI.vavg_T=v_T;
ROI.D_fit=D_fit;
ROI.MSD=MSD;
ROI.meanT=mean(angle(trans.vec.E),1); %T stands for theta
ROI.stdT=std(angle(trans.vec.E),1);

trans.ROI=ROI;
        
end