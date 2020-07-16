function MSDscaling=analyzeMSD(transSlice,Delta)
MSD_mat=NaN(length(transSlice),size(transSlice{1}.image.MSD,1),size(transSlice{1}.image.MSD,2));
  
for j=1:length(transSlice)
    %MSD_mat(j,:,:)=transSlice{j}.image.MSD;
    MSD_ROI(j)=transSlice{j}.ROI.MSD;
end
log10MSD=log10(MSD_ROI);
log10Delta=log10(Delta);
y=log(MSD_ROI);
x=log(Delta);

%bounds and constraints for parameters
Aeq=[];
Beq=[];
lb=[-inf 0];
ub=[0 2];


%fit model 
options=optimoptions('fmincon');
options.Algorithm='sqp';
options.Display='off';
options.MaxFunctionEvaluations = 10000; %modern settings;
options.OptimalityTolerance = 1e-8;
options.ConstraintTolerance = 1e-8;


ss = inf;
for k=1:100
    %generate initial values of parameters
    param_guess(1)=min(log10MSD)*1/rand(1);
    param_guess(2)=rand(1)*2;
    %fit model
    paramhat_=fmincon(@(param)sumofsquares(log10Delta, log10MSD,param),param_guess,[],[],[],[],lb,ub,[],options);
    
    % compute residual sum of squares
    logMSDmodel_ = paramhat_(1)+paramhat_(2)*log10Delta;
    ss_ = sum( (log10MSD - logMSDmodel_).^2);
    
    if ss_ < ss
        paramhat = paramhat_;
        ss = ss_;
        logMSDmodel=logMSDmodel_;
    end
end

% also calculate coefficients based on
% http://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html 
%also exlplained at http://mathworld.wolfram.com/LeastSquaresFitting.html
MSDscaling.alphaMAT=(length(transSlice)*sum(x.*y)-sum(x)*sum(y))/(length(transSlice)*sum(x.^2)-sum(x)^2);


MSDscaling.C=paramhat(1);
MSDscaling.alpha=paramhat(2);
MSDscaling.logDelta=log10Delta;
MSDscaling.logMSD=log10MSD;
MSDscaling.logMSDmodel=logMSDmodel;
MSDscaling.ss=ss;

        
    
