%% plot velocimetry acuracy
close all
clear all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%% Subtract no flow  y/n = 1/0?
load('velocityData.mat')
[B,I]=sort(velocity.actualVdot);

D=1.2E-9;
Delta=25E-3;
min=sqrt(2*D/Delta)*1E6;
SNR=70;
minarray=[1/(sqrt(2)*SNR),pi/2]*min;
qrad=velocity.params{1}.q_rad;
qstar=qrad*sqrt(2*D*Delta);
%%%old
% D=1.2E-9;
% Delta=25E-3;
% min=sqrt(2*D/Delta)*1E6;
% SNR=70;
% minarray=[exp(1/2)/SNR,1]*min;
% qrad=velocity.params{1}.q_rad;
% qstar=qrad*sqrt(2*D*Delta);
%%%%q^* for beadpack (D_0=1.2E-9)
% % %     0.416
% % %     0.665
% % %     0.998
% % %     1.248
% % %     1.497
%%%%q^* for bulk region (D_0=1.89E-9)
% % %     0.521
% % %     0.835
% % %     1.253
% % %     1.566
% % %     1.879
%qtext={'$q^*=0.49$', '$q^*=0.78$', '$q^*=1.18$', '$q^*=1.46$', '$q^*=1.75$'} ;
qtext={'$q_1$', '$q_2$', '$q_3$', '$q_4$', '$q_5$'} ; %


for j=1:2
    subNoFlow=j-1;%%
    if subNoFlow==1
        load('velocityData_subNoFlow.mat')
        velocity=velocity_subNoFlow;
    elseif subNoFlow==0
        load('velocityData.mat')
    end
  
    for k=1:length(I)
        vMean(j,k,:)=velocity.beads.vMean{I(k)}(:)*1E6;
        vSD(j,k,:)=velocity.beads.vSD{I(k)}(:)*1E6;
        actualv(k)=velocity.beads.actualv{I(k)}(1)*1E6;
        vError(j,k,:)=velocity.beads.vMean{I(k)}(:)/velocity.beads.actualv{I(k)}(1);

    end
end
%% plot accuracy
col=get(groot,'DefaultAxesColorOrder');
COLORSpatch=[255 255 234; 234 255 249; 244 240 255]/255;

fig                 = figure();
fig.Units           = 'centimeters';
fig.PaperUnits      = 'centimeters';
fig.Position        = [0 0 8 20];
fig.PaperPosition   = fig.Position;

FontName            = 'helvetica';
FontSize            = 7;
FontWeight          = 'normal';

for q=1:5

    %h=subplot(5,3,(q-1)*3+1);
    h=subplot(5,1,q);

        hold on
    h.YLim=[0 660];
    h.XLim=[0 660];
    h.YLabel.String='measured $v\ \mathrm{[\mu m /s]}$';
    h.XLabel.String='';
    h.FontName          = FontName;
    h.FontSize          = FontSize;
    h.FontWeight        = FontWeight;
     h.Units      = 'centimeters';

    h.Position(4)          = 3;
    pos_(q,:)= h.Position;
  %  h.Position(3)          = 3.5;

    h.YTick             = [0 200 400 600];
   h.XTick             = [0 200 400 600];
   h.XTickLabel= {'' '' '' ''};
% h.YLabel.Units      = 'centimeters';
% h.YLabel.Position(1) = -0.8;

h.YMinorTick = 'on';
h.XMinorTick = 'on';
h.TickLength = [.02 .02] ;
h.TickDir    = 'out';
h.Box               = 'on';

ypatch_=h.YLim+[5 -5];
ypatch=[ypatch_ fliplr(ypatch_)];

hf = patch([2 2 minarray(1) minarray(1)],ypatch ,COLORSpatch(1,:));
hf.LineStyle='none';
hf = patch([minarray(1) minarray(1) minarray(2) minarray(2)] ,ypatch, COLORSpatch(2,:));
hf.LineStyle='none';
hf = patch([minarray(2) minarray(2) h.XLim(2)-2 h.XLim(2)-2],ypatch, COLORSpatch(3,:));
hf.LineStyle='none';

    hl=line([-600,1000],[-600,1000]);
    hl.LineStyle='-';
    hl.Color= [0 0 0];
    hl=errorbar(actualv,vMean(1,:,q),vSD(1,:,q));
    hl.MarkerSize=6;
        hl.Marker='d';
        hl.LineStyle='none';
        hl.Color= col(6,:)/3;
        hl.MarkerFaceColor=col(6,:);
% 
%     hl=errorbar(actualv,vMean(2,:,q),vSD(2,:,q));
%     hl.Marker='d';
%         hl.MarkerSize=3;
%     hl.LineStyle='none';
%     hl.Color= col(3,:)/2;
%     hl.MarkerFaceColor=col(3,:);

    %annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', 'q')
    text(500,150,qtext{q},'FontSize',9)
end
    h.XLabel.String='actual $v\ \mathrm{[\mu m /s]}$';
       h.XTickLabel= {'0' '200' '400' '600'};

    for q=1:5

        h=axes();

        hold on
        h.YLim=[-5 70];
        h.XLim=[-5 70];
        h.YLabel.String='';
        h.XLabel.String='';
        h.FontName          = FontName;
        h.FontSize          = FontSize;
        h.FontWeight        = FontWeight;
        h.Units      = 'centimeters';

        h.Position         = pos_(q,:)+[ .5 1.2 -3 -1];

        h.YTick             = [0 50];
        h.XTick             = [0 50];


    h.YMinorTick = 'on';
    h.XMinorTick = 'on';
    h.TickLength = [.04 .04] ;
    h.TickDir    = 'in';
    h.Box               = 'on';
    h.Layer='top'

   
    ypatch_=h.YLim+[1 -1];
    ypatch=[ypatch_ fliplr(ypatch_)];
    hf = patch([h.XLim(1) h.XLim(1) minarray(1) minarray(1)],ypatch ,COLORSpatch(1,:));
    hf.LineStyle='none';
    hf = patch([minarray(1) minarray(1) minarray(2) minarray(2)] ,ypatch, COLORSpatch(2,:));
    hf.LineStyle='none';
    hf = patch([minarray(2) minarray(2) h.XLim(2) h.XLim(2)],ypatch, COLORSpatch(3,:));
    hf.LineStyle='none';

        hl=line([-600,1000],[-600,1000]);
        hl.LineStyle='-';
        hl.Color= [0 0 0];
        hl=errorbar(actualv,vMean(1,:,q),vSD(1,:,q));
        hl.MarkerSize=6;
        hl.Marker='d';
        hl.LineStyle='none';
        hl.Color= col(6,:)/3;
        hl.MarkerFaceColor=col(6,:);

end
%% standard deviations of velocities in ROI at each q
squeeze(vSD(1,1,:))
mean(squeeze(vSD(1,1:5,:)),1)
   %% Save figure.
% fig.Renderer='Painters';
%print(fig,'williamson_figure1_Eq_vmin_infSNR.png','-dpng')
print(fig,'williamson_figure_velocityAccuracy.eps','-depsc')
%print(fig,'williamson_figure1_Eq_vmin_infSNR.tif','-dtiff','-r300')
