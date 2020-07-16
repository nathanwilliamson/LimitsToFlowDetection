%% plot velocimetry acuracy
close all
clear all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%% Subtract no flow  y/n = 1/0?
load('velocityData.mat')
[B,I]=sort(velocity.actualVdot);
%ind=1:13;
%ind=[2,3,6,4,9];

D=1.2E-9;
Delta=25E-3;
min=sqrt(2*D/Delta)*1E9*pi*velocity.Rtube^2*velocity.porosity*60;
SNR=70;
minarray=[1/(sqrt(2)*SNR),pi/2]*min;
qrad=velocity.params{1}.q_rad;
qstar=qrad*sqrt(2*D*Delta);
%%%%q^* for beadpack (D_0=1.2E-9)
% % %     0.4869
% % %     0.7794
% % %     1.1690
% % %     1.4615
% % %     1.7539
%%%%q^* for bulk region (D_0=1.89E-9)
% % %     0.6111
% % %     0.9781
% % %     1.4671
% % %     1.8341
% % %     2.2011
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
        VdotMeanbeads(j,k,:)=velocity.beads.VdotMean{I(k)}(:);
        VdotSDbeads(j,k,:)=velocity.beads.VdotSD{I(k)}(:);
        VdotMeanbulk(j,k,:)=velocity.bulk.VdotMean{I(k)}(:);
        VdotSDbulk(j,k,:)=velocity.bulk.VdotSD{I(k)}(:);
        actualVdot(k)=velocity.actualVdot(I(k));

        vErrorbeads(j,k,:)=velocity.beads.VdotMean{I(k)}(:)/velocity.actualVdot(I(k));
        vError(j,k,:)=velocity.bulk.VdotMean{I(k)}(:)/velocity.actualVdot(I(k));

        porosity(j,k,:)=VdotMeanbulk(j,k,:)./VdotMeanbeads(j,k,:);
        
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
    h.YLim=[0 310];
    h.XLim=[0 310];
    h.YLabel.String='measured $\mathrm{\dot{V}}\ \mathrm{[\mu l /min]}$';
    h.XLabel.String='';
    h.FontName          = FontName;
    h.FontSize          = FontSize;
    h.FontWeight        = FontWeight;
     h.Units      = 'centimeters';

     ypatch_=h.YLim+[5 -5];
ypatch=[ypatch_ fliplr(ypatch_)];

hf = patch([h.XLim(1)+.7 h.XLim(1)+.7 minarray(1) minarray(1)],ypatch ,COLORSpatch(1,:));
hf.LineStyle='none';
hf = patch([minarray(1) minarray(1) minarray(2) minarray(2)] ,ypatch, COLORSpatch(2,:));
hf.LineStyle='none';
hf = patch([minarray(2) minarray(2) h.XLim(2)-2 h.XLim(2)-2],ypatch, COLORSpatch(3,:));
hf.LineStyle='none';
     
     
    h.Position(4)          = 3;
    pos_(q,:)= h.Position;
  %  h.Position(3)          = 3.5;

    h.YTick             = [0 100 200 300];
   h.XTick             = [0 100 200 300];
   h.XTickLabel= {'' '' '' ''};
% h.YLabel.Units      = 'centimeters';
% h.YLabel.Position(1) = -0.8;

h.YMinorTick = 'on';
h.XMinorTick = 'on';
h.TickLength = [.02 .02] ;
h.TickDir    = 'out';
h.Box               = 'on';


    hl=line([-600,1000],[-600,1000]);
    hl.LineStyle='-';
    hl.Color= [0 0 0];
    hl=errorbar(actualVdot,VdotMeanbeads(1,:,q),VdotSDbeads(1,:,q));
    hl.MarkerSize=6;
        hl.Marker='d';
        hl.LineStyle='none';
        hl.Color= col(6,:)/3;
        hl.MarkerFaceColor=col(6,:);
    
    hl=errorbar(actualVdot,VdotMeanbulk(1,:,q),VdotSDbulk(1,:,q));
    hl.Marker='o';
    hl.MarkerSize=4;
        hl.LineStyle='none';
        hl.Color= col(7,:)/2;
        hl.MarkerFaceColor=col(7,:);
            
    
    %annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', 'q')
    text(235,120,qtext{q},'FontSize',9)
end
    h.XLabel.String='actual $\mathrm{\dot{V}}\ \mathrm{[\mu l /min]}$';
       h.XTickLabel= {'0' '100' '200' '300'};
for q=1:5
        
    h=axes();
    
    hold on
    h.YLim=[-10 35];
    h.XLim=[-5 35];
    h.YLabel.String='';
    h.XLabel.String='';
    h.FontName          = FontName;
    h.FontSize          = FontSize;
    h.FontWeight        = FontWeight;
    h.Units      = 'centimeters';

    h.Position         = pos_(q,:)+[ .5 1.2 -3 -1];

    h.YTick             = [0 20];
    h.XTick             = [0 20];
% h.YLabel.Units      = 'centimeters';
% h.YLabel.Position(1) = -0.8;

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

%     hl=line([1 1]*minarray(1),[-1000 1000]);
%     hl.LineStyle=':';
%     hl.LineWidth=1;
%     hl.MarkerFaceColor=col(4,:);
%     
%     hl=line([1 1]*minarray(2),[-1000 1000]);
%     hl.LineStyle=':';
%     hl.LineWidth=1;
%     hl.MarkerFaceColor=col(5,:);

    hl=line([-600,1000],[-600,1000]);
    hl.LineStyle='-';
    hl.Color= [0 0 0];
    hl=errorbar(actualVdot,VdotMeanbeads(1,:,q),VdotSDbeads(1,:,q));
    hl.MarkerSize=6;
        hl.Marker='d';
        hl.LineStyle='none';
        hl.Color= col(6,:)/3;
        hl.MarkerFaceColor=col(6,:);
    
    hl=errorbar(actualVdot,VdotMeanbulk(1,:,q),VdotSDbulk(1,:,q));
    hl.Marker='o';
    hl.MarkerSize=4;
        hl.LineStyle='none';
        hl.Color= col(7,:)/2;
        hl.MarkerFaceColor=col(7,:);
            
end
   % h.XLabel.String='actual $v\ \mathrm{[\mu m /s]}$';
   %% Save figure.
% fig.Renderer='Painters';

%print(fig,'williamson_figure1_Eq_vmin_infSNR.png','-dpng')
print(fig,'williamson_figure_flowRateAccuracy.eps','-depsc')
%print(fig,'williamson_figure1_Eq_vmin_infSNR.tif','-dtiff','-r300')
