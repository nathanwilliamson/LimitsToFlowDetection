%% plot velocity image
%%% need to increase length to width ratio by 1.5 to be true image.
%%% 1st column images are bigger than second by 3.5mm (25.5 vs 21.8)
clear all
%close all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%ind=[1,3,4,9];

%ind=[9,5,13,12];
%ind=[8,2,3,4];


%%%slow flows
ind=[1,2,3,4];

%%%fast flows
%ind=[9,5,13,12];

load('velocityData.mat')
magmin=0;
magmax=0;
vmin=0;
vmax=0;
q=4;
%vrange=[0,700];
vrange=[-20,70];
col=get(groot,'DefaultAxesColorOrder');


fig                 = figure();
fig.Units           = 'centimeters';
fig.PaperUnits      = 'centimeters';
fig.Position        = [0 0 16 20];
fig.PaperPosition   = fig.Position;

FontName            = 'helvetica';
FontSize            = 7;
FontWeight          = 'normal';

pixelSize=15/64; %(FOV/#pts)
sliceSize=22/32;
x=[0:(length(velocity.ImageMag{1}(:,1,1))-1)]*pixelSize;
y=[0:(length(velocity.ImageMag{1}(1,:,1))-1)]*sliceSize;
%vtext={'$\mathrm{\dot{V}}=60\ \mathrm{\mu l /min}$', '$\mathrm{\dot{V}}=100\ \mathrm{\mu l /min}$', '$\mathrm{\dot{V}}=200\ \mathrm{\mu l /min}$', '$\mathrm{\dot{V}}=300\ \mathrm{\mu l /min}$'} ;
vtext={'$\mathrm{\dot{V}}=0\ \mathrm{\mu l /min}$', '$\mathrm{\dot{V}}=3\ \mathrm{\mu l /min}$', '$\mathrm{\dot{V}}=10\ \mathrm{\mu l /min}$', '$\mathrm{\dot{V}}=30\ \mathrm{\mu l /min}$'} ;

%imgtext1={'velocity images' 'phase-corrected\n velocity images'};
%imgtext2={'velocity profile' 'phase-corrected\n velocity profile'};

% for k=1:length(ind)
%     ax1=subplot(length(ind),5,(k-1)*5+1);
%     visImage=velocity.ImageMag{ind(k)}(:,:,q);
%     imagesc(flipud(visImage'),'XData', [0 length(velocity.ImageMag{1}(:,1,1))*pixelSize],'YData', [0 length(velocity.ImageMag{1}(:,1,1))*sliceSize]) ;
%     caxis(ax1,[0 2E6])
%     colormap(ax1,gray)
%     
% 
% end
    
j=1
    subNoFlow=j-1;%%
    if subNoFlow==1
        load('velocityData_subNoFlow.mat')
        velocity=velocity_subNoFlow;
    elseif subNoFlow==0
        load('velocityData.mat')
    end
    
    for k=1:length(ind)
        h=subplot(length(ind),4,(k-1)*4+1);
        hold on
        h.FontName          = FontName;
        h.FontSize          = FontSize;
        h.FontWeight        = FontWeight;
        h.Units      = 'centimeters';
        h.XLabel.String='';
        h.YLabel.String='';
        h.YLim=[0-.5 length(velocity.ImageMag{1}(1,:,1))+.5];
        h.XLim=[0-.5 length(velocity.ImageMag{1}(:,1,1))+.5];

         h.XTick=[];
         h.YTick=[];

         h.Position(3)          = h.Position(3)+0.4; 
         h.Position(4)          = h.Position(4)+0.5;
%        
        visImage=velocity.Image{ind(k)}(:,:,q)*1E6;
        imAlpha=ones(size(visImage));
        imAlpha(isnan(visImage))=0;
        imagesc(visImage,'XData', [0+.5 length(visImage(1,:))-.5],'YData', [0+.5 length(visImage(:,1))-.5],'AlphaData',imAlpha) ;
        h.XLim=[0 length(visImage(1,:))];
        h.YLim=[0 length(visImage(:,1))]; 
        set(gca,'color',0*[1 1 1]);
        %imagesc(visImage') ;
        colormap(h,jet);
        caxis(h,vrange);
         
         
        hc=colorbar;
        set( hc, 'YDir', 'reverse' );
%        if k==1
%            text(2,5,imgtext1{j},'FontSize',FontSize)
%        end
       rect=rectangle('Position',[ 0 1 24 3]);
       rect.LineStyle='--';
       rect.LineWidth=2;
       rect.EdgeColor=[1 1 1];
       
       rect=rectangle('Position',[ 0 1 24 3]);
       rect.LineStyle='--';
       rect.LineWidth=1;
       rect.EdgeColor=col(6,:);
       
       rect=rectangle('Position',[0 7  24 3]);
       rect.LineStyle='--';
       rect.LineWidth=2;
       rect.EdgeColor=col(7,:);

            ha=annotation('rectangle');
            ha.Units='centimeters';
            ha.Position(1)=h.Position(1);
            ha.Position(2)=h.Position(2)+3.3+0.6;
            ha.Position(3)=12.45;
            ha.Position(4)=0.01;
            ha.Color=[.8 0 0];
            ha.FaceColor=[.8 0 .0];

            ha=annotation('textbox');
            ha.Interpreter='latex';
            ha.String=vtext{k};
            ha.FontSize=9;
            ha.Units='centimeters';
            %ha.Position(3)=0.3
            ha.Position(1)=h.Position(1)+12.45-2.5;
            ha.Position(2)=h.Position(2)+3.3;
            ha.Position(3)=2.5;
            ha.Position(4)=.6;
            ha.Color=[0 0 0];
            ha.BackgroundColor=[.95 .95 .95];
            ha.EdgeColor= [.8 0 0];
            
        %if min(min(visImage))<vmin
            vmin(j,k)=min(min(visImage));
       % end
       % if max(max(visImage))>vmax
            vmax(j,k)=max(max(visImage));
       % end
    end


%% plot profiles
Rtube=2.5; % mm
rshift=-2.8;
rsim=linspace(-Rtube,Rtube,40);
r=([0:velocity.beads.dims(2)-1]*pixelSize+rshift);



    subNoFlow=j-1;%%
    if subNoFlow==1
        load('velocityData_subNoFlow.mat')
        velocity=velocity_subNoFlow;
    elseif subNoFlow==0
        load('velocityData.mat')
    end
    
    for k=1:length(ind)
        h=subplot(length(ind),4,(k-1)*4+2);
        hold on
        h.FontName          = FontName;
        h.FontSize          = FontSize;
        h.FontWeight        = FontWeight;
        h.Units      = 'centimeters';
        h.XLabel.String='';
        h.Position(4)          = h.Position(4)+0.5;
        
        %if j==1
        h.YLabel.String='$v\ \mathrm{[\mu m /s]}$ ';
        %end
        %      h.XTick=[];
       %h.YTick=[-1:6]*20;
        h.YLim=vrange;
        set( h, 'YDir', 'reverse' );
        hold on
        hl=errorbar(r,velocity.bulk.vprofileMean{ind(k)}(:,q)*1E6,velocity.bulk.vprofileSD{ind(k)}(:,q)*1E6);
        hl.Marker='o';
        hl.LineStyle='none';
        hl.Color= [0 0 0];
        hl.MarkerFaceColor=col(7,:);
            

        
        hl=errorbar(r,velocity.beads.vprofileMean{ind(k)}(:,q)*1E6,velocity.beads.vprofileSD{ind(k)}(:,q)*1E6);
        hl.Marker='d';
        hl.LineStyle='none';
        hl.Color= [0 0 0];
        hl.MarkerFaceColor=col(6,:);
        v_r = 2*velocity.bulk.actualv{ind(k)}.*(1-((rsim)./Rtube).^2)*1E6 ;
        hl=line(rsim,v_r);
        hl.LineWidth=1;
        hl.LineStyle='-';
        hl.Color=col(7,:);
        v_r = velocity.beads.actualv{ind(k)}*1E6.*ones(40,1);
        hl=line(rsim,v_r);
        hl.LineWidth=1;
        hl.LineStyle='-';
        hl.Color=col(6,:);
        
    end
            h.XLabel.String='radius [mm]';


%% plot axial images

  for k=1:length(ind)
        h=subplot(length(ind),4,(k-1)*4+3);
        hold on
        h.FontName          = FontName;
        h.FontSize          = FontSize;
        h.FontWeight        = FontWeight;
        h.Units      = 'centimeters';
        h.XLabel.String='';
        h.YLabel.String='';
        h.YLim=[0-.5 length(velocity.bulk.ImageSlice{1}(1,:,1))+.5];
        h.XLim=[0-.5 length(velocity.bulk.ImageSlice{1}(:,1,1))+.5];

         h.XTick=[];
         h.YTick=[];


             h.Position(2)          = h.Position(2)+.3;

         h.Position(4)          = h.Position(3);
        
        visImage=velocity.bulk.ImageSlice{ind(k)}(:,:,q)*1E6;
        imAlpha=ones(size(visImage));
        imAlpha(isnan(visImage))=0;
        imagesc(visImage,'XData', [0+.5 length(visImage(1,:))-.5],'YData', [0+.5 length(visImage(:,1))-.5],'AlphaData',imAlpha) ;
        h.XLim=[0 length(visImage(1,:))];
        h.YLim=[0 length(visImage(:,1))];
        set(gca,'color',0*[1 1 1]);
        colormap(h,jet);
        caxis(h,vrange);
         
       rect=rectangle('Position',[0 length(visImage(:,1))/2-.5  length(visImage(:,1)) 1]);
       rect.LineStyle='--';
       rect.LineWidth=2;
       rect.EdgeColor=col(7,:);
  end

  for k=1:length(ind)
        h=subplot(length(ind),4,(k-1)*4+4);
        hold on
        h.FontName          = FontName;
        h.FontSize          = FontSize;
        h.FontWeight        = FontWeight;
        h.Units      = 'centimeters';
        h.XLabel.String='';
        h.YLabel.String='';
        h.YLim=[0-.5 length(velocity.beads.ImageSlice{1}(1,:,1))+.5];
        h.XLim=[0-.5 length(velocity.beads.ImageSlice{1}(:,1,1))+.5];

         h.XTick=[];
         h.YTick=[];


% 
             h.Position(2)          = h.Position(2)+.3;
% 
         h.Position(4)          = h.Position(3);
%        
        visImage=velocity.beads.ImageSlice{ind(k)}(:,:,q)*1E6;
        imAlpha=ones(size(visImage));
        imAlpha(isnan(visImage))=0;
        imagesc(visImage,'XData', [0+.5 length(visImage(1,:))-.5],'YData', [0+.5 length(visImage(:,1))-.5],'AlphaData',imAlpha) ;
        h.XLim=[0 length(visImage(1,:))];
        h.YLim=[0 length(visImage(:,1))];
        set(gca,'color',0*[1 1 1]);
        colormap(h,jet);
        caxis(h,vrange);
         
         
       rect=rectangle('Position',[0 length(visImage(:,1))/2-.5  length(visImage(:,1)) 1]);
       rect.LineStyle='--';
       rect.LineWidth=2;
       rect.EdgeColor=col(6,:);
  end  
  
%% Save figure.
% fig.Renderer='Painters';
%set(gcf, 'InvertHardcopy', 'off')

print(fig,'williamson_figure_velocityImgNew_q4.png','-dpng')
print(fig,'williamson_figure_velocityImgNew_q4.eps','-depsc')

%print(fig,'williamson_figure_velocityImgNew_q2_fastFlows.png','-dpng')
%print(fig,'williamson_figure_velocityImgNew_q2_fastFlows.eps','-depsc')
%print(fig,'williamson_figure1_Eq_vmin_infSNR.tif','-dtiff','-r300')