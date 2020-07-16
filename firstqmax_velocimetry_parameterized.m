%% location of first imag(E(q)) maximum 
% Nathan Williamson 3/20/2018
%https://www.mathworks.com/help/symbolic/solve-an-algebraic-equation.html

close all
clear all
clc
%% define parameters
Delta = 1; %s
D = 2E-9;                    % Diffusion coefficient [m^2/s] 

v=linspace(0.6,200,100)'*1E-6; %m/s
firstq=zeros(length(v),1);
Ei_firstq=zeros(length(v),1);
ub=2E4;

q_fast = pi./(2*v.*Delta);  %prediction based on velocity only
q_slow = 1./sqrt(2*D*Delta)*ones(length(v),1);% prediction in the limit that v goes to 0

Ei_firstq_slow=sin(q_slow.*v.*Delta).*exp(-q_slow.^2.*D.*Delta)*sqrt(2); %dividing by sqrt(2) to take into account the 2-point method
Ei_firstq_fast=sin(q_fast.*v.*Delta).*exp(-q_fast.^2.*D.*Delta)*sqrt(2);

for i=1:length(v)
    syms f(q)
    f(q) = v(i)/(2*q*D) - tan(q*v(i)*Delta);

%     analyticalq_lb=1/sqrt(2*D*Delta)
%     analyticalq_ub=pi/(2*v*Delta)

    sol=vpasolve(f,q,[0,ub]);
    Ei_firstq(i)=sin(sol*v(i)*Delta)*exp(-sol^2*D*Delta)*sqrt(2);   %dividing by sqrt(2) to take into account the 2-point method
    firstq(i)=sol;
    ub=sol;
    clear sol f(q)
end


%% plot
 
FontName='helvetica';
FontSize=7;
FontWeight='normal';
LineWidth = 1;

fig                 = figure();
fig.Units           = 'centimeters';
fig.PaperUnits      = 'centimeters';
fig.Position        = [0 0 8 8];
fig.PaperPosition   = fig.Position;

%subplot(2,1,1) 



h                   = axes();
h.Units             = 'centimeters';
h.FontName          = FontName;
h.FontSize          = FontSize;
h.FontWeight        = FontWeight;
h.Position          = [1.25 5 6.5 2.5];
h.Box               = 'on';
%h.YLabel.String = '$q \sqrt{2 D\Delta}$ [rad]';
%h.XLabel.String = '$v /\sqrt{2 D/\Delta}$';
h.YLabel.String = '$q \times l_D$  [rad]';
h.XLabel.String = '$v/v_{SF} $';
h.YLabel.Interpreter = 'latex';
h.XLabel.Interpreter = 'latex';

h.YLim = [0 1.75];
h.XLim = [0 2];
% h.YTick =[];
% h.XTick =[];
%h.XLim = [-2 4]*1E-4;

hold on
ha=annotation('textbox');
        ha.Interpreter='latex';
        ha.String='a';
        ha.FontSize=9;
        ha.Units='centimeters';
        %ha.Position(3)=0.3
        ha.Position(1)=h.Position(1)-0.85;
        ha.Position(2)=h.Position(2)+2.5;
        ha.Position(3)=0.5;
        ha.Position(4)=.5;
        ha.Color=[0 0 0];
        ha.BackgroundColor=[.95 .95 .95];
        ha.EdgeColor= [.8 0 0];
        
COLORSpatch=[255 255 234; 234 255 249; 244 240 255]/255;
COLORS = 1/255 * [  55  80  162 ; ...
    93  187 70  ; ...
    241 156 31  ; ...
    237 28  36  ; ...
    129 41  134 ];

ypatch_=h.YLim+[5 -5];
ypatch=[ypatch_ fliplr(ypatch_)];


hf = patch([.002 .002 1 1] ,ypatch, COLORSpatch(2,:));
hf.LineStyle='none';
hf = patch([1 1 h.XLim(2)-0.002 h.XLim(2)-0.002],ypatch, COLORSpatch(3,:));
hf.LineStyle='none';        

%x = [ .35, .25];
x = [ .3, .22];
y = [.67, .67];
%a = annotation('arrow',x,y);
a = annotation('textarrow',x,y,'String','slow-flow regime');
a.FontSize= 7;
a.Interpreter='latex';
%x = [ .5, .6];
x = [ .81, .89];
y = [.67, .67];
%a = annotation('arrow',x,y);
a = annotation('textarrow',x,y,'String','fast-flow regime ');
a.FontSize= 7;
a.Interpreter='latex';


hl1 = line(v/(pi/2*sqrt(2*D/Delta)),firstq*sqrt(2*D*Delta));
hl1.LineStyle = '-';
hl1.LineWidth=LineWidth;
% hl.Marker           = '.';
% hl.MarkerSize       = 8;
hl1.Color = 'k';
hl1.DisplayName='maximum';


hl2 = line(v/(pi/2*sqrt(2*D/Delta)),q_fast*sqrt(2*D*Delta));
hl2.LineStyle = '-';
hl2.LineWidth=LineWidth;
% hl.Marker           = '.';
% hl.MarkerSize       = 8;
hl2.DisplayName='$q_{fast}$';
hl2.Color = COLORS(1,:);

hl3 = line(v/(pi/2*sqrt(2*D/Delta)),q_slow*sqrt(2*D*Delta));
hl3.LineStyle = '-';
hl3.LineWidth=LineWidth;
% hl.Marker           = '.';
% hl.MarkerSize       = 8;
hl3.DisplayName='$q_{slow}$';
hl3.Color = COLORS(2,:);

h.YLabel.Position(1) = -.15;


  legend([hl1 hl2 hl3],'Box','off')
 
%% plot signal intensity of imaginary channel
h                   = axes();
h.Units             = 'centimeters';

h.FontName          = FontName;
h.FontSize          = FontSize;
h.FontWeight        = FontWeight;
h.Position          = [1.25 1.25 6.5 2.5];
h.Box               = 'on';
h.YLim = [0 1.25];
h.XLim = [0 2];
h.YLabel.String = '$1/\mathrm{SNR}$';
%h.XLabel.String = '$v /\sqrt{2 D/\Delta}$';
h.XLabel.String = '$v/v_{SF}$';
h.YLabel.Interpreter = 'latex';
h.XLabel.Interpreter = 'latex';
h.YLabel.Position(1) = -.15;
hold on

ha=annotation('textbox');
        ha.Interpreter='latex';
        ha.String='b';
        ha.FontSize=9;
        ha.Units='centimeters';
        %ha.Position(3)=0.3
        ha.Position(1)=h.Position(1)-0.85;
        ha.Position(2)=h.Position(2)+2.5;
        ha.Position(3)=0.5;
        ha.Position(4)=.5;
        ha.Color=[0 0 0];
        ha.BackgroundColor=[.95 .95 .95];
        ha.EdgeColor= [.8 0 0];
        
ypatch_=h.YLim+[5 -5];
ypatch=[ypatch_ fliplr(ypatch_)];


hf = patch([.002 .002 1 1] ,ypatch, COLORSpatch(2,:));
hf.LineStyle='none';
hf = patch([1 1 h.XLim(2)-0.002 h.XLim(2)-0.002],ypatch, COLORSpatch(3,:));
hf.LineStyle='none';
        
hl = line(v/(pi/2*sqrt(2*D/Delta)),Ei_firstq);
hl.LineStyle = '-';
hl.LineWidth=LineWidth;
% hl.Marker           = '.';
% hl.MarkerSize       = 8;
hl.Color = 'k';

hl = line(v/(pi/2*sqrt(2*D/Delta)),Ei_firstq_slow);
hl.LineStyle = '-';
hl.LineWidth=LineWidth;
% hl.Marker           = '.';
% hl.MarkerSize       = 8;
hl.Color = COLORS(2,:);

hl = line(v/(pi/2*sqrt(2*D/Delta)),1./(2.3./v.*sqrt(D/Delta))*sqrt(2)); %taking sqrt(2) into account for 2-point method
hl.LineStyle = '--';
hl.LineWidth=LineWidth;
% hl.Marker           = '.';
% hl.MarkerSize       = 8;
hl.Color = COLORS(2,:);

hl = line(v/(pi/2*sqrt(2*D/Delta)),Ei_firstq_fast);
hl.LineStyle = '-';
hl.LineWidth=LineWidth;
% hl.Marker           = '.';
% hl.MarkerSize       = 8;
hl.Color = COLORS(1,:);

print(fig,'williamson_figure2_firstqmax.eps','-depsc')
%print(fig,'williamson_figure2_firstqmax_bw.eps','-deps')
print(fig,'williamson_figure2_firstqmax.tif','-dtiff','-r300')
