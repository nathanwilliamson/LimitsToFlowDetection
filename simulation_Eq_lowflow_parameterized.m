%Test fast fourier transform of q-space data for diffusion and flow (from Bloch-Torrey) 
% to find propagator
%simulates data with the possibility of a g offset, like we do in our
%experiments.
% Nathan Williamson 10/18/2017
close all
clear all
clc
%% define parameters
delta = 0.001   ;            % delta
Delta = 1    ;           % Delta
Gamma  = 2.6752219E8;          % Gammamagnetic ratio [rad S^-1 T^-1]
gstart = 0.000;               % gradient start. We assume this value to be 0 in the F.T.
gmax = 4;                % gradient end
gnum = 10000;              % number of gradient pts. An odd number will lead to an even number for the FT, and vice-versa, but this doesn't seem to matter. We assume linear spacing in the FT.    
g = linspace(gstart,gmax,gnum); %gradient vecor
q = g *Gamma*delta/ (2*pi) ;      % q-space. q needs to be in units of [1/m], rather than [rad/m] for the the fast fourier transform. a cycle in the fft is 0 to 1. Note the units are [radians/m] for q-space in callaghan. In the S.T. Relation it needs to be rad/s because by euler's formula we are really talking about rotations in radians. 
b_=Gamma^2.*(g.^2).*delta^2 .* (Delta-delta/3);

%% simulate data based on bloch torrey diffusion and flow
%v= linspace(0,60E-6,4);
D = 2E-9;                    % Diffusion coefficient [m^2/s] 

v=[0.01 0.1 1 10]*(pi/2*sqrt(2*D/Delta));
vAxis=[0 1 2 3];
%v = 200E-6;                   % mean velocity [m/s]
for j=1:length(v);
E(j,:) = exp(1i*(2*pi*q).*v(j).*Delta-(2*pi*q).^2*D*(Delta-delta/3)); % multiply q by 2*pi because In the S.T. relation it needs to be rad/s because by euler's formula we are really talking about rotations in radians.  
%  v = -0.002;                   % mean velocity [m/s]
% D = 2E-9;                    % Diffusion coefficient [m^2/s] 
%  E =E + exp(1i*(2*pi*q)*v*Delta-(2*pi*q).^2*D*(Delta-delta/3)); 
SNR=100;
E(j,:) = E(j,:) + randn(1,gnum)*abs(E(j,1))/SNR+1i*randn(1,gnum)*abs(E(1))/SNR;                                
end


%% For FT 
 for j=1:length(v)
 %E(j,:) = E(j,:).*exp(-1i*angle(E(j,1)));        % phase data so that the first E(q) point is entirely real.     

 % E = E.*exp(-1i*angle(E(1))).*exp(+1i*pi/2);        % phase data so that the first E(q) point is entirely imaginary.     
 y(j,:)=[-fliplr(q(2:end)),q(1:end-1)];
 xdim = length(y);                  % number of points
 E_ft(j,:)=[conj(fliplr(E(j,2:end))),E(j,1:end-1)];   % we do the conjugate of E for negative q because the imaginary signal has opposite sign after reflection.
 end
 
FontName='helvetica';
FontSize=13;
FontWeight='normal';
LineWidth = 1;
 
fig = figure();
fig.Units = 'centimeters';
fig.PaperUnits = 'centimeters';
fig.Position = [0 0 18 17];
fig.PaperPosition = fig.Position;
%subplot(2,1,1) 
h = axes();
h.Units='centimeters';
h.FontName=FontName;
h.FontSize=FontSize;
h.FontWeight=FontWeight;
h.Position = [2 2 15 15];
%h.XLabel.String = '$q \sqrt{2 D\Delta}$ [rad]';
h.XLabel.String = '$q\times l_D$ [rad]';
%h.XLabel.Position= [-1 -1.0036e+04 -1.5334];
h.XAxisLocation = 'origin';
%h.YLabel.String = '$v /\sqrt{2 D/\Delta}$';
h.YLabel.String = '$v/v_{SF}$';
%h.YLabel.Position= [ 8.5 0 -0.4];
h.ZLabel.String = '$E(q)$';
h.YLabel.Interpreter = 'latex';
h.XLabel.Interpreter = 'latex';
h.ZLabel.Interpreter = 'latex';
 h.YTick = vAxis/(pi/2*sqrt(2*D/Delta));
 h.ZLim = [-1,1.05];

h.YTickLabel =[0.01,0.1,1,10];
%h.YTickLabel =[.6,6,60,600];

h.XLim = [-4 4];


hold on
yMat=repmat(vAxis/(pi/2*sqrt(2*D/Delta)),numel(y(1,:)),1)';
plot3(y'*2*pi*sqrt(2*D*Delta),yMat',real(fliplr(E_ft)'),'b','LineWidth',LineWidth)

plot3(y'*2*pi*sqrt(2*D*Delta),yMat',imag(fliplr(E_ft)'),'r','LineWidth',LineWidth)
plot3(y'*2*pi*sqrt(2*D*Delta),yMat',zeros(length(y),length(v)),'k','LineStyle',':','LineWidth',1)

grid;
%xlabel('x'); ylabel('y'); zlabel('z');
view(40,40); %// Adjust viewing angle so you can clearly see data

%% annotation
COLORSpatch=[255 255 234; 234 255 249; 244 240 255]/255;

%x = [ .35, .25];
x=[0.744753086419753 0.572283950617284];
y=    [0.28945643153527 0.172199170124481];
a = annotation('arrow',x,y);
%x = [ .5, .6];
x=[0.784 0.891743827160494];
y=    [0.311 0.387966804979253];
%a = annotation('arrow',x,y);
a = annotation('arrow',x,y);
%annotation(fig,'arrow',[0.556321713870733 0.472352941176471],...
%    [0.158751037344399 0.0954356846473029]);
annotation(fig,'textbox',...
    [0.745 0.335 0.080 0.031],...
    'String',{'fast-flow regime'}, 'FitBoxToText','on','Interpreter','latex','BackgroundColor',COLORSpatch(3,:),'FaceAlpha',1);
annotation(fig,'textbox',...
    [0.61,0.23,0.1,0.046],...
    'String',{'slow-flow regime'}, 'FitBoxToText','on','Interpreter','latex','BackgroundColor',COLORSpatch(2,:),'FaceAlpha',1);
%annotation(fig,'textbox',...
 %   [0.45 0.12 0.19 0.04],...
  %  'String',{'unmeasurable flow'}, 'FitBoxToText','off','Interpreter','latex','BackgroundColor',COLORSpatch(1,:),'FaceAlpha',1);

%% Save figure.
 
 fig.Renderer='Painters';
%print(fig,'williamson_figure1_Eq_vmin_infSNR.png','-dpng')
print(fig,'williamson_figure1_Eq_vmin.eps','-depsc')
%print(fig,'williamson_figure1_Eq_vmin_infSNR.tif','-dtiff','-r300')
 %E_ft=abs(E_ft);
% % X = (1:xdim) - xdim/2;      % x-domain variable
%  xshift = 0;                  % not necessary since q-space is always max at q=0, 
%  y = X - (xshift+1);         % zero in the middle of the x-domain

