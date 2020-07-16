function params=parameters(delta,Delta,grad)

Gradient=1.49530722*grad/100;  %grad/100 is a fraction from zero to 1 which specifies the gradient amplitude on the spectrometer
Gamma = 2.6752219E8; %gyromagnetic ratio [rad S^-1 T^-1]
q=Gradient*delta*Gamma/ (2*pi); %q needs to be in units of [1/m], rather than [rad/m] for the the fast fourier transform. a cycle in the fft is 0 to 1.
y=[-flipud(q(2:end));q(1:end-1)];
N = length(y);                  % N=2*L-1; number of points
qinv = ((1:N)' - N/2)/(y(end)-y(1));
dqinv = qinv(2) - qinv(1);
b=Gamma^2.*Gradient.^2.*delta^2.* (Delta-delta/3); %b is the gamma^2 g^2 delta^2 (Delta - Delta/3)

params.delta=delta;
params.Delta=Delta;
params.Gradient=Gradient;
params.Gamma=Gamma;
params.q=q;
params.q_rad=q*2*pi;  % needs to be in radians/m for the velocity from phase shift calculation
params.N=N;
params.qinv=qinv;
params.dqinv=dqinv;
params.b=b;
params.venc=pi./(q*2*pi*Delta);
end