close all
clear all 
clc 

% Comparison between three grid-search approaches
% 1- Conventional Beamformer
% 2- Capon's Beamformer
% 3- spectrum-MUSIC algorithm

M = 8;      % sensors
K = 2;      % number of signals (sources)
N = 5;     % number of observations
d = 0.5;    % Distance between elements in wavelengths
Pn = .09;    % Noise power
sig_pr = [.9 .9];    % signals' power
DoA = [-20 0];
Resolution = 2*180/pi/M;
iter = 1;
j = 1;

% Compute steering Matrix and vectors
angles=(-90:1:90);       % for grid search

% far-field assumption 
a1 = exp(-1i*2*pi*d*(0:M-1)'*sin([angles(:).']*pi/180));  
A  = generate_steering_matrix(M,d,DoA);
S = diag(sqrt(sig_pr))*(randn(K,N)+1j*randn(K,N));
Noise = sqrt(Pn/2)*( randn(M,N) + 1j*randn(M,N) );
X = A*S + Noise;

%% MUSIC Algorithm

% comparing 3 different techniques
spectrum_music = zeros(iter,length(angles));
spectrum_bf = zeros(iter,length(angles));
spectrum_cap = zeros(iter,length(angles));


R = X*X'./N ;
[Q ,D]  = eig(R);
[D ,I]  = sort(diag(D),1,'descend');   %Find K largest eigenvalues
D = diag(D);
Q = Q (:,I);       % Sort the eigenvectors to put signal eigenvectors first 

Qs = Q (:,1:K);       % Get the signal eigenvectors
Qn = Q(:,K+4:M);    % Get the noise eigenvectors

%% Grid search (or we can solve for the roots)
for k=1:length(angles)
    % Compute MUSIC spectrum (spatial-spectrum)
    search_music(k) = abs( 1/(a1(:,k)'*Qn*Qn'*a1(:,k)) ); 
    search_bf(k) = abs(a1(:,k)'*R*a1(:,k)/(a1(:,k)'*a1(:,k)));
    search_cap(k) = abs(20/( a1(:,k)'*inv(R)*a1(:,k) ));
%     search_cap_cn = 
end
spectrum_music(j,:) = search_music;
spectrum_bf(j,:)    = search_bf;
spectrum_cap(j,:) = search_cap;


plot(angles,spectrum_music,Color=[0 .2 .9]); hold on; grid on;


xline(DoA,'-.',LineWidth=.9)
legend('MUSIC','Actual DoAs')
title(['Spatial Spectrum : SNR = ',num2str(10*log10(max(sig_pr)/Pn)),'dB : Resolution = ',num2str(Resolution),' degrees'])
xlabel('Angle [degrees]')

% figure 
% plot(angles,spectrum_bf,Color=[0.8, 0.3, 0.1]);
% plot(angles,spectrum_cap,Color=[0.5, 0.2, 0.6]);

% figure
% plot(angles,spectrum_bf,Color=[0.8, 0.3, 0.1]);hold on; grid on;
% plot(angles,spectrum_cap,Color=[0.5, 0.2, 0.6]);
% xline(DoA,'-.',LineWidth=.9)
% 
% legend('BF','Capon BF','Actual DoAs')
% title(['Spatial Spectrum : SNR = ',num2str(10*log10(max(sig_pr)/Pn)),'dB : Resolution = ',num2str(360/(M)),' degrees'])
% xlabel('Angle in degrees')

[~,index] = maxk(spectrum_music,K);
estimated_angles_music = angles(index)
%% Analytical estimation

A_est = X*S'*inv(S*S');      % estimate of the steering matrix
P = S*S'./N;
[Es,Ls] = eig(A*P*A');
[Es_hat,Ls_hat] = nnzmat(Es,Ls);
Rn = Noise*Noise'./N;
[En,Ln] = eig(Rn);
Ln = diag(Ln);
Ln_hat = Ln(K+1:end);
Ln_hat = diag(Ln_hat);
En_hat = En(:,K+1:end);
disp('R_hat = A*P*A'' + Cn')
R_hat = A*P*A' + Rn;
[Q_hat,D_hat] = eig(R_hat);
Ln = diag(Ln);

D_hat = sort(diag(D_hat),1,'descend');
% Ls = sort(diag(Ls),1,'descend')
% Ln = sort(diag(Ln),1,'descend')
% D = diag(D)
disp(['norm(R) = ',num2str(norm(R,"fro"))])
disp(['norm(R - R_hat) = ',num2str(norm(R-R_hat,"fro"))])

%% Root MUSIC

Ftr1 = fftshift( sum( abs(fft(Qn,3*M)).^2 ,2) );
Ftr2 = fftshift( abs(fft(sum(Qn,2).^2,3*M)) );      % we have to take the fft of each column indiv. 
[mini,I] = mink(Ftr1,K);
electric_angles = (I-1).*2*pi/M;
figure
plot(1./Ftr1); hold on; grid on;
plot(1./Ftr2); 
legend('original','trial')
% AoA = asin( electric_angles./(d*2*pi) )*180/pi