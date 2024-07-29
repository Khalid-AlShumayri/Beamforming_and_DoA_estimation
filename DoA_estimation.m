close all
clear all 
clc 

% Comparison between three grid-search approaches
% 1- Conventional Beamformer
% 2- Capon's Beamformer
% 3- spectrum-MUSIC algorithm

M = 5;      % sensors
K = 3;      % number of signals (sources)
N = 10;     % number of observations
d = 0.5;    % Distance between elements in wavelengths
Pn = .09;    % Noise power
sig_pr = 0.9.*ones(1,K);    % signals' power
DoA = -88 + 178*rand(1,K);
Resolution = 2*180/pi/M;
iter = 1;
j = 1;

% Compute steering Matrix and vectors
angles=(-90:.5:90);       % for grid search

% far-field assumption 
a1 = exp(-1i*2*pi*d*(0:M-1)'*sin([angles(:).']*pi/180));  
A  = generate_steering_matrix(M,d,DoA);
S = diag(sqrt(sig_pr./2))*(randn(K,N)+1j*randn(K,N));

Noise = sqrt(Pn/2)*( randn(M,N) + 1j*randn(M,N) );
X = A*S + Noise;

% generate coherent sources 
correlation_matrix = toeplitz([1 .2 .1]);
L = chol(correlation_matrix,"lower");
S_corr = L*S;
X_corr = A*S_corr + Noise;
R_corr = X_corr*X_corr'./N;

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
Qn = Q(:,K+1:M);    % Get the noise eigenvectors

%% Grid search (or we can solve for the roots)
for alpha=1:length(angles)
    % Compute MUSIC spectrum (spatial-spectrum)
    search_music(alpha) = abs( 1/(a1(:,alpha)'*Qn*Qn'*a1(:,alpha)) ); 
    search_bf(alpha) = abs(a1(:,alpha)'*R*a1(:,alpha)/(a1(:,alpha)'*a1(:,alpha)));
    search_cap(alpha) = abs(20/( a1(:,alpha)'*inv(R)*a1(:,alpha) ));
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
% xline(DoA,'-.',LineWidth=.9)

% figure
% plot(angles,spectrum_bf,Color=[0.8, 0.3, 0.1]);hold on; grid on;
% plot(angles,spectrum_cap,Color=[0.5, 0.2, 0.6]);
% xline(DoA,'-.',LineWidth=.9)
% 
% legend('BF','Capon BF','Actual DoAs')
% title(['Spatial Spectrum : SNR = ',num2str(10*log10(max(sig_pr)/Pn)),'dB : Resolution = ',num2str(360/(M)),' degrees'])
% xlabel('Angle in degrees')

%% Peak detection and error

% Step 1: Initial peak detection
[peaks, locs] = findpeaks(spectrum_music);

% Step 2: Determine a dynamic threshold
threshold = mean(spectrum_music);

% Step 3: Identify significant peaks
AoA_peaks = peaks(peaks > threshold);
AoA_music = sort(angles( locs(peaks > threshold) ))
DoA = sort(DoA)
RMSE = sqrt(sum((DoA-AoA_music).^2)/M)
%% Analytical estimation (Validation?)

A_est = X*S'*inv(S*S');      % estimate of the steering matrix


P = S*S'./N;
[Us,Ls] = eig(A*P*A');
[Ls ,Is]  = sort(diag(Ls),1,'descend');   %Find K largest eigenvalues
Us = Us (:,Is);
Us = Us(:,1:K);
Ls = Ls(1:K);
Ls = diag(Ls);


Rn = Noise*Noise'./N;
[Un,Ln] = eig(Rn);
[Ln ,In]  = sort(diag(Ln),1,'descend');   
Un = Un (:,In);
Un = Un(:,K+1:end);
Ln = Ln(K+1:end);
Ln = diag(Ln);

disp('R_sum = A*P*A'' + Cn')
R_hat = A*P*A' + Rn;
[Q_sum,D_sum] = eig(R_hat);
[D_sum,I_sum] = sort(diag(D_sum),1,'descend');
Q_sum = Q_sum(:,I_sum);
Q_sum_s = Q_sum(:,1:K);
Q_sum_n = Q_sum(:,K+1:end);


disp(['norm(R) = ',num2str(norm(R,"fro"))])
disp(['norm(R - R_hat) = ',num2str(norm(R-R_hat,"fro"))])

% AoA = asin( electric_angles./(d*2*pi) )*180/pi
%% Root MUSIC (fft approach)

spectrum_root = fftshift( sum( abs(fft(Qn,3*M)).^2 ,2) );
% trial = fftshift( abs(fft(sum(Qn,2).^2,3*M)) );      % we have to take the fft of each column indiv. 
[mini,I] = mink(spectrum_root,K);
electric_angles = (I-1).*2*pi/M;
% figure
% plot(1./spectrum_root); hold on; grid on;
% plot(1./trial); 
