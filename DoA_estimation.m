close all
clear all 
clc 

% Comparison between three grid-search DoA approaches
% 1- Conventional Beamformer
% 2- Capon's Beamformer
% 3- spectrum-MUSIC algorithm

M = 64;      % sensors
K = 3;      % number of signals (sources)
N = 1;     % number of observations
d = 0.5;    % Distance between elements in wavelengths
Pn = .09;    % Noise power
sig_pr = 0.9.*ones(1,K);    % signals' power
DoA = [-30 20 35];
Resolution = 2*180/pi/M;
iter = 1;
j = 1;

% Compute steering Matrix and vectors
angles=(-90:.1:90);       % for grid search

% far-field assumption 
a1 = exp(-1i*2*pi*d*(0:M-1)'*sin([angles(:).']*pi/180));  
A  = generate_steering_matrix(M,d,DoA);

S = diag(sqrt(sig_pr./2))*(randn(K,N)+1j*randn(K,N));

Noise = sqrt(Pn/2)*( randn(M,N) + 1j*randn(M,N) );
X = A*S + Noise;
R = X*X'./N ;
%% MUSIC Algorithm

% comparing the three different techniques
spec_music      = zeros(iter,length(angles));
spec_music_coh = zeros(iter,length(angles));
spec_bf = zeros(iter,length(angles));
spec_cap = zeros(iter,length(angles));

% Uncorrelated signals
[Q ,D]  = eig(R);
[D ,I]  = sort(diag(D),1,'descend');   %Find K largest eigenvalues
D = diag(D);
Q = Q (:,I);       % Sort the eigenvectors to put signal eigenvectors first 
Qs = Q (:,1:K);       % Get the signal eigenvectors
Qn = Q(:,K+1:M);    % Get the noise eigenvectors


%% Grid-search 
for i=1:length(angles)
    % Compute MUSIC spectrum (spatial-spectrum)
    srch_music(i)     = abs( 1/(a1(:,i)'*Qn*Qn'*a1(:,i)) ); 
    srch_bf(i) = abs(a1(:,i)'*R*a1(:,i)/(a1(:,i)'*a1(:,i)));
    srch_cap(i) = abs(1/( a1(:,i)'*inv(R)*a1(:,i) ));
%     search_cap_cn = 
end
spec_music(j,:)  = srch_music;
spec_bf(j,:)    = srch_bf;
spec_cap(j,:) = srch_cap;


plot(angles,spec_music,Color=[0 .2 .9]); hold on; grid on;
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

% % Peak detection and error

% Step 1: Initial peak detection
[peaks, locs] = findpeaks(spec_music);

% Step 2: Determine a dynamic threshold
threshold = mean(spec_music) + 3*std(spec_music);

% Step 3: Identify significant peaks
AoA_peaks = peaks(peaks > threshold);
AoA_music = sort(angles( locs(peaks > threshold) ))
DoA = sort(DoA)
RMSE = sqrt(sum((DoA-AoA_music).^2)/M)
%% Analytical estimation (Validation?)

% A_estm = X*S'*inv(S*S');      % estimate of the steering matrix
% theta_estm = angle(A_estm(2,:));

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

R_hat = A*P*A' + Rn;
[Q_sum,D_sum] = eig(R_hat);
[D_sum,I_sum] = sort(diag(D_sum),1,'descend');
Q_sum = Q_sum(:,I_sum);
Q_sum_s = Q_sum(:,1:K);
Q_sum_n = Q_sum(:,K+1:end);


disp(['norm(R) = ',num2str(norm(R,"fro"))])
disp(['norm(R - R_hat) = ',num2str(norm(R-R_hat,"fro"))])

% AoA = asin( electric_angles./(d*2*pi) )*180/pi
