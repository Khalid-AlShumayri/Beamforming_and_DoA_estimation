close all
clear all 
clc 

M = 128;      % sensors
K = 3;      % number of signals (sources)
N = 1;     % number of observations
d = 0.5;    % Distance between elements in wavelengths
Pn = .9;    % Noise power
sig_pr = 0.9.*ones(1,K);    % signals' power
DoA = [-30 20 35];
L = floor(M / 2);  % Length of subarrays
snr = min(sig_pr)/Pn;
Resolution = (2/(L*snr))*180/pi;
iter = 1;
j = 1;

% Compute steering Matrix and vectors
angles=(-90:.1:90);       % for grid search

% far-field assumption 
A  = generate_steering_matrix(M,d,DoA);
a_srch = exp(-1i*2*pi*d*(0:L-1)'*sin([angles(:).']*pi/180));

S = diag(sqrt(sig_pr./2))*(randn(K,N)+1j*randn(K,N));

Noise = sqrt(Pn/2)*( randn(M,N) + 1j*randn(M,N) );
X = A*S + Noise;
% R_main = 
%% generate coherent sources (Needs modification)
corr_matrix = [1   .2   0;
               .2  1    0;
               0   0    1];

Lr = chol(corr_matrix,"lower");
S_corr = Lr*S;
X_corr = A*S_corr + Noise;
R_coh = X_corr*X_corr'./N;

%% MUSIC Algorithm

% comparing the three different techniques
spec_music      = zeros(iter,length(angles));
spec_music_coh = zeros(iter,length(angles));

% Compute spatially smoothed covariance matrix
R = zeros(L, L);
for k = 1:(M - L + 1)
    x_sub = X(k:k + L - 1, :);
    R = R + (x_sub * x_sub') ./ (M - L + 1);
end

%% Uncorrelated signals
[Q ,D]  = eig(R);
[D ,I]  = sort(diag(D),1,'descend');   %Find K largest eigenvalues
D = diag(D);
Q = Q (:,I);       % Sort the eigenvectors to put signal eigenvectors first 
Qs = Q (:,1:K);       % Get the signal eigenvectors
Qn = Q(:,K+1:L);    % Get the noise eigenvectors


%%  Correlated signals (coherent sources) (Needs revision)
% [Q_coh ,D_coh]  = eig(R_coh);
% [D_coh ,I_coh]  = sort(diag(D_coh),1,'descend');   %Find K largest eigenvalues
% D_coh = diag(D_coh);
% Q_coh = Q_coh(:,I);       % Sort the eigenvectors to put signal eigenvectors first 
% Qs_coh = Q_coh(:,1:K);       % Get the signal eigenvectors
% Qn_coh = Q_coh(:,K+1:M);    % Get the noise eigenvectors

%% Grid-search 
for i=1:length(angles)
    % Compute MUSIC spectrum (spatial-spectrum)
    srch_music(i)     = abs( 1/(a_srch(:,i)'*Qn*Qn'*a_srch(:,i)) ); 
    srch_music_diag(i)     = abs( 1/(a_srch(:,i)'*Qn_d*Qn_d'*a_srch(:,i)) );
%     srch_music_corr(i)= abs( 1/(A_srch(:,i)'*Qn_coh*Qn_coh'*A_srch(:,i)) );

end
spec_music  = srch_music;
spec_music_diag = srch_music_diag;
% spec_music_coh = srch_music_corr;

% Plot 
plot(angles,spec_music,Color=[0 .2 .9]); hold on; grid on;
plot(angles,spec_music_diag,Color=[.9 .2 0])
% plot(angles,spec_music_coh,"Color",[.8 .2 .1]) 
xline(DoA,'-.',LineWidth=.9)

title(['Spatial Spectrum : SNR = ',num2str(10*log10(min(sig_pr)/Pn)),'dB : Resolution = ',num2str(Resolution),' degrees'])
xlabel('Angle [degrees]')

%% Peak detection and error

% Step 1: Initial peak detection
[peaks, locs] = findpeaks(spec_music);

% Step 2: Determine a dynamic threshold
threshold = mean(spec_music) + 1.5*std(spec_music);
yline(threshold,'--',LineWidth=.9)
legend('MUSIC','MUSIC w DL','Actual DoAs','threshold')

% Step 3: Identify significant peaks
ind = sort(angles( locs(peaks > threshold) ));
AoA_music = zeros(1,K);
AoA_music(1,1:length(ind)) = ind

DoA = sort(DoA)
RMSE = sqrt(sum((DoA-AoA_music).^2)/K)


% AoA = asin( electric_angles./(d*2*pi) )*180/pi

%% Root MUSIC (fft approach)

N_dft = 2*L;
spectrum_root =  sum( abs(fft(Qn,N_dft)).^2 ,2);
% trial = fftshift( abs(fft(sum(Qn,2).^2,3*M)) );      % we have to take the fft of each column indiv. 

% % plot
freq = [-N_dft/2:1:N_dft/2-1]/N_dft;
index = asin(freq/d)*180/pi;

% % Peak detection and error
 
% Step 1: Initial peak detection
[peaks2, locs2] = findpeaks(1./spectrum_root);

% Step 2: Determine a dynamic threshold
threshold2 = mean(1./spectrum_root);

% Step 3: Identify significant peaks
ind2 = sort(index( locs2(peaks2 > threshold2) ));

AoA_fft = zeros(1,K);
AoA_fft(1,1:length(ind2)) = ind2
RMSE_root = sqrt(sum((DoA-AoA_fft).^2)/K)

figure
plot(index,1./spectrum_root); hold on; grid on;
yline(threshold2,'--')
xline(DoA,'--')
title('FFT')
xlabel('Angles [Degrees]')
ylabel('PSD')
%% Analytical estimation (Validation?)

% A_estm = X*S'*inv(S*S');      % estimate of the steering matrix
% theta_estm = angle(A_estm(2,:));

% P = S*S'./N;
% [Us,Ls] = eig(A*P*A');
% [Ls ,Is]  = sort(diag(Ls),1,'descend');   %Find K largest eigenvalues
% Us = Us (:,Is);
% Us = Us(:,1:K);
% Ls = Ls(1:K);
% Ls = diag(Ls);
% 
% Rn = Noise*Noise'./N;
% [Un,Ln] = eig(Rn);
% [Ln ,In]  = sort(diag(Ln),1,'descend');   
% Un = Un (:,In);
% Un = Un(:,K+1:end);
% Ln = Ln(K+1:end);
% Ln = diag(Ln);
% 
% R_hat = A*P*A' + Rn;
% [Q_sum,D_sum] = eig(R_hat);
% [D_sum,I_sum] = sort(diag(D_sum),1,'descend');
% Q_sum = Q_sum(:,I_sum);
% Q_sum_s = Q_sum(:,1:K);
% Q_sum_n = Q_sum(:,K+1:end);
% 
% disp(['norm(R) = ',num2str(norm(R,"fro"))])
% disp(['norm(R - R_hat) = ',num2str(norm(R-R_hat,"fro"))])