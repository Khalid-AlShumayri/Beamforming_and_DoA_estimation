close all
clear all 
clc 

M =32;      % sensors
K = 3;      % number of signals (sources)
N = 100;     % number of observations
d = 0.5;    % Distance between elements in wavelengths
Pn = .9;    % Noise power
sig_pr = 0.9.*ones(1,K);    % signals' power
DoA = [-30 20 60];

L = floor(M / 2);  % Length of subarrays
% L =M;
snr = min(sig_pr)/Pn;
Resolution = (2/(L*snr))*180/pi;
iter = 1;
j = 1;

% Compute steering Matrix and vectors
angles=(-90:1:90);       % for grid search

% far-field assumption 
A  = generate_steering_matrix(M,d,DoA);
a_srch = exp(-1i*2*pi*d*(0:L-1)'*sin([angles(:).']*pi/180));

S = diag(sqrt(sig_pr./2))*(randn(K,N)+1j*randn(K,N));

Noise = sqrt(Pn/2)*( randn(M,N) + 1j*randn(M,N) );
X = A*S + Noise;

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

% Uncorrelated signals
[Q ,D]  = eig(R);
[D ,I]  = sort(diag(D),1,'descend');   %Find K largest eigenvalues
D = diag(D);
Q = Q (:,I);       % Sort the eigenvectors to put signal eigenvectors first 
Qs = Q (:,1:K);       % Get the signal eigenvectors
Qn = Q(:,K+1:L);    % Get the noise eigenvectors

% Grid search
srch_music = zeros(1, length(angles));
for i = 1:length(angles)
    srch_music(i) = abs(1 / (a_srch(:, i)' * Qn * Qn' * a_srch(:, i)));
end

figure
plot(angles,10*log10(abs(srch_music))); hold on; grid on;
xline(DoA, '-.', 'LineWidth', .9);

title(['SNR = ', num2str(10 * log10(min(sig_pr) / Pn)), 'dB : Resolution = ', num2str(Resolution), ' degrees']);
xlabel('Angle [degrees]');
ylabel('Spatial Spectrum [dB]')


%%% ---- Peak detection ---- %%%
[peaks, locs] = findpeaks(srch_music);
threshold = mean(srch_music) + std(srch_music);
significant_peaks = angles(locs(peaks > threshold));
detected_DoAs = sort(significant_peaks);
DoA_music = zeros(1, K);
DoA_music(1, 1:length(detected_DoAs)) = detected_DoAs;
%%% ------------------------ %%%

yline(10*log10(threshold), '--', 'LineWidth', .9);
legend('MUSIC', 'Actual DoAs', 'Threshold');


RMSE_srch = calc_rmse(DoA,DoA_music)
%% Root MUSIC

phi = exp(1i*d*2*pi*sin(DoA*pi/180));
root_Qn = roots_cols(flipud(Qn));       % flipping the signal up side down

rlr = real(root_Qn);
mgr = imag(root_Qn);
data = [rlr(:), mgr(:)];

labels = dbscan_clustering(data,.3,floor(3)); % points within the cluster is M-K

figure
scatter(real(root_Qn),imag(root_Qn),'o'); hold on;
scatter(real(phi),imag(phi),'x','black')

figure
scatter(real(root_Qn(labels==-1)),imag(root_Qn(labels==-1)),'o',color=[.4 .2 .6]); hold on;
scatter(real(root_Qn(labels==0)),imag(root_Qn(labels==0)),'o',color=[.2 .4 .6]); hold on;
scatter(real(phi),imag(phi),'x','black')% points within the cluster is M-K
%% Root MUSIC (fft approach)

% N_dft = L;
% spectrum_root =  sum( abs(fft(Qn,N_dft)).^2 ,2);
% % trial = fftshift( abs(fft(sum(Qn,2).^2,3*M)) );      % we have to take the fft of each column indiv. 
% 
% % % plot
% freq = [-N_dft/2:1:N_dft/2-1]/N_dft;
% index = asin(freq/d)*180/pi;
% 
% % % Peak detection and error
%  
% % Step 1: Initial peak detection
% [peaks2, locs2] = findpeaks(1./spectrum_root);
% 
% % Step 2: Determine a dynamic threshold
% threshold2 = mean(1./spectrum_root);
% 
% % Step 3: Identify significant peaks
% ind2 = sort(index( locs2(peaks2 > threshold2) ));
% 
% AoA_fft = zeros(1,K);
% AoA_fft(1,1:length(ind2)) = ind2
% DoA
% RMSE_root = sqrt(sum((DoA-AoA_fft).^2)/K)
% 
% figure
% plot(index,1./spectrum_root); hold on; grid on;
% yline(threshold2,'--')
% xline(DoA,'--')
% title('FFT')
% xlabel('Angles [Degrees]')
% ylabel('PSD')
% 
% %% Validation
% 
% [Un,Ln] = eig(Noise*Noise'./N);
% P = S*S'/N;
% [Us,Ls] = eig(A*P*A');