close all
clear all 
clc 

M =16;      % sensors
K=3;
K_root = K;
alpha = 0.1; % threshold to classify the eigenvalues
N = 25;     % number of observations
d = 0.5;    % Distance between elements in wavelengths
Pn = 1;     % Noise power
sig_pr = 0.9;
sig_corr = sig_pr*ones(1,K);    % signals' power
DoA = [-30 20 25];
stepsize = 0.01;

% L = floor(M / 2);  % Length of subarrays
L = M;

snr = sig_pr/Pn;
Resolution = (2/(sqrt(N)*L*snr))*180/pi;    % Resolution of the MUSIC algorithm. Source: Chatgpt!
iter = 1;
j = 1;

% Compute steering Matrix and vectors
angles=(-90:stepsize:90);       % for grid search

% far-field assumption 
A  = generate_steering_matrix(M,d,DoA);
a_srch = exp(-1i*2*pi*d*(0:L-1)'*sin([angles(:).']*pi/180));

S = diag(sqrt(sig_corr./2))*(randn(K,N)+1j*randn(K,N));

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
R = R/N;

% Uncorrelated signals
[Q ,D]  = eig(R);
D = D./trace(D);
[D ,I]  = sort(diag(D),1,'descend');   %Find K largest eigenvalues
K = length(D(D>alpha));
D = diag(D);
Q = Q (:,I);       % Sort the eigenvectors to put signal eigenvectors first 
Us = Q (:,1:K);       % Get the signal eigenvectors
Un = Q(:,K+1:L);    % Get the noise eigenvectors
Un_root = Q(:,K+1:L);

% Grid search
srch_music = zeros(1, length(angles));
for i = 1:length(angles)
    srch_music(i) = abs(1 / (a_srch(:, i)' *Un*Un' * a_srch(:, i)));
end

figure
plot(angles,10*log10(abs(srch_music))); hold on; grid on;
xline(DoA, '-.', 'LineWidth', .9);

title(['step-size = ',num2str(stepsize),': SNR = ', num2str(10 * log10(min(sig_pr) / Pn)), 'dB : Resolution = ', num2str(Resolution), ' degrees']);
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
RMSE_srch = calc_rmse(DoA,DoA_music)

yline(10*log10(threshold), '--', 'LineWidth', .9);
legend('MUSIC', 'Actual DoAs', 'Threshold');



%% Root MUSIC

phi = exp(1i*d*2*pi*sin(DoA*pi/180));
root_Un = roots_cols(flipud(Un_root));       % flipping the signal up side down

rlr = real(root_Un);
imgr = imag(root_Un);
data = [rlr(:), imgr(:)];

% labels = dbscan_clustering(data,.05,floor(L-K-3)); % points within the cluster is M-K
[label,corepts] = dbscan(data,3/M,floor(M-K_root-1)); % distance ~ 1/M

% plot the roots of the columns of the matrix Un 
figure
scatter(data(:,1),data(:,2),'o'); hold on;
scatter(real(phi),imag(phi),'x','red',linewidth=1.2)

% plot the clustered roots
% CLR = ['blue','orange'];
figure
% h1 = scatter(data(label==-1,1),data(label==-1,2),'o',color=[.4 .2 .6]); hold on;
scatter(real(phi),imag(phi),'x','red',linewidth=1.2);hold on;
gscatter(data(:,1),data(:,2),label); 
% for i=1:K
%     scatter(data(label==i,1),data(label==i,2),'o');
% end

%%
% Step 2: Calculate the number of clusters
uniqueClusters = unique(label(label > 0));  % Exclude noise (clusters labeled as 0)
numClusters = length(uniqueClusters);

% Step 3: Split large clusters until the number of clusters matches knownNumClusters
while numClusters < K
    % Find the largest cluster
    clusterSizes = histcounts(label, uniqueClusters);
    [~, maxClusterIdx] = max(clusterSizes);
    maxClusterLabel = uniqueClusters(maxClusterIdx);
    
    % Extract points belonging to the largest cluster
    clusterPoints = data(label == maxClusterLabel, :);
    
    % Split the largest cluster using K-means
    k = 2; % Start by splitting into 2 clusters, adjust as needed
    [newIdx, ~] = kmeans(clusterPoints, k, 'Replicates', 5);
    
    % Update the original clustering result
    maxLabel = max(label); % Find the maximum label in the current idx
    for i = 1:k
        label(label == maxClusterLabel & newIdx == i) = maxLabel + i;
    end
    
    % Recalculate the number of clusters
    uniqueClusters = unique(label(label > 0));  % Exclude noise (clusters labeled as 0)
    numClusters = length(uniqueClusters);
end

figure
gscatter(rlr,imgr,label)
%%

phi_dbscan = center_median(data,label,K);
DoA_dbscan = asin(angle(phi_dbscan)/(2*pi*d))*180/pi;
rmse_dbscan = calc_rmse(DoA,DoA_dbscan)

% points within the cluster is M-K
h2 = scatter(real(phi),imag(phi),'x','red') ;
h3 = scatter(real(phi_dbscan),imag(phi_dbscan),'x','black');
legend([h1,h2,h3],{'Noise','Ground Truth','Predicted'})

title(['Roots of the columns of the matrix Un'])
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
%% Validation

Rn = Noise*Noise'/N;
Rv = inv(Rn);
[Qn,Dn] = eig(Rn);
P = S*S'/N;
[Qs,Ds] = eig(A*P*A');

test = Rv - Rv*A*inv(eye(k,k) + P*A'*Rv*A)*P*A'*Rv;