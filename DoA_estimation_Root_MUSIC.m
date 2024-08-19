close all
clear all 
clc 

% Parameters
DoA = [0 40 69];
M = 8;      % sensors
K = 3;      % number of signals (sources)
N = 20;     % number of observations (snapshots)
d = 0.5;    % Distance between elements in wavelengths
snr = 10;   %dB
snr_linear = 10^(snr/10);
sig_pr = 1 .* ones(1, K);    % signals' power
max_sig_pr = max(sig_pr);
Pn = max_sig_pr*10^(-snr/10);    % Noise power
num_iterations = 100;  % Number of Monte Carlo iterations
stepsize = 1;
Res = (2/M) * 180 / pi;

L = floor(M / 2);  % Length of subarrays (spatial smoothing)

% Generate correlated sources
P_ch =  reshape([1 .6 .1 .6  1 .1 .1 .1 1],K,K);
Lr = chol(P_ch, "lower");

% Far-field assumption 
angles = -90:stepsize:90;  % for grid search 
a = exp(-1i * 2 * pi * d * (0:M-1)' * sin([angles(:).'] * pi / 180));
a_sps = exp(-1i * 2 * pi * d * (0:L-1)' * sin([angles(:).'] * pi / 180));
A = generate_steering_matrix(M, d, DoA);

% Initialize results
rmse_values         = zeros(1, num_iterations);
rmse_values_ch      = zeros(1, num_iterations);
rmse_values_ch_sps  = zeros(1, num_iterations);  

% Monte Carlo simulation
for iter = 1:num_iterations
    % Generate steering matrix and signals
    
    S = diag(sqrt(sig_pr ./ 2)) * (randn(K, N) + 1j * randn(K, N));
    Noise = sqrt(Pn / 2) * (randn(M, N) + 1j * randn(M, N));
    X = A * S + Noise;
    R = X*X'./N;
    
    % Generate correlated sources
    S_ch = Lr * S;
    X_ch = A * S_ch + Noise;
    R_ch = X_ch * X_ch' ./ N;

    % Compute spatially smoothed covariance matrix
    R_sps = zeros(L, L);
    R_ch_sps = zeros(L, L);
    for k = 1:(M - L + 1)
        x_sub = X(k:k + L - 1, :);
        x_ch = X_ch(k:k+L-1,:);
        R_sps = R_sps + (x_sub * x_sub') ./ (N + M - L + 1);
        R_ch_sps = R_ch_sps + (x_ch * x_ch') ./ (N + M - L + 1);
    end

    % MUSIC Algorithm for uncorrelated signals
    [Q, D] = eig(R);
    [D, I] = sort(diag(D), 1, 'descend');
    D = diag(D);
    Q = Q(:, I);
    Qs = Q(:, 1:K);
    Qn = Q(:, K+1:end);

    % MUSIC Algorithm for correlated signals
    [Q, D] = eig(R_ch);
    [D, I] = sort(diag(D), 1, 'descend');
    D = diag(D);
    Q = Q(:, I);
    Qn_ch = Q(:, K+1:end);

    % MUSIC Algorithm for correlated signals
    [Q, D] = eig(R_ch_sps);
    [D, I] = sort(diag(D), 1, 'descend');
    D = diag(D);
    Q = Q(:, I);
    Qn_ch_sps = Q(:, K+1:L);


    % Grid search
    srch_music    = zeros(1, length(angles));
    srch_ch     = zeros(1, length(angles));
    srch_ch_sps = zeros(1, length(angles));

    for i = 1:length(angles)
        srch_music(i)    = abs(1 / (a(:, i)' * Qn * Qn' * a(:, i)));
        srch_ch_sps(i) = abs(1 / (a_sps(:, i)' * Qn_ch_sps * Qn_ch_sps' * a_sps(:, i)));
        srch_ch(i)     = abs(1 / (a(:, i)' * Qn_ch* Qn_ch' * a(:, i)));
    end

    %%% ---- Peak detection ---- %%%
    [peaks, locs1] = findpeaks(srch_music);
    threshold = mean(srch_music) + 0.5*std(srch_music);
    significant_peaks = angles(locs1(peaks > threshold));
    detected_DoAs = sort(significant_peaks);
    DoA_music = zeros(1, K);
    DoA_music(1, 1:length(detected_DoAs)) = detected_DoAs;

    [peaks, locs2] = findpeaks(srch_ch);
    threshold = mean(srch_ch);
    significant_peaks = angles(locs2(peaks > threshold));
    detected_DoAs = sort(significant_peaks);
    DoA_ch = zeros(1, K);
    DoA_ch(1, 1:length(detected_DoAs)) = detected_DoAs;

    [peaks, locs3] = findpeaks(srch_ch_sps);
    threshold = mean(srch_ch_sps) ;
    significant_peaks = angles(locs3(peaks > threshold));
    detected_DoAs = sort(significant_peaks);
    DoA_ch_sps = zeros(1, K);
    DoA_ch_sps(1, 1:length(detected_DoAs)) = detected_DoAs;

    %%% ------------------------ %%%

    % Calculate RMSE for this iteration
    rmse_values(iter)        = calc_rmse(DoA,DoA_music);
    rmse_values_ch(iter)     = calc_rmse(DoA,DoA_ch);
    rmse_values_ch_sps(iter) = calc_rmse(DoA,DoA_ch_sps);
end

% Calculate and display mean RMSE
rmse_sps    = mean(rmse_values);
rmse_ch     = mean(rmse_values_ch);
rmse_ch_sps = mean(rmse_values_ch_sps);

disp(['Mean RMSE (uncorrelated signals): ', num2str(rmse_sps)]);
disp(['Mean RMSE (correlated signals): ', num2str(rmse_ch)]);
disp(['Mean RMSE (correlated signals + SPS): ', num2str(rmse_ch_sps)]);
disp(['Step-size = ',num2str(stepsize),' degrees']);
disp(['Number of iterations = ',num2str(num_iterations)])

% Plot example of the last iteration
figure;
hold on; grid on;
h1 = plot(angles, srch_music, 'Color', [0 .2 .9]); 
h2 = plot(angles, srch_ch); 
h3 = plot(angles,srch_ch_sps);
h4 = xline(DoA, '-.', 'LineWidth', .9,Color='black');
% h4 = yline(threshold, '--', 'LineWidth', .9);
legend([h1 h2 h3 h4(1)],{'Uncorrelated', 'Correlated','Correlated + SPS','Actual DoAs'},Location="best");
title(['MUSIC algorithm : SNR = ', num2str(10 * log10(min(sig_pr) / Pn)), 'dB : Resolution = ', num2str(Res), ' degrees']);
xlabel('Angle [degrees]');
ylabel('Spatial Spectrum')



%% ----------------------------------------------------------
%  ------------------------ Root MUSIC ----------------------
%  ----------------------------------------------------------

rmse_dbscan = zeros(1,iter);

% Monte-Carlo 
for iter = 1:num_iterations
    % Generate steering matrix and signals
    
    S = diag(sqrt(sig_pr ./ 2)) * (randn(K, N) + 1j * randn(K, N));
    Noise = sqrt(Pn / 2) * (randn(M, N) + 1j * randn(M, N));
    X = A * S + Noise;
    R = X*X'./N;

    % Generate correlated sources
    S_ch = Lr * S;
    X_ch = A * S_ch + Noise;
    R_ch = X_ch * X_ch' ./ N;

    % Compute spatially smoothed covariance matrix
    R_sps = zeros(L, L);
    R_ch_sps = zeros(L, L);
    for k = 1:(M - L + 1)
        x_sub = X(k:k + L - 1, :);
        x_ch = X_ch(k:k+L-1,:);
        R_sps = R_sps + (x_sub * x_sub') ./ (N + M - L + 1);
        R_ch_sps = R_ch_sps + (x_ch * x_ch') ./ (N + M - L + 1);
    end

    % MUSIC Algorithm for uncorrelated signals
    [Q, D] = eig(R);
    [D, I] = sort(diag(D), 1, 'descend');
    D = diag(D);
    Q = Q(:, I);
    Qs = Q(:, 1:K);
    Qn = Q(:, K+1:end);


    % MUSIC Algorithm for correlated signals
    [Q, D] = eig(R_ch);
    [D, I] = sort(diag(D), 1, 'descend');
    D = diag(D);
    Q = Q(:, I);
    Qn_ch = Q(:, K+1:end);
    Qn_root = Qn_ch;

    % MUSIC Algorithm for correlated signals
    [Q, D] = eig(R_ch_sps);
    [D, I] = sort(diag(D), 1, 'descend');
    D = diag(D);
    Q = Q(:, I);
    Qn_ch_sps = Q(:, K+1:L);


    phi = exp(1i*d*2*pi*sin(DoA*pi/180));
    % flipping the signal up side down and find the roots of each column
    roots_Qn = roots_cols(flipud(Qn_root));       
    
    rlr = real(roots_Qn);
    rlr = rlr(:);
    imgr = imag(roots_Qn);
    imgr = imgr(:);
    data = [rlr(:), imgr(:)];
    
    % labels = dbscan_clustering(data,.05,floor(L-K-3)); % points within the cluster is M-K
    [label,corepts] = dbscan(data,round(log(M)/log(2) - 1)/M,floor(M-K-2)); % distance ~ 1/M
    

    % plot the roots of the columns of the matrix Un 
%     figure
%     scatter(data(:,1),data(:,2),'o'); hold on;
%     scatter(real(phi),imag(phi),'x','red',linewidth=1.2)
%     legend('Roots','True AoAs')
%     title('Roots of matrix Q_n columns')
%     
%     figure
%     scatter(real(phi),imag(phi),'x','red',linewidth=1.2);hold on;
%     gscatter(data(:,1),data(:,2),label)
%     title('Before applying KNN')

    
    % Calculate the number of clusters
    uniqueClusters = unique(label(label > 0));  % Exclude noise (clusters labeled as 0)
    numClusters = length(uniqueClusters);

    % Split large clusters until the number of clusters matches knownNumClusters
    while numClusters < K
        disp(['Number of clusters: ',num2str(numClusters)])
        % Find the largest cluster
        clusterSizes = histcounts(label(label>0), 'BinWidth',0.5);
        [~, maxClusterIdx] = max(clusterSizes);
        maxClusterLabel = uniqueClusters(maxClusterIdx);
        
        % Extract points belonging to the largest cluster
        clusterPoints = data(label == maxClusterLabel, :);
        idx_max = find(label==maxClusterLabel);
        % Split the largest cluster using K-means
        k = 2; % Start by splitting into 2 clusters, adjust as needed
        [newIdx, ~] = kmeans(clusterPoints, k, 'Replicates', 5);
    
        % Update the original clustering result
        maxLabel = max(label); % Find the maximum label in the current idx
        for i = 1:k
            label(idx_max(newIdx == i)) = maxLabel + i;
        end
        
        % Recalculate the number of clusters
        uniqueClusters = unique(label(label > 0));  % Exclude noise (clusters labeled as 0)
        numClusters = length(uniqueClusters);
    end

%     figure
%     p1 = scatter(real(phi),imag(phi),'x','red',linewidth=1.2); hold on;
%     gscatter(data(:,1),data(:,2),label)
%     % legend([p],{'True AoAs'})
%     title('After applying KNN to break big clusters')
    
    % Calculate RMSE for this iteration
    phi_dbscan = center_median(data,label,K);   % find the median of each cluster
    DoA_dbscan = asin(angle(phi_dbscan)/(2*pi*d))*180/pi;
    rmse_dbscan(iter) = calc_rmse(DoA,DoA_dbscan);

end

% plot of the last iteration
figure
p1 = scatter(real(phi),imag(phi),'x','red',linewidth=1.2); hold on;
gscatter(data(:,1),data(:,2),label)
% legend([p],{'True AoAs'})
title('After applying KNN to break big clusters')


% Calculate and display mean RMSE
rmse_root_dbscan    = mean(rmse_dbscan);
% rmse_ch     = mean(rmse_values_ch);
% rmse_ch_sps = mean(rmse_values_ch_sps);

disp(['Mean RMSE for Root-MUSIC(uncorrelated signals): ', num2str(rmse_root_dbscan)]);
% disp(['Mean RMSE (correlated signals): ', num2str(rmse_ch)]);
% disp(['Mean RMSE (correlated signals + SPS): ', num2str(rmse_ch_sps)]);
% disp(['Step-size = ',num2str(stepsize),' degrees']);
% disp(['Number of iterations = ',num2str(num_iterations)])


% 
% 
% % plot the roots of the columns of the matrix Un 
% figure
% scatter(data(:,1),data(:,2),'o'); hold on;
% scatter(real(phi),imag(phi),'x','red',linewidth=1.2)
% legend('Roots','True AoAs')
% title('Roots of matrix Q_n columns')
% 
% figure
% scatter(real(phi),imag(phi),'x','red',linewidth=1.2);hold on;
% gscatter(data(:,1),data(:,2),label)
% title('Before applying KNN')
% 
% 
% Calculate the number of clusters
% uniqueClusters = unique(label(label > 0));  % Exclude noise (clusters labeled as 0)
% numClusters = length(uniqueClusters);
% 
% Split large clusters until the number of clusters matches knownNumClusters
% while numClusters < K
%     disp(['Number of clusters: ',num2str(numClusters)])
%     Find the largest cluster
%     clusterSizes = histcounts(label(label>0), 'BinWidth',0.5);
%     [~, maxClusterIdx] = max(clusterSizes);
%     maxClusterLabel = uniqueClusters(maxClusterIdx);
%     
%     Extract points belonging to the largest cluster
%     clusterPoints = data(label == maxClusterLabel, :);
%     idx_max = find(label==maxClusterLabel);
%     Split the largest cluster using K-means
%     k = 2; % Start by splitting into 2 clusters, adjust as needed
%     [newIdx, ~] = kmeans(clusterPoints, k, 'Replicates', 5);
% 
%     Update the original clustering result
%     maxLabel = max(label); % Find the maximum label in the current idx
%     for i = 1:k
%         label(idx_max(newIdx == i)) = maxLabel + i;
%     end
%     
%     Recalculate the number of clusters
%     uniqueClusters = unique(label(label > 0));  % Exclude noise (clusters labeled as 0)
%     numClusters = length(uniqueClusters);
% end
% 
% figure
% p2 = scatter(real(phi),imag(phi),'x','red',linewidth=1.2); hold on;
% gscatter(data(:,1),data(:,2),label)
% legend([p],{'True AoAs'})
% title('After applying KNN to break big clusters')
% 
% phi_dbscan = center_median(data,label,K);   % find the median of each cluster
% DoA_dbscan = asin(angle(phi_dbscan)/(2*pi*d))*180/pi;
% rmse_dbscan = calc_rmse(DoA,DoA_dbscan)


% points within the cluster is M-K
% figure
% % gscatter(rlr,imgr,label)
% h1 = scatter(data(label==-1,1),data(label==-1,2)); hold on;
% h2 = scatter(real(phi),imag(phi),'x','red') ; 
% h3 = scatter(real(phi_dbscan),imag(phi_dbscan),'x','black');
% legend([h1,h2,h3],{'Noise','Ground Truth','Predicted'})
% 
% title(['Roots of the columns of the matrix Un'])
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

% Rn = Noise*Noise'/N;
% Rv = inv(Rn);
% [Qn,Dn] = eig(Rn);
% P = S*S'/N;
% [Qs,Ds] = eig(A*P*A');

% test = Rv - Rv*A*inv(eye(k,k) + P*A'*Rv*A)*P*A'*Rv;