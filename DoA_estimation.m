close all
clear all 
clc 

% Comparison between three grid-search DoA approaches
% 1- Conventional Beamformer
% 2- Capon's Beamformer
% 3- spectrum-MUSIC algorithm

M = 5;      % sensors
K = 3;      % number of signals (sources)
N = 10;      % number of observations
d = 0.5;    % Distance between elements in wavelengths
Pn = .09;    % Noise power
sig_pr = .9.*ones(1,K);    % signals' power
snr = min(sig_pr/Pn);
DoA = [-30 -20 35];
Resolution = 2*180/(pi*snr*M*sqrt(N));
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
R_cap = R + .05*trace(R)*eye(M,M);
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
    srch_cap(i) = abs(1/( a1(:,i)'*inv(R_cap)*a1(:,i) ));
%     search_cap_cn = 
end
spec_music(j,:)  = srch_music;
spec_bf(j,:)    = srch_bf;
spec_cap(j,:) = srch_cap;


plot(angles,spec_music,Color=[0 .2 .9]); hold on; grid on;
xline(DoA,'-.',LineWidth=.9)

legend('MUSIC','Actual DoAs')
title(['Spatial Spectrum : SNR = ',num2str(10*log10(snr)),'dB : Resolution = ',num2str(Resolution),' degrees'])
xlabel('Angle [degrees]')

figure 
subplot(211)
plot(angles,spec_bf,Color=[0.8, 0.3, 0.1]); grid on; hold on;
xline(DoA,'-.',LineWidth=.9)
title('Conventional beamformer')
xlabel('Angle [degrees]')

subplot(212)
plot(angles,spec_cap,Color=[0.5, 0.2, 0.6]);grid on; hold on;
xline(DoA,'-.',LineWidth=.9)
title('Capon''s beamformer')
xlabel('Angle [degrees]')


% Step 1: Initial peak detection
[peaks, locs] = findpeaks(spec_music);

% Step 2: Determine a dynamic threshold
threshold = mean(spec_music) + 1.5*std(spec_music);


% Step 3: Identify significant peaks
ind = sort(angles( locs(peaks > threshold) ));
AoA_music = zeros(1,K);
AoA_music(1,1:length(ind)) = ind

DoA = sort(DoA)
RMSE = sqrt(sum((DoA-AoA_music).^2)/K)

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

%% Root MUSIC

phi_xdomain = angle(exp(1i*d*pi*sin(DoA*pi/180)));
root_Qn = roots_matrix(Qn);
root_mag= abs(root_Qn);
root_arg= angle(root_Qn);
[root_mag,indexMatrix] = sort(root_mag);

% Create a linear index array for each column
% linearIndices = indexMatrix + (0:M:(M*(M-K-1))); %works for square matrix

%%% ---------- This block sort the mag. and arg of the roots ----- %%%

%  Get the number of rows and columns
[numRows, numCols] = size(root_arg);

% Create row and column indices
[rowIdx, colIdx] = ndgrid(1:numRows, 1:numCols);

% Convert the row indices based on the index matrix
rowIdx = indexMatrix(:);

% Convert the column indices to match the structure of the matrix
colIdx = colIdx(:);

% Calculate the linear indices
linearIndices = sub2ind(size(root_arg), rowIdx, colIdx);

% Reorder the matrix using the linear indices
root_arg = reshape(root_arg(linearIndices), numRows, numCols)

%%% ------------------------------------------------ %%%

%% Root MUSIC (fft approach)

N_dft = 5*M;
spectrum_root =  sum( abs( fft(Qn,N_dft)).^2 ,2) ;
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

figure
subplot(211)
plot(angles,spec_music,Color=[0 .2 .9]); hold on; grid on;
xline(DoA,'-.',LineWidth=.9)
legend('MUSIC','Actual DoAs')
title(['Spatial Spectrum : SNR = ',num2str(10*log10(snr)),'dB : Resolution = ',num2str(Resolution),' degrees'])
xlabel('Angle [degrees]')

subplot(212)
plot(index,1./spectrum_root); hold on; grid on;
yline(threshold2,'--')
xline(DoA,'--')
title('FFT')
xlabel('Angles [Degrees]')
ylabel('PSD')