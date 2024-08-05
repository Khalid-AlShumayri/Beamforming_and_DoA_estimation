close all
clear all 
clc 

M = 20;      % sensors
K = 3;      % number of signals (sources)
N = 5;     % number of observations
d = 0.5;    % Distance between elements in wavelengths
Pn = .9;    % Noise power
sig_pr = 0.9.*ones(1,K);    % signals' power
DoA = [-30 20 60];
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

%% Root MUSIC

phi = exp(1i*d*2*pi*sin(DoA*pi/180));
root_Qn = roots_cols(flipud(Qn));
root_mag= abs(root_Qn);
root_arg= angle(root_Qn);
[A,indexMatrix] = sort(root_mag);

figure
polarplot(angle(root_Qn),abs(root_Qn),"o"); hold on;
polarplot(angle(phi),abs(phi),'x',Color='black',LineWidth=1.2)

%%% ---------- This block sort the mag. and arg of the roots -------- %%%

[numRows, numCols] = size(root_arg); %  Get the number of rows and columns

% Create row and column indices
[rowIdx, colIdx] = ndgrid(1:numRows, 1:numCols);

% Convert the row indices based on the index matrix
rowIdx = indexMatrix(:);

% Convert the column indices to match the structure of the matrix
colIdx = colIdx(:);

% Calculate the linear indices
linearIndices = sub2ind(size(root_arg), rowIdx, colIdx);

% Reorder the matrix using the linear indices
B = reshape(root_arg(linearIndices), numRows, numCols);
%
%%% ----------------------------------------------------------------- %%%



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

% figure
% plot(index,1./spectrum_root); hold on; grid on;
% yline(threshold2,'--')
% xline(DoA,'--')
% title('FFT')
% xlabel('Angles [Degrees]')
% ylabel('PSD')