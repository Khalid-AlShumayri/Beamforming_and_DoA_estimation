close all
clear all 
clc 

% Comparison between three grid-search approaches
% 1- Conventional Beamformer
% 2- Capon's Beamformer
% 3- MUSIC algorithm

M = 6;      % sensor
K = 3;      % number of signals (sources)
N = 10;     % number of observations
d = 0.5;    % Distance between elements in wavelengths
Pn = .1;    % Noise power
sig_pr = [.9 .7 .5];    % signals' power
DoA = [-30 20 35];
iter = 3;      % This is to eliminate variability between different 
                % simulation runs

% Compute steering Matrix and vectors

angles=(-90:0.05:90);       % for grid search

% far-field assumed
a1 = exp(-1i*2*pi*d*(0:M-1)'*sin([angles(:).']*pi/180));  
A  = generate_steering_matrix(M,d,DoA);

spectrum_music = zeros(iter,length(angles));
spectrum_bf = zeros(iter,length(angles));
spectrum_cap = zeros(iter,length(angles));

for j=1:1:iter

    S = diag(sqrt(sig_pr))*randn(K,N);
    Noise = sqrt(Pn)*randn(M,N);
    X = A*S + Noise;
    R = X*X'/N ;
    
    [Q ,D]  = eig(R);
    [D ,I]  = sort(diag(D),1,'descend');   %Find K largest eigenvalues
    D = diag(D);
    Q = Q (:,I);       % Sort the eigenvectors to put signal eigenvectors first 
    
    Qs = Q (:,K);       % Get the signal eigenvectors
    Qn = Q(:,K+1:M);    % Get the noise eigenvectors
    
    
    for k=1:length(angles)
        % Compute MUSIC spectrum
        search_music(k) = abs( 1/(a1(:,k)'*Qn*Qn'*a1(:,k)) ); 
        search_bf(k) = abs(a1(:,k)'*R*a1(:,k)/(a1(:,k)'*a1(:,k)));
        search_cap(k) = abs(20/( a1(:,k)'*inv(R)*a1(:,k) ));
    end
    spectrum_music(j,:) = search_music;
    spectrum_bf(j,:)    = search_bf;
    spectrum_cap(j,:) = search_cap;
end

spectrum_music = mean(spectrum_music);
spectrum_bf    = mean(spectrum_bf);
spectrum_cap = mean(spectrum_cap);

plot(angles,spectrum_music,Color=[0 .2 .9]); hold on; grid on;
plot(angles,spectrum_bf,Color=[0.8, 0.3, 0.1]);
plot(angles,spectrum_cap,Color=[0.5, 0.2, 0.6]);

xline(DoA,'-.',LineWidth=.9)
legend('MUSIC','BF','Capon BF','Actual DoAs')
title(['Spatial Spectrum : SNR = ',num2str(10*log10(max(sig_pr)/Pn)),'dB'])
xlabel('Angle in degrees')

figure
plot(angles,spectrum_bf,Color=[0.8, 0.3, 0.1]);hold on; grid on;
plot(angles,spectrum_cap,Color=[0.5, 0.2, 0.6]);
xline(DoA,'-.',LineWidth=.9)

legend('BF','Capon BF','Actual DoAs')
title(['Spatial Spectrum : SNR = ',num2str(10*log10(max(sig_pr)/Pn)),'dB'])
xlabel('Angle in degrees')
%% Analytical estimation

A_est = X*S'*inv(S*S');      % estimate of the steering matrix