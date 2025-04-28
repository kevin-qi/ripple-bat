function [PSD_sm,f_fft] = Theta_mod_idx_AF_v1(spikes,time_bin, max_lag)
%% Function for calculating the theta modulation index

%=== Autocorrelogram
[cross_corr,~] = cross_correlogram_AF_v0(spikes,spikes,max_lag,time_bin);     % Calculate autocorrelogram

n_fft = 2^nextpow2(numel(cross_corr));                                                  % Optimal number of points for FFT
Y_fft = fft(cross_corr,n_fft);                                                          % Calculate fft
P_fft = abs(Y_fft/n_fft).^2;                                                            % Power density at all frequences (positive and negative)
f_fft = (1/time_bin)*(0:n_fft/2)/n_fft;                                                 % Fs is the sampling frequency
PSD = P_fft(1:n_fft/2+1);                                                               % Power spectral density
PSD_sm = smoothdata(PSD,'movmedian',n_fft/(1/time_bin));                                % Smooth
    
end

