function TMI = Theta_mod_idx_AF_v0(spikes,time_bin,theta_range, max_lag)
%% Function for calculating the theta modulation index

%=== Autocorrelogram
[cross_corr,bin_centers] = cross_correlogram_AF_v0(spikes,spikes,max_lag,time_bin);     % Calculate autocorrelogram
cross_corr = cross_corr-mean(cross_corr);                                               % Subtract mean


n_fft = 2^nextpow2(numel(cross_corr));                                                  % Optimal number of points for FFT
Y_fft = fft(cross_corr,n_fft);                                                          % Calculate fft
P_fft = abs(Y_fft/n_fft).^2;                                                            % Power density at all frequences (positive and negative)
f_fft = (1/time_bin)*(0:n_fft/2)/n_fft;                                                 % Fs is the sampling frequency
PSD = P_fft(1:n_fft/2+1);                                                               % Power spectral density
PSD_sm = smoothdata(PSD,'movmedian',n_fft/(1/time_bin));                                % Smooth
PSD_th = PSD_sm;    PSD_th([1:knnsearch(f_fft',theta_range(1)),knnsearch(f_fft',theta_range(2)):numel(f_fft)])=0;   % Keep only the theta range
[~,idx_tmp] = max(PSD_th);                                                              % Find maximum in the theta range
n_1Hz = round(1/mean(diff(f_fft)));                                                     % Number of samples corresponding to 1Hz interval
if idx_tmp==1
    TMI = NaN;
else
   TMI = mean(PSD_th(idx_tmp+[-n_1Hz:n_1Hz]))/mean(PSD_th(1:knnsearch(f_fft',50)));% Theta modulation index 
end
    
end

