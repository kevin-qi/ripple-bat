%% Simple code for simulating the generation of theta sequences

clr;
t_sim = 5;              % Seconds to simulate, 5
f_hi = 1e3;             % Sampling freq, 1e3
f_lo = 200;             % Sampling for bins, 200
dt_sim = 1/f_lo;        % Time duration (bins)
N = 50;                 % Number of neurons, 50
t = [0:1/f_hi:t_sim];   % Time
T = numel(t);           % Number of time samples
f_wbt = 8;              % Wingbeat frequency, 8
w = 2*pi*f_wbt;         % Omega (wingbeat)
r = 2*pi*10;            % Omega (intrinsic)
isc = cos(r*t);         % Wingbeat signal
delta = 0.0;            % Contribution of the intrinsic oscillator
bin_TMI = 0.010;        % Time bin for autocorrelation
max_lag_TMI = 0.5;      % Max lag for autocorrelation
N_spatial_bins = 50;    % Number of spatial bins for decoding

pk_fir = 300;          % Proportional to peak firing rate, 150
beta = 0.7;             % Contribution of the wigbeat, 0.8
alpha = 0.0;            % Asymmetry factor 0.9
sigma = 0.3;            % Duration of the firing field, 0.4
dt_adapt = 0.1;        % Adaptation time 0.08  
min_spk_adapt = 3;      % Min spikes for adaptation, 2
wbd = 0;                % Delay wingbeat signal

wbt = cos(w*(t-wbd));   % Wingbeat signal
N_skw = 4*f_hi;         % Asymmetry size (1 narrow, 4 wide)
krnel = flip(circshift(gampdf(linspace(0,100*alpha,N_skw),2,2),N_skw/2));

%=== Initialize relevant vectors
time_bins = [t(1):1/f_lo:t(T)];
time_ctrs = time_bins(1:end-1)+0.5/f_lo;
N_bins = numel(time_ctrs);
rate = zeros(N,T);
t0 = linspace(t(1)+0.1,t(end)-0.1,N);

%=== Generate poisson rate from the product of wigbeat and place fields
s = cell(N,1);
tic;
for i=1:N
    
    adapted = zeros(1,T);
    %=== Gaussian firing rate
    rate(i,:) = exp(-(t-t0(i)).^2/sigma^2);
    
    if alpha>0
        base=zeros(1,T);    base(knnsearch(t',t0(i)))=1;
        rate(i,:) = normalize(conv(base,krnel,'same'),'range');
    end
    
    for j=1:N_bins
        smpl = knnsearch_fast_AF_v0(t',time_ctrs(j),0.01);
        %lambda = rate(i,smpl)*(1-beta*wbt(max(smpl-round(0.012*f_hi),1)))*(1-delta*isc(smpl));
        lambda = rate(i,smpl)*(1-beta*wbt(smpl))*(1-delta*isc(smpl));
        if ~adapted(j)
            Num_sim = poissrnd(pk_fir*lambda*dt_sim);
            s{i,1} = [s{i,1}; (time_bins(j)+rand(Num_sim,1)*dt_sim)];
            if Num_sim>=min_spk_adapt
                n_adapt = round(f_lo*dt_adapt);
                adapted(j:min(j+n_adapt,T)) = 1;
            end
        else
            Num_sim = 0;
        end
            
    end
    s{i,1} = sort(s{i,1});
end
toc;

%=== Calculate wigbeat phase of each cell
wbt_phase_hi = wrapTo2Pi(angle(hilbert(wbt)))-pi;
spk_phase = cellfun(@(x) wbt_phase_hi(ceil(x*f_hi)),s,'UniformOutput',false);

%=== Plot decoded activity
t_dec_bin_dur = 0.030;                                                      % Time bin duration for decoding (def 25ms)
t_dec_bin_ovl = t_dec_bin_dur-0.005;                                        % Time bin overlap for decoding (t_bin_dur-0.005)
smooth_f = [1 .3];                                                          % [space bin,time bin]
prc_lim = [1 98.5];
n_dec_bins = floor((t(T)-t_dec_bin_dur)/(t_dec_bin_dur-t_dec_bin_ovl))+1;   % Number of time bins
st_times = (0:n_dec_bins-1)*(t_dec_bin_dur - t_dec_bin_ovl);                % Start times
ed_times = st_times + t_dec_bin_dur;                                        % End times
ct_times = st_times + t_dec_bin_dur/2;                                      % Center times
t_d = zeros(1, 2 * n_dec_bins);                                             % All times (initialize)
t_d(1:2:end) = st_times;    t_d(2:2:end) = ed_times;                        % All times (define)
n_vector = zeros(1,N);                                                      % Initialize vector containing the number of spikes fired by each place cell
p_x = ones(1,N_spatial_bins);                                               % Probability of being on a given bin of the flight cluster
f_x = zeros(numel(p_x),N);                                                  % Vector with the firing rates of place cells in a given bin of the flight cluster
for nc=1:N, f_x(:,nc) = interp1(rate(nc,:)',linspace(1,T,N_spatial_bins)');end
p_dec_flight = zeros(size(p_x,2),n_dec_bins);
p_dec_shifted = zeros(size(p_x,2),n_dec_bins);
spk_dsty = zeros(1,n_dec_bins);
bin_pos_real = round(linspace(round(t0(1)/t(T)*N_spatial_bins),round(t0(end)/t(T)*N_spatial_bins),n_dec_bins));

%=== Loop across bins
for i_dec = 1:n_dec_bins
    for nc=1:N,n_vector(nc) = histcounts(s{nc,1},[t_d(2*i_dec-1),t_d(2*i_dec)]);end
    p_x = double(p_x>-1);   % Uniform Prior
    p_dec_flight(:,i_dec) = decode_1Dpos_AF_v0(n_vector,t_dec_bin_dur,p_x,f_x);
    p_dec_shifted(:,i_dec) = circshift(p_dec_flight(:,i_dec)  ,-bin_pos_real(i_dec)+round(N_spatial_bins/2));
    spk_dsty(i_dec) = sum(n_vector);
end

%=== PLOTTING ==========================================================================================

%=== Show place fields
figure('units','normalized','outerposition',[0 .7 .5 .3]);
cmap = jet(N);
for i=1:2:N
    plot(t,rate(i,:),'Color',cmap(i,:));    hold on;
end

%=== Plot firing fields
figure('units','normalized','outerposition',[0 .1 .5 .3]);
for nc= 1:N
    plot(s{nc}, nc*ones(size(s{nc})), 'k|','MarkerSize', max(1,round(N*0.1)));   hold on;          % Raster for each cell
end
hold on;    plot(t,normalize(wbt,'range',[0 N]));


%=== Plot phase vs position
figure('units','normalized','outerposition',[0 .5 1 .5]);
tiledlayout(7,10,'TileSpacing','compact');
for nc = 1:min(N,100)
    nexttile;   scatter([s{nc,1};s{nc,1}],[spk_phase{nc,1}';spk_phase{nc,1}'+2*pi],6,'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.5);    hold on;
    %yticks(pi*[-1 0 1 2 3]);  yticklabels({'-180', '0', '180', '360', '540'});  
    axis square;    xticks([]); yticks([]);
end

%=== Plot decoded trace
figure('units','normalized','outerposition',[0 .4 .5 .3]);
imagesc([t(1) t(T)],[-1 1],imgaussfilt(p_dec_flight,smooth_f),prctile(p_dec_flight, prc_lim,'all')');
colormap(flipud(gray));     hold on;    plot(t,wbt,'r');
ylabel('Spatial bin'); axis off;  set(gca,'YDir','normal'); xlabel('Temporal bin');

%=== Plot autocorrelograms
% figure('units','normalized','outerposition',[0 0 1 .4]);
% tiledlayout(10,10,'TileSpacing','tight');
ac = [];
for nc = 1:N
    [ac_TMI_tmp,ac_TMI_bins] = cross_correlogram_AF_v0(s{nc,1},s{nc,1},max_lag_TMI,bin_TMI);
%     nexttile;   area(ac_TMI_bins,ac_TMI_tmp,'FaceColor','k','EdgeColor','none','FaceAlpha',0.5);
%     yticks([]); hold on;    plot(repmat(1/f_wbt*[-3:3],2,1),repmat(ylim,7,1)','k--');   xlim([0 max_lag_TMI]);
    ac = [ac; ac_TMI_tmp'];
end

%=== Average autocorrelogram
figure('units','normalized','outerposition',[.5 .4 .15 .3]);
area(ac_TMI_bins,mean(ac,'omitnan'),'FaceColor','r','EdgeColor','none','FaceAlpha',0.5); yticks([]); xlabel('Time lag (s)'); hold on;
plot(repmat(1/f_wbt*[-3:3],2,1),repmat(ylim,7,1)','k--');
xlim([0 max_lag_TMI]);

%=== Params and initialize relevant variables
n_reshape = 30;                 % Default 30
binned_wbt = interp1(t,wbt,ct_times);
wbt_phase = wrapTo2Pi(angle(hilbert(binned_wbt)))-pi;

%=== Cut at wingbeat maxima
single_SWP = table();
counter=1;
warning('off');
bin_size_1D = 0.15;

zero_phs_idx = find(wbt_phase(1:end-1).* wbt_phase(2:end)<0 & diff(wbt_phase)<0);  % Segment based on wingbeat
%plot(wbt_phase);   hold on; stem(zero_phs_idx, ones(size(zero_phs_idx )));  plot(wbt);
sweep_strt = zero_phs_idx(1:end-1); sweep_stop = zero_phs_idx(2:end);
spt_bin_ids = [1:N_spatial_bins]';

for ss=1:numel(sweep_strt)
    
    single_SWP.sft_posterior(counter) = {imgaussfilt(p_dec_shifted(:,sweep_strt(ss):sweep_stop(ss)),smooth_f)};
    single_SWP.rsp_posterior(counter) = {imresize(single_SWP.sft_posterior{counter,1},[size(single_SWP.sft_posterior{counter,1},1),n_reshape])};
    
    single_SWP.wbt(counter) = {binned_wbt(sweep_strt(ss):sweep_stop(ss))};
    single_SWP.rsz_wbt(counter) = {interp1(single_SWP.wbt{counter,1},linspace(1,numel(single_SWP.wbt{counter,1}),n_reshape)')};
    
    single_SWP.spk_dsty(counter) = {spk_dsty(sweep_strt(ss):sweep_stop(ss))};
    single_SWP.rsz_spk_dsty(counter) = {interp1(single_SWP.spk_dsty{counter,1},linspace(1,numel(single_SWP.spk_dsty{counter,1}),n_reshape)')};
    
    
    single_SWP.wbt_phase(counter) = {wbt_phase(sweep_strt(ss):sweep_stop(ss))};
    single_SWP.rsz_wbt_phase(counter) = {interp1(single_SWP.wbt_phase{counter,1},linspace(1,numel(single_SWP.wbt_phase{counter,1}),n_reshape)')};
    
    single_SWP.mean_spk_dsty(counter) = mean(single_SWP.spk_dsty{counter,1});
    cnt_mass = spt_bin_ids'*single_SWP.sft_posterior{counter,1};
    single_SWP.med_jmp_distance(counter) = median(abs(diff(cnt_mass)));
    single_SWP.est_dist(counter) = {interp1((cnt_mass-N_spatial_bins/2)*bin_size_1D,linspace(1,numel(cnt_mass),n_reshape)')};
    
    counter=counter+1;
end
warning('on');

%=== Extract subset using defined criteria (exclude flight tails, epochs of low firing and flat sweeps)
SWP_sst = single_SWP;
%SWP_sst = single_SWP(single_SWP.mean_spk_dsty>prctile(single_SWP.mean_spk_dsty,10) & single_SWP.med_jmp_distance>0.1,:);
%SWP_sst = single_SWP(single_SWP.mean_spk_dsty>prctile(single_SWP.mean_spk_dsty,5),:);

N_sweeps = size(SWP_sst,1);

%=== Calculate averages and STDs
est_dist = zeros(size(SWP_sst.est_dist{1,1},1),size(SWP_sst,1));
rsz_spk_dsty = zeros(size(SWP_sst.rsz_spk_dsty{1,1},1),size(SWP_sst,1));
rsz_wbt = zeros(size(SWP_sst.rsz_wbt{1,1},1),size(SWP_sst,1));
rsz_wbt_phase = zeros(size(SWP_sst.rsz_wbt{1,1},1),size(SWP_sst,1));
for i=1:N_sweeps
    est_dist(:,i) = SWP_sst.est_dist{i,1};
    rsz_wbt(:,i) = SWP_sst.rsz_wbt{i,1};
    rsz_spk_dsty(:,i) = SWP_sst.rsz_spk_dsty{i,1};
    rsz_wbt_phase(:,i) = SWP_sst.rsz_wbt_phase{i,1};
end
avg_est_dist = mean(est_dist,2);            sem_est_dist = std(est_dist,[],2)/sqrt(N_sweeps);
avg_rsz_wbt = mean(rsz_wbt,2);              sem_rsz_wbt = std(rsz_wbt,[],2)/sqrt(N_sweeps);
avg_rsz_wbt_phase = mean(rsz_wbt_phase,2);  sem_rsz_wbt_phase = std(rsz_wbt_phase,[],2)/sqrt(N_sweeps);
avg_rsz_spk_dsty = mean(rsz_spk_dsty,2);    sem_rsz_spk_dsty = std(rsz_spk_dsty,[],2)/sqrt(N_sweeps);

%=== Show averages
figure('units','normalized','outerposition',[.8 .4 .17 .3]);
tiledlayout(2,3,'TileSpacing','compact');
nexttile;
plotWinterval_AF_v0(linspace(-180,180,n_reshape),avg_est_dist,avg_est_dist-sem_est_dist,avg_est_dist+sem_est_dist,'k');
hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
xticks([-180 0 180]);    title('Average Decoding Error');    %yticks([]);
nexttile;
plotWinterval_AF_v0(linspace(-180,180,n_reshape),avg_rsz_spk_dsty,avg_rsz_spk_dsty-sem_rsz_spk_dsty,avg_rsz_spk_dsty+sem_rsz_spk_dsty,'r');
hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
xticks([-180 0 180]);    title('Average Spike Density');    %yticks([]);
nexttile;
plotWinterval_AF_v0(linspace(-180,180,n_reshape),avg_rsz_wbt_phase,avg_rsz_wbt_phase-sem_rsz_wbt_phase,avg_rsz_wbt_phase+sem_rsz_wbt_phase,'m');
hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
xticks([-180 0 180]);    title('Wingbeat Phase');    %yticks([]);
nexttile;
plotWinterval_AF_v0(linspace(-180,180,n_reshape),avg_rsz_wbt,avg_rsz_wbt-sem_rsz_wbt,avg_rsz_wbt+sem_rsz_wbt,'k');
hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
xticks([-180 0 180]);    ylabel('Accelerometer');    %yticks([]);
nexttile;
plotWinterval_AF_v0(linspace(-180,180,n_reshape),avg_rsz_wbt,avg_rsz_wbt-sem_rsz_wbt,avg_rsz_wbt+sem_rsz_wbt,'k');
hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
xticks([-180 0 180]);    ylabel('Accelerometer');    %yticks([]);
nexttile;
plotWinterval_AF_v0(linspace(-180,180,n_reshape),avg_rsz_wbt,avg_rsz_wbt-sem_rsz_wbt,avg_rsz_wbt+sem_rsz_wbt,'k');
hold on;    plot(0*[1 1],ylim,'k--'); xlim('tight');
xticks([-180 0 180]);    ylabel('Accelerometer');    %yticks([]);
