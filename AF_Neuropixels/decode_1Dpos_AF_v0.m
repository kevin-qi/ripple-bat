function p_dec = decode_1Dpos_AF_v0(n_vec,t_bin,p_x,f_x)
%% FUNCTION FOR DECODING THE POSTERIOR PROBABILITIES OF A LINEARIZED POSITION, GIVEN SPIKES, PRIOR AND PLACE MAPS (A. Forli, Feb 2024)

% Reference articles: 
% Zhang, K., Ginzburg, I., McNaughton, B. L., & Sejnowski, T. J. (1998). 
% Interpreting neuronal population activity by reconstruction: unified framework with application to hippocampal place cells. 
% Silva, D., Feng, T., & Foster, D. J. (2015). 
% Trajectory events across hippocampal place cells require previous experience. 

% --- INPUTS ---
% n_vec:    vector of number of spikes from all the cells withing the considered time bin
% t_bin:    duration of the time bin (s)
% p_x:      prior probabilty of being within that spatial bin
% f_x:      matrix of firing rate as a function of the spatial bin for all the cells (place fields)
%---------------

%--- OUTPUTS ---
% p_dec:    posterior probabilities of being within that spatial bin
%---------------

%=== Extract number of cells and number of positional bins
n_cells = numel(n_vec);
n_pos = size(f_x,1);

%=== Make sure inputs have the correct size
if ~iscolumn(n_vec),  n_vec=n_vec';         end
if size(f_x,2)~=n_cells, f_x = f_x';        end
if ~iscolumn(p_x),       p_x = p_x';        end

%=== Calculate decoded position
%p_dec = p_x.*(prod(f_x.^(n_vec'),2).*exp(-t_bin*sum(f_x,2)));   % Full version
p_dec = prod(f_x.^(n_vec'),2).*exp(-t_bin*sum(f_x,2));          % Uniform prior


%=== OLD: Loop decoding across positions, replaced by the compact formula above
% p_dec = zeros(size(p_x));   % Initialize decoded position
% for i=1:n_pos
%     p_dec(i) = p_x(i)*prod(f_x(i,:).^(n_vec'))*exp(-t_bin*sum(f_x(i,:)));
% end

%=== Normalize
p_dec = p_dec./sum(p_dec);

end

