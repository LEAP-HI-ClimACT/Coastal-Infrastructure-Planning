% observation likelihoods for two climate models
% observations: SLR

load('n_s_slr_245.mat') % simulations of SLR with SSP245
n_s_245 = n_s;

load('n_s_slr_585.mat')% simulations of SLR with SSP585
n_s_585 = n_s;

slr_states = size(n_s_245, 1);
interval = 1;
total_horizon = 60;
rate = total_horizon/interval; % 8
core_states = 5544;
states = core_states*rate;
total_states = (states)*4+1;
total_states_f = total_states*2; %2 climate models
n_obs = 77*2 + 1;
n_samples = 100*1e6;

prob_slr_245 = zeros(slr_states,rate);
prob_slr_585 = zeros(slr_states,rate);

rate_ = 1;
for i=2:interval:total_horizon
    prob_slr_245(:,rate_) = n_s_245(:,i)/n_samples;
    prob_slr_585(:,rate_) = n_s_585(:,i)/n_samples;  
    rate_ = rate_ + 1;
end

n_models = 2;
obs_models = zeros(slr_states, n_models);

for i = 1 : rate
    for j = 1:slr_states
        p_obs1 = prob_slr_245(j,i);
        p_obs2 = prob_slr_585(j,i);
        if (p_obs1 + p_obs2) < 1e-8
            p_obs1_n = 0;
            p_obs2_n = 0;
        else
            p_obs1_n = p_obs1/(p_obs1+p_obs2);
            p_obs2_n = p_obs2/(p_obs1+p_obs2);
        end
        obs_models(j, :) = [p_obs1_n p_obs2_n];
    end
    p_obs_models{1,i} = obs_models;
    
end

%% expand the observables to SLR_states*n_models + 1(terminal)

for t = 1:rate
    p_obs_models_ = p_obs_models{t};
% Define the size of the expanded matrix
    expanded_matrix_size = [size(p_obs_models_, 1), size(p_obs_models_, 1) * size(p_obs_models_, 2) +1];

    % Initialize the expanded matrix
    expanded_matrix = zeros(expanded_matrix_size);

    % Populate the expanded matrix
    for i = 1:size(p_obs_models_, 1)
        expanded_matrix(i,i) = p_obs_models_(i,1);
        expanded_matrix(i, 77+i) = p_obs_models_(i,2);
    end
    p_obs_models_slr{t} = expanded_matrix;
end

%% expansion for surge states
% Define the number of times to repeat each row
num_surge_states = 72;
for t = 1:rate
    p_obs_ = p_obs_models_slr{t};
% Preallocate space for the expanded matrix
    expanded_matrix = zeros(size(p_obs_, 1) * num_surge_states, size(p_obs_, 2));

% Repeat each row and stack them vertically
    for i = 1:size(p_obs_, 1)
        start_index = (i - 1) * num_surge_states + 1;
        end_index = i * num_surge_states;
        expanded_matrix(start_index:end_index, :) = repmat(p_obs_(i, :), num_surge_states, 1);
    end
    p_obs_surge_models{t} = expanded_matrix;

end
    
%%  expansion for time
p_obs_time = [];

for i = 1:rate
    p_obs_t = p_obs_surge_models{t};
    p_obs_time = [p_obs_time; p_obs_t];
end

%% expansion for systems

n_systems = 4;
p_obs_sys = [];
for i = 1:n_systems
    p_obs_sys = [p_obs_sys; p_obs_time];
end

% 1st terminal state
p_obs_sys(total_states,n_obs) = 1.0;
%% expansion for models

n_models = 2;
p_obs_models_final = [];

for i =1:n_models
    p_obs_models_final = [p_obs_models_final; p_obs_sys];
end

%% 
n_actions = 4;
for i= 1:n_actions
    inp.O{1,i} = p_obs_models_final;
end
