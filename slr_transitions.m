%% calculate transition probabilities for sea level rise states

action_set = {'do nothing', 'construct only_dike1','construct only_dike2','construct both dike1_and_dike2'};
system_set = {'no dikes', 'only dike1','only dike2','both dikes'};

% transitions same for every action and for every system 

tic;
N = 131;
del_t = 1;
n_states= 77; 
n_samples = 100*1e6;%100 million
load('simulations_slr.mat'); % simulations of SLR trajectories
values(:,1) = 7; % from the observed last year data: needed?

states_real = zeros(n_samples, N);
parfor i = 1:N
    for j = 1:n_samples
        states_real(j,i) = convert_to_states(values(j,i));
    end
end

n_s = zeros(n_states, N-1);
states = 1:n_states;
parfor i = 1:N
    for j = 1:n_states
        n_s(j,i) = sum(states_real(:,i) == states(j));
    end
end


t_s_slr = zeros(n_states, n_states, N);
t_s_slr(:,:,1) = eye(77,77);
n_t = zeros(n_states, n_states, N-1);
n_years = 80;

for i = 2:n_years
    for j = 1:n_states
        n_s_1 = n_s(j,i-1);
        if n_s_1 ~= 0
            for k = 1:n_states
                n_t(k,j,i) = 0;
                for l = 1:n_samples
                    if states_real(l,i) == states(k) && states_real(l,i-1) == states(j)
                       n_t(k,j,i) =  n_t(k,j,i)+1;
                    end
                end
            t_s_slr(k,j,i)= n_t(k,j,i)/n_s_1;
            end
        end
    end
    
end
toc;

%% smoothing of transition probabilities
% average of 5 years

t_s_avg = zeros(n_states, n_states, N);
t_s_avg(:,:,1) = eye(77,77);

for i = 2:N-5
    for j = 1:n_states
        for k = 1:n_states
            t_s_avg(k,j,i) = mean(t_s_slr(k,j,i:i+5));
        end
    end

end

for i = 127:131
    for j = 1:n_states
        for k = 1:n_states
            t_s_avg(k,j,i) = mean(t_s_slr(k,j,i:end));
        end
    end
end
