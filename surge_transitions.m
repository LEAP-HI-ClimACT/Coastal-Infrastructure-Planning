% calculation of transition probabilities for storm surge states (GEV)

% GEV parameters for storms in metres- done for this model for now!!

% baseline parameters: model 2: 2.6296 m
mu=0.936; %location
sigma=0.206; %scale
k=0.232; %shape

figure(1)
hold on
x =linspace(-1,7,10000000);%m
plot(x*100, gevpdf(x, k, sigma, mu), 'linewidth', 2)
title('GEV distribution for storm surge')
xlabel('surge values (cm)')
ylabel('pdf')


figure(2)
plot(x*100, cdf('gev',x,k, sigma, mu),'linewidth',2)
title('GEV distribution for storm surge')
xlabel('surge values (cm)')
ylabel('cdf')


surge_values = [0:10:700]; %cm
states_gev = [1:72];
n_states_gev = length(states_gev);
prob_states = zeros(1,length(states_gev));

% gev transitions : same for every action and system
for i = 1:length(surge_values)-1
    prob_states(i+1) = cdf('gev',surge_values(i+1)/100,k,sigma,mu) - cdf('gev',surge_values(i)/100,k,sigma,mu);
end
prob_states(72) = 1 - cdf('gev',7, k, sigma, mu);
prob_states(1) = cdf('gev', 0, k, sigma, mu);

figure
plot(states_gev, prob_states, 'linewidth', 2)
xlabel('storm surge states')
ylabel('probability of state')
title('probability of storm surge states')

t_s_surge = repmat(prob_states, n_states_gev, 1);
t_s_surge = t_s_surge';%to match the format of transition slr
