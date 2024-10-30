% calculate transition probabilities for combined states (SLR + storm
% surge)

action_set = {'do nothing', 'construct only_dike1','construct only_dike2','construct both dike1_and_dike2'};
system_set = {'no dikes', 'only dike1','only dike2','both dikes'};

% transitions same for every action and for every system 

horizon1 = 80; %from 2020 to 2099
slr_states = 77;
surge_states = 72;
total_states = slr_states*surge_states;

t_total1 = cell(1,horizon1);
tic

for t =1:horizon1
    Ts = zeros(total_states, total_states);
    s=1;
    for i =1:slr_states
        for j=1:surge_states
            T_m=zeros(slr_states,surge_states);
            for k=1:slr_states
                for l=1:surge_states
                    T_m(k,l)=t_s_avg(i,k,t)*t_s_surge(j,l);
                end
            end
            Ts(s,:)=reshape(T_m',1,[]);
            s=s+1;
        end
    end
    t_total1{1,t} = Ts';% transpose done here
end

time = toc

