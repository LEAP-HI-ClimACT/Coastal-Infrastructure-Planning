% rewards matrix for finite horizon (terminal costs- expected returns over
% next 20 years) 
% one terminal state
% carbon costs accumulated and with increasing scc over time-
% non-stationary rewards model

% Action Rewards - States (end state) dependent = R(s',s,a):

% load('P.mat');
horizon = 40;
core_states = 5544;
states = core_states*rate;
total_states = (states)*4 + 1;% one terminal state extra

[row,col] = find(P);
fst = sparse(row, col, ones(1,length(row)), states, states);
n_systems = 4;

tic
% Action 1 - do nothing
% load('R_dn_two_dikes_d.mat');
R_rate= cell(1,n_systems);

%system A
R_core = cell(1, horizon);
R_t = R_dn{1,1};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end

R1=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R1(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R1(i+5544:i+11087,i+5544:i+11087)= 0.0;
R1 = R1.*fst;

R_rate{1,1} = sparse(R1);

%system B
R_core = cell(1, horizon);
R_t = R_dn{1,2};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end
R2=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R2(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R2(i+5544:i+11087,i+5544:i+11087)= 0.0;
R2 = R2.*fst;

R_rate{1,2} = sparse(R2);

%system C
R_core = cell(1, horizon);
R_t = R_dn{1,3};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end
R3=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R3(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R3(i+5544:i+11087,i+5544:i+11087)= 0.0;
R3 = R3.*fst;

R_rate{1,3} = sparse(R3);

%system D
R_core = cell(1, horizon);
R_t = R_dn{1,4};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end
R4=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R4(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R4(i+5544:i+11087,i+5544:i+11087)= 0.0;
R4 = R4.*fst;

R_rate{1,4} = sparse(R4);

% Full do nothing
Full = sparse(total_states, total_states);
for i= 1:states:total_states-1
    Full(i:i+(states-1), i:i+(states-1)) = R_rate{1,floor(i/states)+1};
end
R_f1 = R_f_dn{1,1};
R_f2 = R_f_dn{1,2};
R_f3 = R_f_dn{1,3};
R_f4 = R_f_dn{1,4};
Full(5544*39+1: states, total_states) = R_f1';
Full(states+5544*39+1: states*2, total_states) = R_f2';
Full(2*states+5544*39+1: states*3, total_states) = R_f3';
Full(3*states+5544*39+1: states*4, total_states) = R_f4';

Full(total_states, total_states) = 0.00;
inp.R{1,1} = sparse(Full);
r1 = toc

% Action 2 - construct dike 1 only
% load('R_1_two_dikes_d.mat');
R_rate= cell(1,n_systems);

%system A
R_core = cell(1, horizon);
R_t = R_1{1,1};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end
R1=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R1(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R1(i+5544:i+11087,i+5544:i+11087)= 0.0;
R1 = R1.*fst;

R_rate{1,1} = sparse(R1);

%system B
R_core = cell(1, horizon);
R_t = R_1{1,2};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end
R2=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R2(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R2(i+5544:i+11087,i+5544:i+11087)= 0.0;
R2 = R2.*fst;

R_rate{1,2} = sparse(R2);

%system C
R_core = cell(1, horizon);
R_t = R_1{1,3};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end
R3=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R3(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R3(i+5544:i+11087,i+5544:i+11087)= 0.0;
R3 = R3.*fst;

R_rate{1,3} = sparse(R3);

%system D
R_core = cell(1, horizon);
R_t = R_1{1,4};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end
R4=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R4(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R4(i+5544:i+11087,i+5544:i+11087)= 0.0;
R4 = R4.*fst;

R_rate{1,4} = sparse(R4);

% Full action 1 reward matrix
Full = sparse(total_states, total_states);

Full(1:1+(states-1),1+(states):1+2*(states)-1) = R_rate{1,1};
Full(1*states+1:2*states,1+1*(states):1+2*(states)-1) = R_rate{1,2};
Full(2*states+1:3*states,1+3*states:1+4*(states)-1) = R_rate{1,3};
Full(3*states+1:4*states,1+3*states:1+4*(states)-1) = R_rate{1,4};

R_f1 = R_f_1{1,1};
R_f2 = R_f_1{1,2};
R_f3 = R_f_1{1,3};
R_f4 = R_f_1{1,4};
Full(5544*39+1: 5544*40, total_states) = R_f1';
Full(states+5544*39+1: states*2, total_states) = R_f2';
Full(2*states+5544*39+1: states*3, total_states) = R_f3';
Full(3*states+5544*39+1: states*4, total_states) = R_f4';

Full(total_states, total_states) = 0.00;

inp.R{1,2} = sparse(Full);
r2 = toc

% Action 3 - construct dike 2 only
% load('R_2_two_dikes_d.mat');
R_rate= cell(1,n_systems);

%system A
R_core = cell(1, horizon);
R_t = R_2{1,1};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end
R1=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R1(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R1(i+5544:i+11087,i+5544:i+11087)= 0.0;
R1 = R1.*fst;

R_rate{1,1} = sparse(R1);

%system B
R_core = cell(1, horizon);
R_t = R_2{1,2};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end
R2=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R2(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R2(i+5544:i+11087,i+5544:i+11087)= 0.0;
R2 = R2.*fst;

R_rate{1,2} = sparse(R2);

%system C
R_core = cell(1, horizon);
R_t = R_2{1,3};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end
R3=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R3(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R3(i+5544:i+11087,i+5544:i+11087)= 0.0;
R3 = R3.*fst;

R_rate{1,3} = sparse(R3);

%system D
R_core = cell(1, horizon);
R_t = R_2{1,4};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end
R4=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R4(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R4(i+5544:i+11087,i+5544:i+11087)= 0.0;
R4 = R4.*fst;

R_rate{1,4} = sparse(R4);

% Full action 2 matrix
Full = sparse(total_states, total_states);

Full(1:states,1+2*(states):1+3*(states)-1) = R_rate{1,1};
Full(1*states+1:2*states,1+3*states:1+4*(states)-1) = R_rate{1,2};
Full(2*states+1:3*states,1+2*(states):1+3*(states)-1) = R_rate{1,3};
Full(3*states+1:4*states,1+3*states:1+4*(states)-1) = R_rate{1,4};

R_f1 = R_f_2{1,1};
R_f2 = R_f_2{1,2};
R_f3 = R_f_2{1,3};
R_f4 = R_f_2{1,4};
Full(5544*39+1: 5544*40, total_states) = R_f1';
Full(states+5544*39+1: states*2, total_states) = R_f2';
Full(2*states+5544*39+1: states*3, total_states) = R_f3';
Full(3*states+5544*39+1: states*4, total_states) = R_f4';

Full(total_states, total_states) = 0.00;

inp.R{1,3} = sparse(Full);
r3 = toc

% Action 4 - construct both dikes
% load('R_12_two_dikes_d.mat');
R_rate= cell(1,n_systems);

%system A
R_core = cell(1, horizon);
R_t = R_12{1,1};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end
R1=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R1(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R1(i+5544:i+11087,i+5544:i+11087)= 0.0;
R1 = R1.*fst;

R_rate{1,1} = sparse(R1);

%system B
R_core = cell(1, horizon);
R_t = R_12{1,2};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end
R2=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R2(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R2(i+5544:i+11087,i+5544:i+11087)= 0.0;
R2 = R2.*fst;

R_rate{1,2} = sparse(R2);

%system C
R_core = cell(1, horizon);
R_t = R_12{1,3};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end
R3=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R3(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R3(i+5544:i+11087,i+5544:i+11087)= 0.0;
R3 = R3.*fst;

R_rate{1,3} = sparse(R3);

%system D
R_core = cell(1, horizon);
R_t = R_12{1,4};
parfor t = 1:horizon
    R = R_t{1,t};
    dum = R_t{1,t};
    for i=1:core_states-1
        R=[R; dum];
    end
    R_core{1,t} = R;
end
R4=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    R4(i:i+5543,i+5544:i+11087)=R_core{krate};
end
current = i
R4(i+5544:i+11087,i+5544:i+11087)= 0.0;
R4 = R4.*fst;

R_rate{1,4} = sparse(R4);

% Full action 3 reward matrix
Full = sparse(total_states, total_states);

Full(1:1+(states-1),1+3*(states):1+4*(states)-1) = R_rate{1,1};
Full(1+1*states:2*states,1+3*(states):1+4*(states)-1) = R_rate{1,2};
Full(1+2*states:3*states,1+3*(states):1+4*(states)-1) = R_rate{1,3};
Full(1+3*states:4*states,1+3*(states):1+4*(states)-1) = R_rate{1,4};

R_f1 = R_f_12{1,1};
R_f2 = R_f_12{1,2};
R_f3 = R_f_12{1,3};
R_f4 = R_f_12{1,4};
Full(5544*39+1: 5544*40, total_states) = R_f1';
Full(states+5544*39+1: states*2, total_states) = R_f2';
Full(2*states+5544*39+1: states*3, total_states) = R_f3';
Full(3*states+5544*39+1: states*4, total_states) = R_f4';

Full(total_states, total_states) = 0.00;
inp.R{1,4} = sparse(Full);
r4= toc


%% call rewards_fp to get the expected costs from the starting state
rewards_fp; %to convert R(s',s,a) to R(s,a) for fast parser R: %d: %d: * %lf