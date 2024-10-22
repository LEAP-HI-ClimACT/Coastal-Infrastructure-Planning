% finite horizon formulation
%% transition matrix

rate = 40;
core_states = 5544;
load('t_total_norm_100_gev2.mat'); %the calculated transitions
jMap=0;
for i=2:rate
jMap=jMap+1;
Phat1{jMap}=t_total_norm{i};
end

T_z = 0;
s_z = [];
s_v = [];

for t=1:rate-1
T = Phat1{t};
T(T<=1.0e-5) = 0.0;
T = round(T,10);

% renormalization
  for i = 1:core_states
      s_v = [s_v ; sum(T(i,:))];
      if sum(T(i,:))~= 0.0
        T(i,:) = T(i,:)/sum(T(i,:));
      else
        T_z = T_z+1;
        s_z =[s_z; i];
        T(i,i) = 1.0;
      end
      
  end
Pcore{1,t}= round(T,10);
end

%
T_nz = 0;
T1_z = 0;
s1_z = [];

T_sum_1 = 0;
T_sum_not1 = 0;
s_v1 = [];

for t = 1:rate-1
    T_n = Pcore{1,t};
    for i = 1:core_states
        % check for zero sum of outgoing transitions
        if (sum(T_n(i,:))- 0.0) >= 1e-10
        T_nz = T_nz + 1;
        else
        T1_z = T1_z+1;
        s1_z =[s1_z; i];
        end
        % check for sum of outgoing transitions to be = 1
        s_v1 = [s_v1; sum(T_n(i,:))];
        if round(sum(T_n(i,:)),5) ~= 1.00
            T_sum_not1 = T_sum_not1 + 1;
        else
            T_sum_1 = T_sum_1 + 1;
        end
            
    end
end       

states = core_states*rate;
% total_states = (states)*4 + 1;% one terminal state extra
P=sparse(states,states);
krate=0;
for i=1:5544:(states - 11087)%5544*2 - 1 = 11087%3*5544+1 = states - 11087 = 16633
    krate=krate+1;
    P(i:i+5543,i+5544:i+11087)=Pcore{krate};
end
current = i
P(i+5544:i+11087,i+5544:i+11087)= 0.0;
T_dn{1,1} = sparse(P);

%% Full do_nothing transition matrix
% 
tic
T_dn{1,1} = sparse(P);
total_states = (states)*4+1;
Full = sparse(total_states, total_states);

for i= 1:states:total_states-1
    Full(i:i+(states-1), i:i+(states-1)) = T_dn{1,1};
end

Full(5544*39+1: states, total_states) = 1.00;
Full(states+5544*39+1: states*2, total_states) = 1.00;
Full(2*states+5544*39+1: states*3, total_states) = 1.00;
Full(3*states+5544*39+1: states*4, total_states) = 1.00;

Full(total_states, total_states) = 1.00;
inp.T{1,1} = sparse(Full);
%
t1 = toc

%% Full construct only dike 1 transition matrix
% 
Full = sparse(total_states, total_states);
Full(1:1+(states-1),1+(states):1+2*(states)-1) = T_dn{1,1};
Full(1*states+1:2*states,1+1*(states):1+2*(states)-1) = T_dn{1,1};
Full(2*states+1:3*states,1+3*states:1+4*(states)-1) = T_dn{1,1};
Full(3*states+1:4*states,1+3*states:1+4*(states)-1) = T_dn{1,1};

Full(5544*39+1: 5544*40, total_states) = 1.00;
Full(states+5544*39+1: states*2, total_states) = 1.00;
Full(2*states+5544*39+1: states*3, total_states) = 1.00;
Full(3*states+5544*39+1: states*4, total_states) = 1.00;

Full(total_states, total_states) = 1.00;
inp.T{1,2}= sparse(Full);
% 
t2 = toc
%% Full construct only dike 2 transition matrix
% 

Full = sparse(total_states, total_states);
Full(1:states,1+2*(states):1+3*(states)-1) = T_dn{1,1};
Full(1*states+1:2*states,1+3*states:1+4*(states)-1) = T_dn{1,1};
Full(2*states+1:3*states,1+2*(states):1+3*(states)-1) = T_dn{1,1};
Full(3*states+1:4*states,1+3*states:1+4*(states)-1) = T_dn{1,1};

Full(5544*39+1: 5544*40, total_states) = 1.00;
Full(states+5544*39+1: states*2, total_states) = 1.00;
Full(2*states+5544*39+1: states*3, total_states) = 1.00;
Full(3*states+5544*39+1: states*4, total_states) = 1.00;

Full(total_states, total_states) = 1.00;

inp.T{1,3}= sparse(Full);
% 
t3 = toc
%% Full construct both dike 1 and dike 2 transition matrix
% 

Full = sparse(total_states, total_states);
Full(1:1+(states-1),1+3*(states):1+4*(states)-1) = T_dn{1,1};
Full(1+1*states:2*states,1+3*(states):1+4*(states)-1) = T_dn{1,1};
Full(1+2*states:3*states,1+3*(states):1+4*(states)-1) = T_dn{1,1};
Full(1+3*states:4*states,1+3*(states):1+4*(states)-1) = T_dn{1,1};

Full(5544*39+1: 5544*40, total_states) = 1.00;
Full(states+5544*39+1: states*2, total_states) = 1.00;
Full(2*states+5544*39+1: states*3, total_states) = 1.00;
Full(3*states+5544*39+1: states*4, total_states) = 1.00;

Full(total_states, total_states) = 1.00;

inp.T{1,4}= sparse(Full);

t4 = toc