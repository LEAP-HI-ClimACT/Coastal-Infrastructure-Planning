% generate the pomdp input file for fast parser

tic
interval = 10;
total_horizon = 40;
rate = total_horizon/interval; % 8
core_states = 5544;
states = core_states*rate;
total_states = (states)*4+1;
total_states_f = total_states*2; %2 climate models
% name = 'pomdp_5_years_245_585';
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extens = '.pomdp';
% filename = join([name extens],'');
head01 = fopen('heading_pomdp_245_585_10y_model2_high.txt','w');
% head01 = fopen(filename,'w'); %opening the file

% parsing;

[rpa,clpa] = size(inp.T{1});

start=zeros(1,rpa);
start_state = 290;
start(1,start_state)=0.20; 
start(1,total_states + start_state) = 0.80; 

% Heading
inp.discRate = 0.97;
inp.act = 4;
inp.states = rpa; % 887041

[rpa,n_obs] = size(inp.O{1});
inp.obs = n_obs;
inp.b0 = start; 

fprintf(head01,'\n');
fprintf(head01,'%s %f\n','discount:',inp.discRate);
fprintf(head01,'%s %1s\n','values:','reward');
fprintf(head01,'%s %d\n','states:',rpa);
fprintf(head01,'%s %d\n','actions:',inp.act);
fprintf(head01,'%s %d\n','observations:',inp.obs);

%%% Initial state definition
n_lines = 5;
fprintf(head01,'\n');
% fprintf(head01,'%s %d\n','start:', start_state-1);
fprintf(head01,'%s','start:');
% fprintf(head01,'\n');
for i=1:rpa
%    if mod(i, 1e4) == 1 && i > 1
%        fprintf(head01, '\n');
%    end
   fprintf(head01,' %2.1f',inp.b0(i)); 
end

% fprintf(head01,'%s %s %d %d\n','start', 'include:', start_state-1, total_states + start_state-1);

fclose(head01); % Closing the file
toc
%% % Transition matrices
tic
parfor j = 1:length(inp.T)
    fid = fopen( sprintf( 'T_245_585_10y%i.txt',j),'w');
    [Tstart, Tend, T] = find(inp.T{j});
    fprintf(fid,'\n');
    for i = 1:length(T)
        fprintf(fid,'%s %-6d %s %-6d %s %-6d %11.10f\n',...
            'T:',j-1,':',Tstart(i)-1,':',Tend(i)-1,T(i));%start to end
    end
    fclose(fid);
end
system('copy /b T_245_585_10y1.txt+T_245_585_10y2.txt+T_245_585_10y3.txt+T_245_585_10y4.txt T_245_585_10y.txt')
t = toc
    
%% % Emission matrices
tic
parfor j = 1:length(inp.O)
    fid = fopen( sprintf( 'O_245_585_10y%i.txt',j),'w');
    [Ostart, Oend, O] = find(inp.O{j});
    fprintf(fid,'\n');
    for i = 1:length(O)
        fprintf(fid,'%s %-6d %s %-6d %s %-6d %11.10f\n',...
            'O:',j-1,':',Ostart(i)-1,':', Oend(i)-1, O(i));
    end
    fclose(fid);
end
system('copy /b O_245_585_10y1.txt+O_245_585_10y2.txt+O_245_585_10y3.txt+O_245_585_10y4.txt O_245_585_10y.txt')
o = toc
    
    
%% % Reward matrices
tic
parfor j = 1:length(inp.R)
    fid = fopen( sprintf( 'R_245_585_10y_t%i.txt',j),'w');
    [Rstart, Rend, R] = find(inp.R{j});
    fprintf(fid,'\n');
    for i = 1:length(R)
            fprintf(fid,'%s %-6d %s %-6d %s %s %s %s %-10.6f\n',...
    'R:',j-1,':',Rstart(i)-1,':','*',':','*', R(i));
    end
    fclose(fid);
end
system('copy /b R_245_585_10y_t1.txt+R_245_585_10y_t2.txt+R_245_585_10y_t3.txt+R_245_585_10y_t4.txt R_245_585_10y_t.txt')
r = toc

%% 
tic
system('copy /b heading_pomdp_245_585_10y_model2_high.txt + T_245_585_10y.txt + O_245_585_10y.txt + R_245_585_10y.txt pomdp_10_years_245_585_final_model2_high.pomdp')
toc