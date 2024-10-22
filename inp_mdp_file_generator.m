% finite horizon input file

name = 'two_dikes';
extens = '.mdp';
filename = join([name extens],'');

horizon = 40;
rate = 40;
core_states = 5544;
states = core_states*rate;
n_systems = 4;
total_states = (states)*4 + 1;% 3 system and one terminal state extra
start=zeros(1,total_states);
start_state = 290;%slr state = 5 and storm state 2, therefore combined state = 290 = 72*4 + 2 for system A
start(1,start_state)=1; 
%

%%% Prescription of utilities
inp.discRate = 0.97;
inp.act = 4;
inp.states = total_states; % 887041
inp.b0 = start; 
inp.s0 = start_state;
tic

head01 = fopen('heading.txt','w'); %opening the file

[rpa,clpa] = size(inp.T{1});
%%% Heading
fprintf(head01,'\n');
fprintf(head01,'%s %f\n','discount:',inp.discRate);
fprintf(head01,'%s %1s\n','values:','reward');
fprintf(head01,'%s %d\n','states:',rpa);
fprintf(head01,'%s %d\n','actions:',inp.act);

%%% Initial state definition
start_state = 290;
fprintf(head01,'\n');
fprintf(head01,'%s %d\n','start:', start_state-1);
fclose(head01); % Closing the file
toc

tic
%% % Transition matrices
parfor j = 1:length(inp.T)
    fid = fopen( sprintf( 'T_%i.txt',j),'w');
    [Tstart, Tend, T] = find(inp.T{j});
    fprintf(fid,'\n');
    for i = 1:length(T)
        fprintf(fid,'%s %-6d %s %-6d %s %-6d %11.10f\n',...
            'T:',j-1,':',Tstart(i)-1,':',Tend(i)-1, T(i));%start to end
    end
    fclose(fid);
%     toc
end
system('copy /b T_1.txt+T_2.txt+T_3.txt+T_4.txt T.txt')
t = toc

tic    
%%
% Reward matrices
parfor j = 1:length(inp.R_fp)
    fid = fopen( sprintf( 'R_%i.txt',j),'w');
    [Rstart, Rend, R] = find(inp.R_fp{j});
    fprintf(fid,'\n');
    for i = 1:length(R)
            fprintf(fid,'%s %-6d %s %-6d %s %s %-10.1f\n',...
    'R:',j-1,':',Rstart(i)-1,':','*',R(i));
    end
    fclose(fid);
end

system('copy /b R_1.txt+R_2.txt+R_3.txt+R_4.txt R.txt')
r = toc

% copy into final combined file
%% 
tic
system('copy /b heading.txt+T.txt+R.txt mdp_file.mdp')
toc