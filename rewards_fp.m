% generate rewards for fast parser

rewards = inp.R;
trans = inp.T;
inp.act = 4;

rewards_ = cell(1,inp.act);

parfor i = 1:inp.act
    [rpa,clpa] = size(rewards{i});
     R_matrix = rewards{i};
     T_matrix = trans{i};

     R = zeros(rpa,1);
     for j = 1:rpa

        R(j) = R_matrix(j,:)*T_matrix(j,:)';
     end

    rewards_{i} = R;
end
inp.R_fp = rewards_;

        