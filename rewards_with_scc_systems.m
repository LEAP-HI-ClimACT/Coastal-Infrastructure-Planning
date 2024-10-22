% correct way to calculate the rewards with the scc accumulation from the
% current time it is emitted

load('scc.mat');% use any scc values
horizon1 = 40;
extra_years = 20; % for terminal costs calculation
horizon = horizon1 + extra_years;
scc(41:60) = scc(40);
discounted_sum_scc = zeros(1, horizon);
discount_factor = 0.97;
for i = 1:horizon
    scc_i = scc(i:end);
    discount = 1;
    for j = 1:length(scc_i)
        discounted_sum_scc(i) = discounted_sum_scc(i) + discount*scc_i(j);
        discount = discount*discount_factor;           
    end
end

n_systems = 4;
s = 0.0085;% slope of city (without seawall)
gamma = 0.97;
Dike_Height = 1.5;%m
cc_dike = -Dike_Construction_Cost1(Dike_Height); %$/m length along shore
cm_dike = -100;%$/m length along shore
f_damage = 0.07;% 20 billions
val_z = 15.3*1e6;
cf = -f_damage*val_z;
vol_z = (8.50*(1/s)*8.50)*0.5;

[cc_carbon_dike, ] = - DikeConstructionMaterialGHGCost1(Dike_Height,1);
[cm_carbon_dike, ] = - AnnualMaintenanceGHGCost(1);

% do nothing
% system A
R_dn = cell(1,n_systems); 
slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards1_dn = zeros(1, n_states);

for i = 1:n_states
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;%m
        area1 = 0.5*(1/s)*totalh1^2;

        vol_f(i) = (area1 - 0);%area in m^2
        c_flood = cf*vol_f(i)/vol_z;
        R_flood1(i) = c_flood; 
end

R_flood_carbon1 = zeros(horizon, n_states);
for t = 1:horizon
    c_carbon_flood1 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood1(i)] = FloodDamageGHGCost(R_flood1(i), discounted_sum_scc(t));      
end
R_flood_carbon1(t,:) = c_carbon_flood1; 
end

%total
rewards1_dn = cell(1, horizon);
for t = 1:horizon
    rewards1_dn{1,t} = (R_flood1 + R_flood_carbon1(t,:))*gamma;
end
R_dn{1,1} = rewards1_dn;

% system B
slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards2_dn = zeros(1, n_states);
b = 0.20;%m
t = 1.70;%m (dike height = 1.5 m)

for i = 1:n_states
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;
        if (totalh1> b) && (totalh1 <= t)
            area1 = 0.5*b*b*(1/s) + (totalh1 - b)*b*(1/s);
        else
            area1 = 0.5*(1/s)*totalh1^2;
        end
       
        vol_f(i) = (area1 - 0);
        c_flood = cf*vol_f(i)/vol_z;
        R_flood2(i) = c_flood;
end

R_flood_carbon2 = zeros(40, n_states);
for t = 1:horizon
    c_carbon_flood2 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood2(i)] = FloodDamageGHGCost(R_flood2(i), discounted_sum_scc(t));      
end
R_flood_carbon2(t,:) = c_carbon_flood2; 
end

R_m = cm_dike * ones(1, n_states);
R_mc = zeros(40, n_states);
for t = 1:horizon
    Rmc(t) = cm_carbon_dike*discounted_sum_scc(t);
    R_mc(t,:) = ones(1, n_states)*Rmc(t);
end

%total
rewards2_dn = cell(1, horizon);
for t = 1:horizon
    rewards2_dn{1,t} = (R_flood2 + R_flood_carbon2(t,:))*gamma + R_m + R_mc(t,:);
end
R_dn{1,2} = rewards2_dn;

% system C
slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards3_dn = zeros(1, n_states);
b = 1.90;%m
t = 3.40;%m

for i = 1:n_states
    
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;
        if (totalh1> b) && (totalh1 <= t)
            area1 = 0.5*b*b*(1/s) + (totalh1 - b)*b*(1/s);
        else
            area1 = 0.5*(1/s)*totalh1^2;
        end
        
        vol_f(i) = (area1 - 0);%area in m^2
        c_flood = cf*vol_f(i)/vol_z;
        R_flood3(i) = c_flood;

end
%risk carbon
R_flood_carbon3 = zeros(40, n_states);
for t = 1:horizon
    c_carbon_flood3 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood3(i)] = FloodDamageGHGCost(R_flood3(i), discounted_sum_scc(t));      
end
R_flood_carbon3(t,:) = c_carbon_flood3; 
end

%maintenance
R_mc = zeros(40, n_states);
R_m = cm_dike * ones(1, n_states);
for t = 1:horizon
    Rmc(t) = cm_carbon_dike*discounted_sum_scc(t);
    R_mc(t,:) = ones(1, n_states)*Rmc(t);
end

%total
rewards3_dn = cell(1, 40);
for t = 1:horizon
    rewards3_dn{1,t} = (R_flood3 + R_flood_carbon3(t,:))*gamma + R_m + R_mc(t,:);
end
R_dn{1,3} = rewards3_dn;

% system D
slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards4_dn = zeros(1, n_states);
b1 = 0.20;%m
t1 = 1.70;%m
b2 = 1.90;%m
t2 = 3.40;%m
for i = 1:n_states
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;
        if totalh1 <= b1
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b1 && totalh1 <=t1
            area1 = 0.5*(1/s)*b1^2 + (totalh1 - b1)*b1*(1/s);
        elseif totalh1>t1 && totalh1<=b2
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b2 && totalh1<=t2
            area1 = 0.5*(1/s)*b2^2 + (totalh1 - b2)*b2*(1/s);
        elseif totalh1>t2
            area1 = 0.5*(1/s)*totalh1^2;
        end
        
        vol_f(i) = (area1 - 0);%area in m^2
        c_flood = cf*vol_f(i)/vol_z;
        R_flood4(i) = c_flood;
        
end

%risk carbon
R_flood_carbon4 = zeros(40, n_states);
for t = 1:horizon
    c_carbon_flood4 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood4(i)] = FloodDamageGHGCost(R_flood4(i), discounted_sum_scc(t));      
end
R_flood_carbon4(t,:) = c_carbon_flood4; 
end

%maintenance
R_mc = zeros(40, n_states);
R_m = (2*cm_dike) * ones(1, n_states);
for t = 1:horizon
    Rmc(t) = (cm_carbon_dike)*discounted_sum_scc(t) + cm_carbon_dike*discounted_sum_scc(t);
    R_mc(t,:) = ones(1, n_states)*Rmc(t);
end

%total
rewards4_dn = cell(1, horizon);
for t = 1:horizon
    rewards4_dn{1,t} = (R_flood4 + R_flood_carbon4(t,:))*gamma + R_m + R_mc(t, :);
end
R_dn{1,4} = rewards4_dn;

% action 2- construct dike 1 only
R_1 = cell(1, n_systems);

% system A
slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards1_dike1 = zeros(1, n_states);
b = 0.20;%m
t = 1.70;%m (dike height = 1.5 m)

for i = 1:n_states
  
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;
        if (totalh1> b) && (totalh1 <= t)
            area1 = 0.5*b*b*(1/s) + (totalh1 - b)*b*(1/s);
        else
            area1 = 0.5*(1/s)*totalh1^2;
        end
        
        vol_f(i) = (area1 - 0);%area in m^2
        c_flood = cf*vol_f(i)/vol_z;
        R_flood1(i) = c_flood;

end

%risk carbon
R_flood_carbon1 = zeros(40, n_states);
for t = 1:horizon
    c_carbon_flood1 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood1(i)] = FloodDamageGHGCost(R_flood1(i), discounted_sum_scc(t));      
end
R_flood_carbon1(t,:) = c_carbon_flood1; 
end

% current construction (carbon gets emitted in the atmosphere)
R_c = cc_dike*ones(1, n_states);
R_const = zeros(40, n_states);
for t = 1:horizon
    Rconst(t) = (cc_carbon_dike)*discounted_sum_scc(t);
    R_const(t,:) = ones(1, n_states)*Rconst(t);
end

%total
rewards1_dike1 = cell(1, 40);
for t = 1:horizon
    rewards1_dike1{1,t} = (R_flood1 + R_flood_carbon1(t,:))*gamma + R_const(t, :) + R_c;
end

R_1{1,1} = rewards1_dike1;

% system B
slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards2_dike1 = zeros(1, n_states);
b = 0.20;%m
t = 1.70;%m (dike height = 1.5 m)

for i = 1:n_states
   
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;
        if (totalh1> b) && (totalh1 <= t)
            area1 = 0.5*b*b*(1/s) + (totalh1 - b)*b*(1/s);
        else
            area1 = 0.5*(1/s)*totalh1^2;
        end
        
        vol_f(i) = (area1 - 0);%area in m^2
        c_flood = cf*vol_f(i)/vol_z; 
        R_flood2(i) = c_flood;
     
end
%risk carbon
R_flood_carbon2 = zeros(40, n_states);
for t = 1:horizon
    c_carbon_flood2 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood2(i)] = FloodDamageGHGCost(R_flood2(i), discounted_sum_scc(t));      
end
R_flood_carbon2(t,:) = c_carbon_flood2; 
end

% current construction (carbon gets emitted in the atmosphere)
R_c = cc_dike*ones(1, n_states);
R_const = zeros(40, n_states);
for t = 1:horizon
    Rconst(t) = (cc_carbon_dike)*discounted_sum_scc(t);
    R_const(t,:) = ones(1, n_states)*Rconst(t);
end

%maintenance
R_m = cm_dike * ones(1, n_states);
R_mc = zeros(40, n_states);
for t = 1:horizon
    Rmc(t) = cm_carbon_dike*discounted_sum_scc(t);
    R_mc(t,:) = ones(1, n_states)*Rmc(t);
end

%total
rewards2_dike1 = cell(1, 40);
for t = 1:horizon
    rewards2_dike1{1,t} = (R_flood2 + R_flood_carbon2(t,:))*gamma + R_const(t, :) + R_c + R_m + R_mc(t,:);
end

R_1{1,2} = rewards2_dike1;

% system C
slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards3_dike1 = zeros(1, n_states);
b1 = 0.20;%m
t1 = 1.70;%m
b2 = 1.90;
t2 = 3.40;

for i = 1:n_states
    
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;
        if totalh1 <= b1
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b1 && totalh1 <=t1
            area1 = 0.5*(1/s)*b1^2 + (totalh1 - b1)*b1*(1/s);
        elseif totalh1>t1 && totalh1<=b2
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b2 && totalh1<=t2
            area1 = 0.5*(1/s)*b2^2 + (totalh1 - b2)*b2*(1/s);
        elseif totalh1>t2
            area1 = 0.5*(1/s)*totalh1^2;
        end
        
        vol_f(i) = (area1 - 0);%area in m^2
        c_flood = cf*vol_f(i)/vol_z;
        R_flood3(i) = c_flood;
end
%risk carbon
R_flood_carbon3 = zeros(40, n_states);
for t = 1:horizon
    c_carbon_flood3 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood3(i)] = FloodDamageGHGCost(R_flood3(i), discounted_sum_scc(t));      
end
R_flood_carbon3(t,:) = c_carbon_flood3; 
end

% current construction (carbon gets emitted in the atmosphere)
R_c = cc_dike*ones(1, n_states);
R_const = zeros(40, n_states);
for t = 1:horizon
    Rconst(t) = (cc_carbon_dike)*discounted_sum_scc(t);
    R_const(t,:) = ones(1, n_states)*Rconst(t);
end

%maintenance
R_m = cm_dike * ones(1, n_states);
R_mc = zeros(40, n_states);
for t = 1:horizon
    Rmc(t) = cm_carbon_dike*discounted_sum_scc(t);
    R_mc(t,:) = ones(1, n_states)*Rmc(t);
end

%total
rewards3_dike1 = cell(1, 40);
for t = 1:horizon
    rewards3_dike1{1,t} = (R_flood3 + R_flood_carbon3(t,:))*gamma + R_const(t, :) + R_c + R_m + R_mc(t,:);
end
R_1{1,3} = rewards3_dike1;

% system D
slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards4_dike1 = zeros(1, n_states);
b1 = 0.20;%m
t1 = 1.70;%m
b2 = 1.90;
t2 = 3.40;
for i = 1:n_states
    
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;
        if totalh1 <= b1
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b1 && totalh1 <=t1
            area1 = 0.5*(1/s)*b1^2 + (totalh1 - b1)*b1*(1/s);
        elseif totalh1>t1 && totalh1<=b2
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b2 && totalh1<=t2
            area1 = 0.5*(1/s)*b2^2 + (totalh1 - b2)*b2*(1/s);
        elseif totalh1>t2
            area1 = 0.5*(1/s)*totalh1^2;
        end
        
        vol_f(i) = (area1 - 0);%area in m^2
        c_flood = cf*vol_f(i)/vol_z;
        R_flood4(i) = c_flood;

end
%risk carbon
R_flood_carbon4 = zeros(40, n_states);
for t = 1:horizon
    c_carbon_flood4 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood4(i)] = FloodDamageGHGCost(R_flood4(i), discounted_sum_scc(t));      
end
R_flood_carbon4(t,:) = c_carbon_flood4; 
end

% current construction (carbon gets emitted in the atmosphere)
R_c = cc_dike*ones(1, n_states);
R_const = zeros(40, n_states);
for t = 1:horizon
    Rconst(t) = (cc_carbon_dike)*discounted_sum_scc(t);
    R_const(t,:) = ones(1, n_states)*Rconst(t);
end
%maintenance
R_m = (cm_dike + cm_dike) * ones(1, n_states);
R_mc = zeros(40, n_states);
for t = 1:horizon
    Rmc(t) = (cm_carbon_dike)*discounted_sum_scc(t) + cm_carbon_dike*discounted_sum_scc(t);
    R_mc(t,:) = ones(1, n_states)*Rmc(t);
end

%total
rewards4_dike1 = cell(1, 40);
for t = 1:horizon
    rewards4_dike1{1,t} = (R_flood4 + R_flood_carbon4(t,:))*gamma + R_const(t, :) + R_c + R_m + R_mc(t,:);
end
R_1{1,4} = rewards4_dike1;

% action 3- construct dike 2 only
R_2 = cell(1,n_systems);

% system A
slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards1_dike2 = zeros(1, n_states);
b = 1.90;%cm
t = 3.40;%cm
for i = 1:n_states
    
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;
        if (totalh1> b) && (totalh1 <= t)
            area1 = 0.5*b*b*(1/s) + (totalh1 - b)*b*(1/s);
        else
            area1 = 0.5*(1/s)*totalh1^2;
        end
        
        vol_f(i) = (area1 - 0);%area in m^2
        c_flood = cf*vol_f(i)/vol_z;
        R_flood1(i) = c_flood;

end
%risk carbon
R_flood_carbon1 = zeros(40, n_states);
for t = 1:horizon
    c_carbon_flood1 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood1(i)] = FloodDamageGHGCost(R_flood1(i), discounted_sum_scc(t));      
end
R_flood_carbon1(t,:) = c_carbon_flood1; 
end

% current construction (carbon gets emitted in the atmosphere)
R_c = cc_dike*ones(1, n_states);
R_const = zeros(40, n_states);
for t = 1:horizon
    Rconst(t) = (cc_carbon_dike)*discounted_sum_scc(t);
    R_const(t,:) = ones(1, n_states)*Rconst(t);
end

%total
rewards1_dike2 = cell(1, 40);
for t = 1:horizon
    rewards1_dike2{1,t} = (R_flood1 + R_flood_carbon1(t,:))*gamma + R_const(t, :) + R_c;
end

R_2{1,1} = rewards1_dike2;

% system B
slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards2_dike2 = zeros(1, n_states);
b1 = 0.20;%m
t1 = 1.70;%m
b2 = 1.90;
t2 = 3.40;
for i = 1:n_states
    
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;
        if totalh1 <= b1
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b1 && totalh1 <=t1
            area1 = 0.5*(1/s)*b1^2 + (totalh1 - b1)*b1*(1/s);
        elseif totalh1>t1 && totalh1<=b2
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b2 && totalh1<=t2
            area1 = 0.5*(1/s)*b2^2 + (totalh1 - b2)*b2*(1/s);
        elseif totalh1>t2
            area1 = 0.5*(1/s)*totalh1^2;
        end
        
        vol_f(i) = (area1 - 0);%area in m^2
        c_flood = cf*vol_f(i)/vol_z;
        R_flood2(i) = c_flood;

end
%risk carbon
R_flood_carbon2 = zeros(40, n_states);
for t = 1:horizon
    c_carbon_flood2 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood2(i)] = FloodDamageGHGCost(R_flood2(i), discounted_sum_scc(t));      
end
R_flood_carbon2(t,:) = c_carbon_flood2; 
end

% current construction (carbon gets emitted in the atmosphere)
R_c = cc_dike*ones(1, n_states);
R_const = zeros(40, n_states);
for t = 1:horizon
    Rconst(t) = (cc_carbon_dike)*discounted_sum_scc(t);
    R_const(t,:) = ones(1, n_states)*Rconst(t);
end

%maintenance
R_m = cm_dike * ones(1, n_states);
R_mc = zeros(40, n_states);
for t = 1:horizon
    Rmc(t) = cm_carbon_dike*discounted_sum_scc(t);
    R_mc(t,:) = ones(1, n_states)*Rmc(t);
end

%total
rewards2_dike2 = cell(1, 40);
for t = 1:horizon
    rewards2_dike2{1,t} = (R_flood2 + R_flood_carbon2(t,:))*gamma + R_const(t, :) + R_c + R_m + R_mc(t,:);
end

R_2{1,2} = rewards2_dike2;

% system C
slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards3_dike2 = zeros(1, n_states);
b = 1.90;%cm
t = 3.40;%cm
for i = 1:n_states
    
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;
        if (totalh1> b) && (totalh1 <= t)
            area1 = 0.5*b*b*(1/s) + (totalh1 - b)*b*(1/s);
        else
            area1 = 0.5*(1/s)*totalh1^2;
        end
        
        vol_f(i) = (area1 - 0);%area in m^2
        c_flood = cf*vol_f(i)/vol_z;
        R_flood3(i) = c_flood;

end
%risk carbon
R_flood_carbon3 = zeros(40, n_states);
for t = 1:horizon
    c_carbon_flood3 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood3(i)] = FloodDamageGHGCost(R_flood3(i), discounted_sum_scc(t));      
end
R_flood_carbon3(t,:) = c_carbon_flood3; 
end

% current construction (carbon gets emitted in the atmosphere)
R_c = cc_dike*ones(1, n_states);
R_const = zeros(40, n_states);
for t = 1:horizon
    Rconst(t) = (cc_carbon_dike)*discounted_sum_scc(t);
    R_const(t,:) = ones(1, n_states)*Rconst(t);
end

%maintenance
R_m = cm_dike * ones(1, n_states);
R_mc = zeros(40, n_states);
for t = 1:horizon
    Rmc(t) = cm_carbon_dike*discounted_sum_scc(t);
    R_mc(t,:) = ones(1, n_states)*Rmc(t);
end

%total
rewards3_dike2 = cell(1, 40);
for t = 1:horizon
    rewards3_dike2{1,t} = (R_flood3 + R_flood_carbon3(t,:))*gamma + R_const(t, :) + R_c + R_m + R_mc(t,:);
end

R_2{1,3} = rewards3_dike2;

% system D
slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards4_dike2 = zeros(1, n_states);
b1 = 0.20;%m
t1 = 1.70;%m
b2 = 1.90;
t2 = 3.40;
for i = 1:n_states
    
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;
        if totalh1 <= b1
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b1 && totalh1 <=t1
            area1 = 0.5*(1/s)*b1^2 + (totalh1 - b1)*b1*(1/s);
        elseif totalh1>t1 && totalh1<=b2
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b2 && totalh1<=t2
            area1 = 0.5*(1/s)*b2^2 + (totalh1 - b2)*b2*(1/s);
        elseif totalh1>t2
            area1 = 0.5*(1/s)*totalh1^2;
        end
        
        vol_f(i) = (area1 - 0);%area in m^2
        c_flood = cf*vol_f(i)/vol_z;
        R_flood4(i) = c_flood;

end
%risk carbon
R_flood_carbon4 = zeros(40, n_states);
for t = 1:horizon
    c_carbon_flood4 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood4(i)] = FloodDamageGHGCost(R_flood4(i), discounted_sum_scc(t));      
end
R_flood_carbon4(t,:) = c_carbon_flood4; 
end

% current construction (carbon gets emitted in the atmosphere)
R_c = cc_dike*ones(1, n_states);
R_const = zeros(40, n_states);
for t = 1:horizon
    Rconst(t) = (cc_carbon_dike)*discounted_sum_scc(t);
    R_const(t,:) = ones(1, n_states)*Rconst(t);
end


%maintenance
R_m = (cm_dike + cm_dike) * ones(1, n_states);
R_mc = zeros(40, n_states);
for t = 1:horizon
    Rmc(t) = (cm_carbon_dike + cm_carbon_dike)*discounted_sum_scc(t);
    R_mc(t,:) = ones(1, n_states)*Rmc(t);
end

%total
rewards4_dike2 = cell(1, 40);
for t = 1:horizon
    rewards4_dike2{1,t} = (R_flood4 + R_flood_carbon4(t,:))*gamma + R_const(t, :) + R_c + R_m + R_mc(t,:);
end
R_2{1,4} = rewards4_dike2;

% action 4: construct both
R_12 = cell(1, n_systems);

% system A

slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards1_dikes12 = zeros(1, n_states);
b1 = 0.20;%m
t1 = 1.70;%m
b2 = 1.90;
t2 = 3.40;

for i = 1:n_states
   
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;
        if totalh1 <= b1
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b1 && totalh1 <=t1
            area1 = 0.5*(1/s)*b1^2 + (totalh1 - b1)*b1*(1/s);
        elseif totalh1>t1 && totalh1<=b2
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b2 && totalh1<=t2
            area1 = 0.5*(1/s)*b2^2 + (totalh1 - b2)*b2*(1/s);
        elseif totalh1>t2
            area1 = 0.5*(1/s)*totalh1^2;
        end
        
        vol_f(i) = (area1 - 0);%area in m^2
        c_flood = cf*vol_f(i)/vol_z;
        R_flood1(i) = c_flood;

end
%risk carbon
R_flood_carbon1 = zeros(40, n_states);
for t = 1:horizon
    c_carbon_flood1 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood1(i)] = FloodDamageGHGCost(R_flood1(i), discounted_sum_scc(t));      
end
R_flood_carbon1(t,:) = c_carbon_flood1; 
end

% current construction (carbon gets emitted in the atmosphere)
R_c = (cc_dike + cc_dike)*ones(1, n_states);
R_const = zeros(40, n_states);
for t = 1:horizon
    Rconst(t) = (cc_carbon_dike + cc_carbon_dike)*discounted_sum_scc(t);
    R_const(t,:) = ones(1, n_states)*Rconst(t);
end

%total
rewards1_dikes12 = cell(1, 40);
for t = 1:horizon
    rewards1_dikes12{1,t} = (R_flood1 + R_flood_carbon1(t,:))*gamma + R_const(t, :) + R_c;
end

R_12{1,1} = rewards1_dikes12;

% system B
slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards2_dikes12 = zeros(1, n_states);
b1 = 0.20;%m
t1 = 1.70;%m
b2 = 1.90;
t2 = 3.40;
for i = 1:n_states
   
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;
        if totalh1 <= b1
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b1 && totalh1 <=t1
            area1 = 0.5*(1/s)*b1^2 + (totalh1 - b1)*b1*(1/s);
        elseif totalh1>t1 && totalh1<=b2
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b2 && totalh1<=t2
            area1 = 0.5*(1/s)*b2^2 + (totalh1 - b2)*b2*(1/s);
        elseif totalh1>t2
            area1 = 0.5*(1/s)*totalh1^2;
        end
        
        vol_f(i) = (area1 - 0);%area in m^2
        c_flood = cf*vol_f(i)/vol_z;
        R_flood2(i) = c_flood;

end
%risk carbon
R_flood_carbon2 = zeros(40, n_states);
for t = 1:horizon
    c_carbon_flood2 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood2(i)] = FloodDamageGHGCost(R_flood2(i), discounted_sum_scc(t));      
end
R_flood_carbon2(t,:) = c_carbon_flood2; 
end

% current construction (carbon gets emitted in the atmosphere)
R_c = (cc_dike+cc_dike)*ones(1, n_states);
R_const = zeros(40, n_states);
for t = 1:horizon
    Rconst(t) = (cc_carbon_dike+cc_carbon_dike)*discounted_sum_scc(t);
    R_const(t,:) = ones(1, n_states)*Rconst(t);
end

%maintenance
R_m = cm_dike * ones(1, n_states);
R_mc = zeros(40, n_states);
for t = 1:horizon
    Rmc(t) = cm_carbon_dike*discounted_sum_scc(t);
    R_mc(t,:) = ones(1, n_states)*Rmc(t);
end

%total
rewards2_dikes12 = cell(1, 40);
for t = 1:horizon
    rewards2_dikes12{1,t} = (R_flood2 + R_flood_carbon2(t,:))*gamma + R_const(t, :) + R_c + R_m + R_mc(t,:);
end

R_12{1,2} = rewards2_dikes12;

% system C
slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards3_dikes12 = zeros(1, n_states);
b1 = 0.20;%m
t1 = 1.70;%m
b2 = 1.90;
t2 = 3.40;
for i = 1:n_states
    
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;
        if totalh1 <= b1
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b1 && totalh1 <=t1
            area1 = 0.5*(1/s)*b1^2 + (totalh1 - b1)*b1*(1/s);
        elseif totalh1>t1 && totalh1<=b2
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b2 && totalh1<=t2
            area1 = 0.5*(1/s)*b2^2 + (totalh1 - b2)*b2*(1/s);
        elseif totalh1>t2
            area1 = 0.5*(1/s)*totalh1^2;
        end
        
        vol_f(i) = (area1 - 0);%area in m^2
        c_flood = cf*vol_f(i)/vol_z;
        R_flood3(i) = c_flood;

end
%risk carbon
R_flood_carbon3 = zeros(40, n_states);
for t = 1:horizon
    c_carbon_flood3 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood3(i)] = FloodDamageGHGCost(R_flood3(i), discounted_sum_scc(t));      
end
R_flood_carbon3(t,:) = c_carbon_flood3; 
end

% current construction (carbon gets emitted in the atmosphere)
R_c = (cc_dike + cc_dike)*ones(1, n_states);
R_const = zeros(40, n_states);
for t = 1:horizon
    Rconst(t) = (cc_carbon_dike+cc_carbon_dike)*discounted_sum_scc(t);
    R_const(t,:) = ones(1, n_states)*Rconst(t);
end

%maintenance
R_m = cm_dike * ones(1, n_states);
R_mc = zeros(40, n_states);
for t = 1:horizon
    Rmc(t) = cm_carbon_dike*discounted_sum_scc(t);
    R_mc(t,:) = ones(1, n_states)*Rmc(t);
end

%total
rewards3_dikes12 = cell(1, 40);
for t = 1:horizon
    rewards3_dikes12{1,t} = (R_flood3 + R_flood_carbon3(t,:))*gamma + R_const(t, :) + R_c + R_m + R_mc(t,:);
end

R_12{1,3} = rewards3_dikes12;

% system D
slr_states = 77;
surge_states = 72;
n_states = slr_states*surge_states;
rewards4_dikes12 = zeros(1, n_states);
b1 = 0.20;%m
t1 = 1.70;%m
b2 = 1.90;
t2 = 3.40;
for i = 1:n_states
    
        [slr1,surge1] = components(i);
        slrh1 = values_slr(slr1);
        surgeh1 = values_surge(surge1);
        totalh1 = (slrh1 + surgeh1)/100;
        if totalh1 <= b1
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b1 && totalh1 <=t1
            area1 = 0.5*(1/s)*b1^2 + (totalh1 - b1)*b1*(1/s);
        elseif totalh1>t1 && totalh1<=b2
            area1 = 0.5*(1/s)*totalh1^2;
        elseif totalh1>b2 && totalh1<=t2
            area1 = 0.5*(1/s)*b2^2 + (totalh1 - b2)*b2*(1/s);
        elseif totalh1>t2
            area1 = 0.5*(1/s)*totalh1^2;
        end
        
        vol_f(i) = (area1 - 0);%area in m^2
        c_flood = cf*vol_f(i)/vol_z;
        R_flood4(i) = c_flood;
end

%risk carbon
R_flood_carbon4 = zeros(40, n_states);
for t = 1:horizon
    c_carbon_flood4 = zeros(1, length(n_states));
for i = 1:n_states
        [~, c_carbon_flood4(i)] = FloodDamageGHGCost(R_flood4(i), discounted_sum_scc(t));      
end
R_flood_carbon4(t,:) = c_carbon_flood4; 
end

% current construction (carbon gets emitted in the atmosphere)
R_c = (cc_dike+cc_dike)*ones(1, n_states);
R_const = zeros(40, n_states);
for t = 1:horizon
    Rconst(t) = (cc_carbon_dike+cc_carbon_dike)*discounted_sum_scc(t);
    R_const(t,:) = ones(1, n_states)*Rconst(t);
end

%maintenance
R_m = (cm_dike + cm_dike) * ones(1, n_states);
R_mc = zeros(40, n_states);
for t = 1:horizon
    Rmc(t) = (cm_carbon_dike + cm_carbon_dike)*discounted_sum_scc(t);
    R_mc(t,:) = ones(1, n_states)*Rmc(t);
end

%total
rewards4_dikes12 = cell(1, 40);
for t = 1:horizon
    rewards4_dikes12{1,t} = (R_flood4 + R_flood_carbon4(t,:))*gamma + R_const(t, :) + R_c + R_m + R_mc(t,:);
end
R_12{1,4} = rewards4_dikes12;
