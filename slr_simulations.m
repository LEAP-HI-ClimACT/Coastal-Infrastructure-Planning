% sample realizations

x = years1;

figure(10)
% historical observations
years = [1900:2021];
p = plot(years(51:122), annual_msl_rel_ref_ipcc(51:122)*100,'linewidth', 2);
p.Color = [0.7 0.7 0.7];
xlim([1950, 2150])
hold on

% projections
% ssp 245
ylow1 = SLR_245(:,8)'/10;
yhigh1 = SLR_245(:,100)'/10;
ymid1 = SLR_245(:, 54)'/10;
xconf = [x x(end:-1:1)] ;         
yconf1 = [yhigh1 ylow1(end:-1:1)];

% figure
p = fill(xconf,yconf1,'blue');
p.FaceColor = [0.7 0.7 1];      
p.EdgeColor = 'none';
p.FaceAlpha = 0.6;

hold on

legend
title('mean sea level projections to 2150 relative to 1995-2014 bsaeline');
xlabel('Year')
ylabel('msl relative to 1995-2014 baseline (cm)')
xticks([1950 1960 1970 1980 1990 2000 2010 2020 2030 2040 2050 2060 2070 2080 2090 2100 2110 2120 2130 2140 2150])
hold on

% sample realizations
color1 = [200, 200, 200] / 255;
n_samples = 5;% 100 million samples
years_sample = [2020:2150];
n_years = length(years_sample);
SLR_year = zeros(n_years, quantiles);

parfor i = 1: quantiles
    for j = 1:n_years
        SLR_year(j,i)= interp1(years1, SLR_245(:,i)'/10, years_sample(j));
    end
end

% sample_q = rand(n_samples,1);% try within 5th and 95th percentile
sample_q = [0.1, 0.4, 0.5, 0.6, 0.9];
values=zeros(n_samples, n_years);

for i = 1:n_samples
    u = sample_q(i);
    for j = 1:n_years
        values(i,j)= interp1(quantiles1, SLR_year(j,:), u)+ 100*0.0318*randn(1,1);
    end
%    figure(4)
   plot(years_sample, values(i,:), 'color', [0.5 0.5 1], 'linewidth', 1)
%    hold on
end