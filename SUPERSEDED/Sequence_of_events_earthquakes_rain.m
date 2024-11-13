close all
clear
clc

%% Define hazard curves

% Resolution
resol = 365;

% Independent hazards
% Hazard curve earthquake
% Zone 923 from Iervolino et al.(2018)
% M_vec_main = 4.15:0.3:7.45;
% lambda_vec_main = (1/resol)*[0.6448 0.2326 0.1334 0.0567 0.0340 0.0255 0.0149 0.0128 0.0071 0.0028 0.0014 0];
M_vec_main = 4.45:0.3:7.45;
lambda_vec_main = (1/resol)*[0.2326 0.1334 0.0567 0.0340 0.0255 0.0149 0.0128 0.0071 0.0028 0.0014 0];

% Successive hazards (Type B)
% Hazard curve aftershock
% lambda_vec_after = (1/resol)*[0.6448 0.2326 0.1334 0.0567 0.0340 0.0255 0.0149 0.0128 0.0071 0.0028 0.0014 0];
lambda_vec_after = (1/resol)*[0.2326 0.1334 0.0567 0.0340 0.0255 0.0149 0.0128 0.0071 0.0028 0.0014 0];

% Parameters of Omori law
a = -1.66;
b = 0.96;
c = 0.03;
p = 0.93;

% Rain event
% IDF from Excel file - Gauge N05 from Tang and Cheung (2011)
Data = xlsread("IDF Curve Data.xlsx");
D_data = Data(3:end,1); % Hours
R_data = Data(3:end,2:end);
T_data = Data(1,2:end);
lambda_data = 1./T_data;
D_data_mat = D_data.*ones(length(D_data),length(T_data));
I_data = R_data./D_data_mat;
lambda_data_mat = lambda_data.*ones(size(D_data_mat));
max_lambda = max(max(lambda_data_mat));
I_vec = linspace(min(min(I_data)),max(max(I_data)),1000);
T_mat = zeros(length(D_data),length(I_vec));
for i = 1:length(D_data)
    for j = 1:length(I_vec)
        I_data_temp = I_data(i,:);
        index = find(I_data_temp >= I_vec(j), 1, 'first');
        if isempty(index) == 0
            if index == 1
                T_mat(i,j) = T_data(1);
            elseif index == length(T_data)
                T_mat(i,j) = T_data(end);
            else
                T_mat(i,j) = (T_data(index + 1) - T_data(index))/2;
            end
        else
            T_mat(i,j) = Inf;
        end
    end
end
lambda_mat_rain = (1/resol)*(1./T_mat);
% Find empirical marginal CDF for Intensity
Vec_I = [I_data(:,9),ones(19,2).*I_data(:,8),ones(19,5).*I_data(:,7),ones(19,10).*I_data(:,6),ones(19,20).*I_data(:,5),ones(19,50).*I_data(:,4),ones(19,100).*I_data(:,3),ones(19,200).*I_data(:,2),ones(19,500).*I_data(:,1)];
Vec_I = reshape(Vec_I,[],1);
[I_cdf,I_vec_cdf] = ecdf(Vec_I);

    
%% Simulate hazards throughout life-cycle
lc = 50*resol; % Days
    
t = 0; % Days
i = 0;
t_vec = [];
typ_vec = [];
M_vec = [];
% Variables stored in vectors [main after normal_after conc]
lambda_min_main = lambda_vec_main(1);
lambda_min_rain = max(max(lambda_mat_rain));
lambda_min_after = 0;
lambda_null = 0;
% Keep track of last main shock
t_main = 0;
m_main = 0;
% Flag that there have been no landslides
flag_LS = 0;

rain_ind = 0;
rain_m = [];

while t < lc
    t;
    % [main after normal_after fire normal_fire]
    lambda_case = [lambda_min_main,lambda_min_after,lambda_null,lambda_min_rain];
    % Simulate occurrence of first hazard/switch
    lambda_comp = sum(lambda_case);
    t_occ = exprnd(1/lambda_comp);
    t = t + t_occ;
    % Compute adjusted rate of aftershock and null event
    if lambda_min_after ~= 0
        lambda_min_after_adj = (10^(a+b*(m_main - M_vec_main(1)))-10^a)/(((t-t_main) + c)^p);
        lambda_null = lambda_min_after - lambda_min_after_adj;
        lambda_case = [lambda_min_main,lambda_min_after_adj,lambda_null,lambda_min_rain];
    end
    % Check which hazard it was
    P_hz = lambda_case/sum(lambda_case);
    P_cum = cumsum(P_hz);
    cas = randsrc(1,1,[1 2 3 4; P_hz]);
    switch cas
        case 1 % Mainshock
            U_M = rand;
            lambda_m = lambda_min_main*(1-U_M);
            Meq = interp1(lambda_vec_main,M_vec_main,lambda_m);
            M_vec = [M_vec Meq];
            typ_vec = [typ_vec 1];
            t_vec = [t_vec, t];
            t_main = t;
            m_main = Meq;
            % Change rate of aftershocks
            lambda_vec_after_case = lambda_vec_after.*(M_vec_main < Meq); % Truncating the rate curve so that aftershocks are smaller than mainshocks
            % Cut vectors to first zero value, and scale vector of rates
            zeroIndex = find(lambda_vec_after_case == 0, 1);
            lambda_vec_after_case = lambda_vec_after_case(1:zeroIndex);
            M_vec_after = M_vec_main(1:zeroIndex);
            lambda_min_after = (10^(a+b*(m_main - M_vec_main(1)))-10^a)/(c^p);
            lambda_vec_after_case = lambda_vec_after_case.*(lambda_min_after/lambda_vec_after_case(1));
            lambda_null = 0;
            % Generate landslides
            % Find PGA (from Campbell 1981)
            R = 20; % Km, fixed distance
            PGA = 0.0159*exp(0.868*Meq)*(R + 0.0606*exp(0.70*Meq))^(-1.09);
            % Probability of landslide (from Parker et al. (2015)
            SL = 35; % Fixed slope, degrees
            NDS = 0; % Normalized stream to ridge distance, fixed
            G = 1; % Lithology category
            %P_LS = 1/(1 - exp(-(-7.0207 + 10.9946*PGA + 0.099*SL + 1.114*NDS)))
            P_LS = 1/(1 + exp(-(-7.0207 + 10.9946*PGA + 0.099*SL + 1.114*NDS)));
            if flag_LS == 1
                P_LS = P_LS*15*exp(-0.12*(t - t_land)/resol);
            end
            rn = rand;
            if rn < P_LS
                t = t + 0.01;
                t_vec = [t_vec, t];
                typ_vec = [typ_vec 5];
                M_vec = [M_vec 0.1];
                t_land = t;
                flag_LS = 1;
            end
        case 2 % Aftershock
            U_M = rand;
            lambda_m = lambda_min_after*(1-U_M);
            Maf = interp1(lambda_vec_after_case,M_vec_after,lambda_m);
            M_vec = [M_vec Maf];
            typ_vec = [typ_vec 2];
            t_vec = [t_vec, t];
            % Recompute rate of aftershocks
            lambda_min_after = (10^(a+b*(m_main - M_vec_main(1)))-10^a)/(((t-t_main)+c)^p);
            lambda_vec_after_case = lambda_vec_after_case.*(lambda_min_after/lambda_vec_after_case(1));
            % If new rate is very small, set it to zero
            if lambda_min_after < 0.0001
                lambda_min_after = 0;
            end
            lambda_null = 0;
            % Generate landslides
            % Find PGA (from Campbell 1981)
            R = 20; % Km, fixed distance
            PGA = 0.0159*exp(0.868*Maf)*(R + 0.0606*exp(0.70*Maf))^(-1.09);
            % Probability of landslide (from Parker et al. (2015)
            SL = 35; % Fixed slope, degrees
            NDS = 0.5; % Normalized stream to ridge distance, fixed
            G = 1; % Lithology category
            P_LS = 1/(1 + exp(-(-7.0207 + 10.9946*PGA + 0.099*SL + 1.114*NDS)));
            if flag_LS == 1
                P_LS = P_LS*15*exp(-0.12*(t - t_land));
            end
            rn = rand;
            if rn < P_LS
                t = t + 0.01;
                t_vec = [t_vec, t];
                typ_vec = [typ_vec 5];
                M_vec = [M_vec 0.1];
                t_land = t;
                flag_LS = 1;
            end
        case 3 % Null event
            % Recompute rate of aftershocks
            lambda_min_after = (10^(a+b*(m_main - M_vec_main(1)))-10^a)/(((t-t_main)+c)^p);
            % If new rate is very small, set it to zero
            if lambda_min_after < 0.0001
                lambda_min_after = 0;
            end
            lambda_null = 0;
        case 4 % Rain
            U_MI = rand;
            MI = interp1(I_cdf,I_vec_cdf,U_MI);
            lambda_d_temp = interp2(I_vec,D_data,lambda_mat_rain,MI,D_data);
            % Make vector unique
            index1 = find(lambda_d_temp == 0, 1, 'last');
            index2 = find(lambda_d_temp == lambda_min_rain, 1, 'first');
            if isempty(index1) == 1
                lambda_d_temp2 = lambda_d_temp(1:index2);
                D_data_temp = D_data(1:index2);
            elseif isempty(index2) == 1
                lambda_d_temp2 = lambda_d_temp(index1:end);
                D_data_temp = D_data(index1:end);
            else
                lambda_d_temp2 = lambda_d_temp(index1:index2);
                D_data_temp = D_data(index1:index2);
            end
            [lambda_d_temp2,ia,ic] = unique(lambda_d_temp2);
            D_data_temp = D_data_temp(ia);
            U_MD = rand;
            lambda_d = max(lambda_d_temp2)*(1-U_MD);
            if length(D_data_temp) == 1
                MD = 2;
            else
                MD = interp1(lambda_d_temp2,D_data_temp,lambda_d);
            end
            if isnan(MD)
                MD = 2;
            end
            Mr = MI*MD;
            rain_ind = rain_ind + 1;
            rain_m(rain_ind,:) = [MI,MD];
            M_vec = [M_vec Mr];
            typ_vec = [typ_vec 4];
            t_vec = [t_vec, t];
            % Generate landslides
            % Critical duration for fixed intensity (From Liu and Wang
            % 2022)
            D_c = 9.00*1000*MI^(-2.149) + 1.61;
            % D_c = 2.61*10000*MI^(-2.030) + 27.04; % after stabilization
            if MD > D_c
                t = t + 0.01;
                t_vec = [t_vec, t];
                typ_vec = [typ_vec 5];
                M_vec = [M_vec 0.1];
                t_land = t;
                flag_LS = 1;
            end
    end
end

figure(1)
% Mainshocks
t_ms_vec = t_vec(typ_vec == 1);
m_ms_vec = M_vec(typ_vec == 1);
% Aftershocks
t_as_vec = t_vec(typ_vec == 2);
m_as_vec = M_vec(typ_vec == 2);
% Rain
t_r_vec = t_vec(typ_vec == 4);
m_r_vec = M_vec(typ_vec == 4);
% Landslide
t_l_vec = t_vec(typ_vec == 5);
m_l_vec = M_vec(typ_vec == 5);

scatter(t_as_vec/resol,ones(1,length(t_as_vec)),exp(m_as_vec)/2,'magenta')
hold on
scatter(t_ms_vec/resol,ones(1,length(t_ms_vec)),exp(m_ms_vec),'b','LineWidth',1.5)
scatter(t_r_vec/resol,ones(1,length(t_r_vec)),(rain_m(:,1)+1),'k')
scatter(t_ms_vec/resol,2*ones(1,length(t_ms_vec)),exp(m_ms_vec),'b','LineWidth',1.5)
scatter(t_as_vec/resol,3*ones(1,length(t_as_vec)),exp(m_as_vec)/2,'magenta')
scatter(t_r_vec/resol,ones(1,length(t_r_vec)),(rain_m(:,2)+1),'k')
scatter(t_r_vec/resol,4*ones(1,length(t_r_vec)),(rain_m(:,1)+1),'k')
scatter(t_r_vec/resol,5*ones(1,length(t_r_vec)),(rain_m(:,2)+1),'k')
scatter(t_l_vec/resol,ones(1,length(t_l_vec)),'+','r','Linewidth',1.5)
scatter(t_l_vec/resol,6*ones(1,length(t_l_vec)),'+','r','Linewidth',1.5)
hold off
xlim([0 50])
yticks([1 2 3 4 5 6])
xlabel('Time [yr]','Interpreter','latex','FontSize',16)
yticklabels({'All','Main Shock','Aftershock','Intensity','Duration','Landslide'})
ax = gca;
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% t_vec contains the times of occurrence of the shocks
% M_vec contains the magnitude of the events
% typ_vec is 1 for main shocks and 2 for aftershocks

% Number of mainshocks and aftershocks and maximum magnitude
typ_vec(t_vec>lc) = [];
M_vec(t_vec>lc)   = [];   
t_vec(t_vec>lc)   = [];   

N_main_shocks = sum(typ_vec == 1);
N_after_shocks = sum(typ_vec == 2);
N_rain_events = sum(typ_vec == 4);
N_landslides = sum(typ_vec == 5);
Max_main = max(M_vec);

save ('scenario_er.mat','t_vec','M_vec','typ_vec');