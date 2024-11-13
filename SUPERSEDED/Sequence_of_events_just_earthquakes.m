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
    
%% Simulate hazards throughout life-cycle
lc = 50*resol; % Days

% for k = 1:1000
    
t = 0; % Days
i = 0;
t_vec = [];
typ_vec = [];
M_vec = [];
% Variables stored in vectors [main after normal_after conc]
lambda_min_main = lambda_vec_main(1);
lambda_min_after = 0;
lambda_null = 0;
% Keep track of last main shock
t_main = 0;
m_main = 0;

while t < lc
    t;
    % [main after normal_after fire normal_fire]
    lambda_case = [lambda_min_main,lambda_min_after,lambda_null];
    % Simulate occurrence of first hazard/switch
    lambda_comp = sum(lambda_case);
    t_occ = exprnd(1/lambda_comp);
    t = t + t_occ;
    % Compute adjusted rate of aftershock and null event
    if lambda_min_after ~= 0
        lambda_min_after_adj = (10^(a+b*(m_main - M_vec_main(1)))-10^a)/(((t-t_main) + c)^p);
        lambda_null = lambda_min_after - lambda_min_after_adj;
        lambda_case = [lambda_min_main,lambda_min_after_adj,lambda_null];
    end
    % Check which hazard it was
    P_hz = lambda_case/sum(lambda_case);
    P_cum = cumsum(P_hz);
    cas = randsrc(1,1,[1 2 3; P_hz]);
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
        case 3 % Null event
            % Recompute rate of aftershocks
            lambda_min_after = (10^(a+b*(m_main - M_vec_main(1)))-10^a)/(((t-t_main)+c)^p);
            % If new rate is very small, set it to zero
            if lambda_min_after < 0.0001
                lambda_min_after = 0;
            end
            lambda_null = 0;
    end
end

figure(1)
% Mainshocks
t_ms_vec = t_vec(typ_vec == 1);
m_ms_vec = M_vec(typ_vec == 1);
% Aftershocks
t_as_vec = t_vec(typ_vec == 2);
m_as_vec = M_vec(typ_vec == 2);

scatter(t_as_vec/resol,ones(1,length(t_as_vec)),exp(m_as_vec)/2,'magenta')
hold on
scatter(t_ms_vec/resol,ones(1,length(t_ms_vec)),exp(m_ms_vec),'b','LineWidth',1.5)
scatter(t_ms_vec/resol,2*ones(1,length(t_ms_vec)),exp(m_ms_vec),'b','LineWidth',1.5)
scatter(t_as_vec/resol,3*ones(1,length(t_as_vec)),exp(m_as_vec)/2,'magenta')
xlim([0 50])
yticks([1 2 3])
xlabel('Time [yr]','Interpreter','latex','FontSize',16)
yticklabels({'All','Main Shock','Aftershock'})
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
Max_main = max(M_vec);

% Stochastic_event_set{k} = [typ_vec;M_vec;t_vec/resol];
% 
% end

% save('scenario_MA.mat','Stochastic_event_set');