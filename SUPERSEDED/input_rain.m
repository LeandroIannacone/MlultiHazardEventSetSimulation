clear all
close all
clc

% Accumulated rainfall vs duration

% Parameters of the GEV distribution

% duration = [1,2,4,6,8,12,18,24,48,72,96,120,168];
% mu_R = [50.4,70.7,93.5,111.8,125.9,143.5,167.7,183.2,233.6,261.0,274.1,287.4,308.4];
% sigma_R = [13.6,18.6,33.5,40.7,45.5,48.9,58.1,62.8,74.7,87.3,96.3,103.3,112.6];

% scatter(duration,mu_R)
% hold on
% scatter(duration,sigma_R)
% hold off

% Import IDF from Excel file - Gauge N05 from Tang and Cheung (2011)
Data = xlsread("IDF Curve Data.xlsx");
D_data = Data(3:end,1); % Hours
R_data = Data(3:end,2:end);
T_data = Data(1,2:end);
lambda_data = 1./T_data;
D_data_mat = D_data.*ones(length(D_data),length(T_data));
I_data = R_data./D_data_mat;
lambda_data_mat = lambda_data.*ones(size(D_data_mat));
max_lambda = max(max(lambda_data_mat));

surf(I_data,D_data_mat,lambda_data_mat)
% Set the Z-axis to logarithmic scale
% set(gca, 'ZScale', 'log')

% xlim([0 100])
% ylim([0 100])

D_vec = linspace(0,800,101);
I_vec = linspace(0,350,101);

I_data_new = [I_data,350*linspace(0,1,size(D_data_mat,1))'];
D_data_new = [D_data_mat,ones(size(D_data_mat,1),1)*800];
lambda_data_new = [lambda_data_mat,zeros(size(lambda_data_mat,1),1)];
figure
surf(I_data_new,D_data_new,lambda_data_new)

% Plot multivariate CDF of I and D
CDF_data = (max_lambda - lambda_data_new)/max_lambda;
figure
surf(I_data_new,D_data_new,CDF_data)

figure
plot(D_data,CDF_data(:,end-1))

[I_vecq,D_vecq] = meshgrid(I_vec,D_vec);
% lambda_data_query = interp2(I_data,D_data_mat,lambda_data_mat,I_vecq,D_vecq);

lambda_data_query = zeros(length(I_vec),length(D_vec));
for i = 1:length(I_vec)
    for j = 1:length(D_vec)
        lambda_data_query(i,j) = interp2(I_data,D_data_mat,lambda_data_mat,I_vec(i),D_vec(j));
    end
end