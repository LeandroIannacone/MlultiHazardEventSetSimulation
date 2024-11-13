clear all
close all
clc

% Import IDF from Excel file - Gauge N05 from Tang and Cheung (2011)
Data = xlsread("IDF Curve Data.xlsx");
D_data = Data(3:end,1); % Hours
R_data = Data(3:end,2:end);
T_data = Data(1,2:end);

figure
surf(T_data,D_data,R_data)
R_vec = linspace(min(min(R_data)),max(max(R_data)),1000);

T_mat = zeros(length(D_data),length(R_vec));
for i = 1:length(D_data)
    for j = 1:length(R_vec)

        R_data_temp = R_data(i,:);
        index = find(R_data_temp >= R_vec(j), 1, 'first');
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

figure
surf(R_vec,D_data,1./T_mat)

%%
lambda_data = 1./T_data;
D_data_mat = D_data.*ones(length(D_data),length(T_data));
I_data = R_data./D_data_mat;
lambda_data_mat = lambda_data.*ones(size(D_data_mat));
max_lambda = max(max(lambda_data_mat));

figure
surf(T_data,D_data,R_data)
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

figure
surf(I_vec,D_data,1./T_mat)

%% Find empirical marginal CDF for Intensity
Vec_I = [I_data(:,9),ones(19,2).*I_data(:,8),ones(19,5).*I_data(:,7),ones(19,10).*I_data(:,6),ones(19,20).*I_data(:,5),ones(19,50).*I_data(:,4),ones(19,100).*I_data(:,3),ones(19,200).*I_data(:,2),ones(19,500).*I_data(:,1)];
Vec_I = reshape(Vec_I,[],1);

[I_cdf,I_vec_cdf] = ecdf(Vec_I);
figure
plot(I_vec_cdf,I_cdf)

%% Simulate events for rainfalls
U_I = rand;
MI = interp1(I_cdf,I_vec_cdf,U_I)
[Ig,Dg] = meshgrid(I_vec,D_data);
lambda_d_temp = interp2(I_vec,D_data,1./T_mat,MI,D_data);
% Make vector unique
index1 = find(lambda_d_temp == 0, 1, 'last');
index2 = find(lambda_d_temp == 0.5, 1, 'first');
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
% if min(lambda_d_temp2) > 0
%     lambda_d_temp2 = [0;lambda_d_temp2];
%     D_data_temp = [745;D_data_temp];
% end
U_D = rand;
lambda_d = max(lambda_d_temp2)*(1-U_D);
MD = interp1(lambda_d_temp2,D_data_temp,lambda_d)
