% Mat_comb_plot = [0.1040,0.9608,0.3459,2.5454;
%     0.1508,2.0008,0.5012,2.9995;
%     0.5675,0.6973,1.8743,0.5792;
%     0.1174,1.2988,0.3859+0.3111,0.0625];
% Mat_comb_plot = log(Mat_comb_fin_mean);
Mat_comb_plot = log(Mat_comb_fin_median);


Mat_comb_plot(:,1) = Vec_count_fin_median;
Mat_comb_plot(:,2) = Vec_count_fin_median;
Mat_comb_plot(:,4) = Vec_count_fin_median;
Mat_comb_plot(:,5) = Vec_count_fin_median;
Mat_comb_plot(3,:) = [];
Mat_comb_plot(:,3) = [];

x_vec = 0:0.01:3.99;
y_vec = 0:0.01:3.99;
x_flo = floor(x_vec) + 1;
y_flo = floor(y_vec) + 1;
z = zeros(length(x_vec),length(y_vec));
for i = 1:length(x_vec)
    for j = 1:length(x_vec)
        haz1 = x_flo(i);
        haz2 = y_flo(j);
        z(i,j) = Mat_comb_plot(haz1,haz2);
    end
end

xs = surf(x_vec,y_vec,z);
colormap(flipud(copper))
shading interp