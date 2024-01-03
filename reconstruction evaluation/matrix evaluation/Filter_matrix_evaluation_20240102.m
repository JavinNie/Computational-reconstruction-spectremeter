clear;
clc;
load('200_0.3.mat')

% 计算矩阵 T 的条件数
condition_number = cond(T);
disp(['矩阵 T 的条件数为：' num2str(condition_number)]);

% 计算每行的自相关系数
autocorr_means = zeros(size(T, 1), 1);

for i = 1:size(T, 1)
    row_autocorr = xcorr(T(i, :), 'coeff');
    autocorr_means(i) = mean(row_autocorr);
end

disp('每行的自相关系数平均值：autocorr_means');
disp(autocorr_means);

% 计算每行之间的互相关系数
inter_row_correlation = corrcoef(T');
disp('每行之间的互相关系数矩阵：inter_row_correlation');
disp(inter_row_correlation);
