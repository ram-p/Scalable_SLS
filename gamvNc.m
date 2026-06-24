% Code to test LQG_MP.m for different number of clusters, i.e.,
% controllers. Memory is limited to a certain number, but could also be a
% tuning parameter.

close all
clear
clc
format compact

T = 10;
k = 5;
Nc_lim = 4;
gam = zeros(1, Nc_lim);
soltim = zeros(1, Nc_lim);

for i = 1:Nc_lim
    i
    [gam(i), soltim(i)] = LQG_MP(k, T, i);
end

figure(1)
set(gcf, 'DefaultLineLineWidth', 2.5)
plot(1:Nc_lim, gam, 'b');
xlabel('Number of clusters')
ylabel('Cost')
grid on
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 24)

figure(2)
set(gcf, 'DefaultLineLineWidth', 2.5)
plot(1:Nc_lim, soltim, 'r');
xlabel('Number of clusters')
ylabel('Solve Time (s)')
grid on
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 24)