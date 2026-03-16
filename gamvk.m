% Code to test LQG_OF.m for different memory. Careful: for larger T (>10),
% limit k to k_lim = 8. Longer memory uses significant computational
% resources.

close all
clear
clc
format compact

T = 30;
k_lim = 8;
gam = zeros(1, k_lim);
soltim = zeros(1, k_lim);

for k = 1:k_lim
    k
    [gam(k), soltim(k)] = LQG_OF(k, T);
end

figure(1)
set(gcf, 'DefaultLineLineWidth', 2.5)
plot(1:k_lim, gam, 'b');
xlabel('Measurement memory')
ylabel('Cost')
grid on
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 24)

figure(2)
set(gcf, 'DefaultLineLineWidth', 2.5)
plot(1:k_lim, soltim, 'r');
xlabel('Measurement memory')
ylabel('Solve Time (s)')
grid on
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 24)