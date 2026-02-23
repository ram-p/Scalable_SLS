close all
clear
clc
format compact

imax = 10;      % Number of runs
T = 10;         % Time horizon
xnorm = zeros(imax, T+1);

disp('Fault-tolerant controller: Iterations:')

for iter = 1:imax
    ctr = 1;    % ctr = 1 is fault-tolerant
    [x, u, K_final, t_fault, N] = L1(ctr, T);
    disp(iter)

    for t = 1:T+1
        xnorm(iter, t) = norm(x(:,t), inf);
    end

    clearvars -except iter xnorm T t_fault N imax
end

obj_mean = max(xnorm, [], 1);
obj_std = std(xnorm);

figure(1)
set(gcf, 'DefaultLineLineWidth', 2.5)
shadedErrorBar(0:T, obj_mean, obj_std, 'lineprops', '-g', 'patchSaturation', 0.33)
hold on
xlabel('Time')
ylabel('$\|x_t\|_{\infty}$', 'Interpreter', 'latex')
grid on
box on
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 20)

disp('Non-fault-tolerant controller: Iterations:')

for iter = 1:imax
    ctr = 3;    % ctr = 3 is memoryless
    [x, u, K_final, t_fault, N] = L1(ctr, T);
    disp(iter)

    for t = 1:T+1
        xnorm(iter, t) = norm(x(:,t), inf);
    end

    clearvars -except iter xnorm T t_fault N
end

obj_mean = max(xnorm, [], 1);
obj_std = std(xnorm);

figure(1)
shadedErrorBar(0:T, obj_mean, obj_std, 'lineprops', 'r')
xline(t_fault-1, 'r--', 'LineWidth', 2.5)
yline(max(N), 'k--', 'LineWidth', 2.5)
legend('With memory', 'Without memory', 'Location', 'best')