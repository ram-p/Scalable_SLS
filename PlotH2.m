close all
clear
clc
format compact

imax = 15;      % Number of runs
T = 10;         % Time horizon
xnorm = zeros(imax, T+1); unorm = xnorm; obj = xnorm;

disp('Fault-tolerant controller: Iterations:')

for iter = 1:imax
    ctr = 1;    % ctr = 1 is fault-tolerant
    [x, u, K_final, t_fault] = H2(ctr, T);
    disp(iter)
    
    for t = 1:T+1
        xnorm(iter, t) = norm(x(:,t));
        unorm(iter, t) = norm(u(:,t));
        obj(iter, t) = xnorm(iter, t) + unorm(iter, t);
    end

    clearvars -except iter xnorm unorm obj T t_fault imax
end

obj_mean = mean(obj, 1);
obj_std = std(obj);

figure(1)
set(gcf, 'DefaultLineLineWidth', 2.5)
shadedErrorBar(0:T, obj_mean, obj_std, 'lineprops', '-g', 'patchSaturation', 0.33)
hold on
xlabel('Time')
ylabel('$\|x_t\|_{Q}^{2} + \|u_t\|_{R}^{2}$', 'Interpreter', 'latex')
grid on
box on
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 20)

disp('Non-fault-tolerant controller: Iterations:')

for iter = 1:imax
    ctr = 2;    % ctr = 2 is non-fault-tolerant
    [x, u, K_final, t_fault] = H2(ctr, T);
    disp(iter)
    
    for t = 1:T+1
        xnorm(iter, t) = norm(x(:,t));
        unorm(iter, t) = norm(u(:,t));
        obj(iter, t) = xnorm(iter, t) + unorm(iter, t);
    end

    clearvars -except iter xnorm unorm obj T t_fault imax
end

obj_mean = mean(obj, 1);
obj_std = std(obj);

figure(1)
shadedErrorBar(0:T, obj_mean, obj_std, 'lineprops', 'r')
xline(t_fault-1, 'r--', 'LineWidth', 2.5)
legend('Fault tolerant control', 'No fault tolerance', 'Location', 'best')
% ylim([0 5])
% hold on
% xline(t_fault-1, 'r--', 'LineWidth', 2.5)
% xlabel('Time')
% ylabel('$\|x_t\|_{Q}^{2} + \|u_t\|_{R}^{2}$', 'Interpreter', 'latex')
% grid on
% box on