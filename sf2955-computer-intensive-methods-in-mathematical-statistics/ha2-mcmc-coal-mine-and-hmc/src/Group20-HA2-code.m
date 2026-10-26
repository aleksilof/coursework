%%%%% HA2 SF2955 %%%%%

%% Problem 1
clear all, close all, clc

tau = load("coal-mine.csv");

tstart = 1851;
tend = 1963;

t_breakp = 1; % VARIERA, NUMBER OF BREAKPOINTS
d = t_breakp + 1; % NUMBER OF INTERVALS


rho = 0.05; % VARIERA, 

N = 10000;

for var_theta = [0.1, 1, 10, 100]
    [theta_vec, lambda_vec, t_vec, accept] = hybrid_MCMC(N, tau, tstart, tend, d, var_theta, rho);

    figure(Name=['var_theta = ', num2str(var_theta)], NumberTitle="off"); 
    subplot(3,1,1); plot(theta_vec); xlabel('Iteration'); ylabel('$\theta$', 'interpreter','latex','fontsize',14); title('theta trace')
    subplot(3,1,2); plot(lambda_vec); xlabel('Iteration'); ylabel('$\lambda$', 'interpreter','latex','fontsize',14); title('lambda trace')
    subplot(3,1,3); plot(t_vec); xlabel('Iteration'); ylabel('Year'); title('breakpoints trace')

    fprintf('\n\n<strong>var_theta = %.1f\n</strong>', var_theta)
    acceptance_rate = accept / N;
    fprintf('Acceptance rate for rho = %.2f: %.2f %%\n', rho, acceptance_rate*100)

    for ii = 1:t_breakp
        mean_t_vec = mean(t_vec);
        fprintf('Mean breakpoint %.0f: %.2f\n', ii, mean_t_vec(ii))
    end
    fprintf('Mean theta: %.2f\n', mean(theta_vec))
    for ii = 1:d
        mean_lambda_vec = mean(lambda_vec);
        fprintf('Mean lambda %.0f: %.2f\n', ii, mean_lambda_vec(ii))
    end
    
end

var_theta = 1; % HYPERPARAMETER

for rho = [0.01, 0.05, 0.1, 0.5, 1, 2]
    [theta_vec, lambda_vec, t_vec, accept] = hybrid_MCMC(N, tau, tstart, tend, d, var_theta, rho);

    figure(Name=['rho = ', num2str(rho)], NumberTitle="off"); 
    subplot(3,1,1); plot(theta_vec); xlabel('Iteration'); ylabel('$\theta$', 'interpreter','latex','fontsize',14); title('theta trace')
    subplot(3,1,2); plot(lambda_vec); xlabel('Iteration'); ylabel('$\lambda$', 'interpreter','latex','fontsize',14); title('lambda trace')
    subplot(3,1,3); plot(t_vec); xlabel('Iteration'); ylabel('Year'); title('breakpoints trace')

    fprintf('\n\n<strong>rho = %.2f\n</strong>', rho)
    acceptance_rate = accept / N;
    fprintf('Acceptance rate for rho = %.2f: %.2f %%\n', rho, acceptance_rate*100)

    for ii = 1:t_breakp
        mean_t_vec = mean(t_vec);
        fprintf('Mean breakpoint %.0f: %.2f\n', ii, mean_t_vec(ii))
    end
    fprintf('Mean theta: %.2f\n', mean(theta_vec))
    for ii = 1:d
        mean_lambda_vec = mean(lambda_vec);
        fprintf('Mean lambda %.0f: %.2f\n', ii, mean_lambda_vec(ii))
    end
    
end
%% Problem 2
clear all, close all, clc
y = load("hmc-observations.csv");
sigma = 2;
Sigmat = [5 0; 0 0.5];
U = @(theta) 0.5/sigma^2 * sum((y - theta(1)^2 - theta(2)^2).^2) + 0.5 * theta * (Sigmat \ theta');
H = @(theta, v) U(theta) + 0.5 * sum(v.^2);
N = 20000;
thetavec = zeros(N, 2);
thetavec(1, :) = [0 -1];
L = 1;
epsilon = 0.05;

accept_HMC = 0;
for n = 2:N
    V = randn(1, 2);
    [theta, v] = leapfrog(thetavec(n - 1, :), V, sigma, Sigmat, y, epsilon, L);
    H_proposed = H(theta, v);
    H_current = H(thetavec(n - 1, :), V);
    alpha = min(1, exp(H_current - H_proposed));
    if rand <= alpha
        thetavec(n, :) = theta;
        accept_HMC =  accept_HMC + 1;
    else
        thetavec(n, :) = thetavec(n - 1, :);
    end
end

acceptance_rate_HMC = accept_HMC / N;
disp("Avg acceptance rate for HMC")
disp(round(acceptance_rate_HMC*100, 2))

% scatter plot
figure;
scatter(thetavec(:,1), thetavec(:,2), 1, 'filled')
xlabel('\theta_1'); ylabel('\theta_2')
title('HMC samples')

figure;
histogram2(thetavec(:,1), thetavec(:,2), 50, 'DisplayStyle', 'tile', 'ShowEmptyBins', 'on')
xlabel('\theta_1'); ylabel('\theta_2')
title('HMC 2D histogram')
colorbar

figure;
[counts, xedges, yedges] = histcounts2(thetavec(:,1), thetavec(:,2), 50);
xcenters = (xedges(1:end-1) + xedges(2:end)) / 2;
ycenters = (yedges(1:end-1) + yedges(2:end)) / 2;
surf(xcenters, ycenters, counts')
xlabel('\theta_1'); ylabel('\theta_2'); zlabel('counts')
title('HMC 3D surface')
colorbar
shading interp

% MH
zeta        =  0.4;
thetaMH     = zeros(N, 2);
thetaMH(1,:) = [0 -1];

acceptMH    = 0;

for n = 2:N
    theta_star = thetaMH(n - 1,:) + zeta * randn(1,2);
    ln_alpha = ln_posterior(theta_star, y, sigma, Sigmat) - ln_posterior(thetaMH(n-1,:), y, sigma, Sigmat);  
    if log(rand) <= ln_alpha
        thetaMH(n,:) = theta_star;
        acceptMH = acceptMH + 1;
    else
        thetaMH(n,:) = thetaMH(n-1,:);
    end
end

acceptance_rate_MH = acceptMH / N;
fprintf('MH acceptance rate: %.4f\n', round(acceptance_rate_MH*100, 2))

% 2D histogram
figure;
histogram2(thetaMH(:,1), thetaMH(:,2), 50, 'DisplayStyle', 'tile', 'ShowEmptyBins', 'on')
xlabel('\theta_1'); ylabel('\theta_2')
title('MH 2D histogram')
colorbar

% scatter with trajectory
figure;
scatter(thetaMH(:,1), thetaMH(:,2), 1, 'filled')
xlabel('\theta_1'); ylabel('\theta_2')
title('MH samples')


% 3D surface
figure;
[counts, xedges, yedges] = histcounts2(thetaMH(:,1), thetaMH(:,2), 50);
xcenters = (xedges(1:end-1) + xedges(2:end)) / 2;
ycenters = (yedges(1:end-1) + yedges(2:end)) / 2;
surf(xcenters, ycenters, counts')
xlabel('\theta_1'); ylabel('\theta_2'); zlabel('counts')
title('MH 3D surface')
colorbar
shading interp

figure;
subplot(2,1,1); 
acf = xcorr(thetavec(:,1) - mean(thetavec(:,1)), 50, 'normalized');
stem(0:50, acf(51:end));
title('HMC ACF \theta_1'); xlabel('lag'); ylabel('ACF')

subplot(2,1,2);
acf = xcorr(thetaMH(:,1) - mean(thetaMH(:,1)), 50, 'normalized');
stem(0:50, acf(51:end));
title('MH ACF \theta_1'); xlabel('lag'); ylabel('ACF')

figure;
subplot(2,1,1); 
acf = xcorr(thetavec(:,2) - mean(thetavec(:,2)), 50, 'normalized');
stem(0:50, acf(51:end));
title('HMC ACF \theta_2'); xlabel('lag'); ylabel('ACF')

subplot(2,1,2);
acf = xcorr(thetaMH(:,2) - mean(thetaMH(:,2)), 50, 'normalized');
stem(0:50, acf(51:end));
title('MH ACF \theta_2'); xlabel('lag'); ylabel('ACF')

%%
%%%%%%%%% Functions %%%%%%%%%

function [theta_T, v_T] = leapfrog(theta_0, v_0, sigma, Sigmat, y, epsilon, L)
    grad_U = @(theta) -2*theta*sum(y- theta(1)^2-theta(2)^2)/sigma^2 + (Sigmat \ theta')';
    theta_T = theta_0;
    v_T = v_0 - 0.5*epsilon*grad_U(theta_T); % v_1/2 = v_0 - 0.5*epsilon*grad_U(theta_0)
    for m = 1:L
        theta_T = theta_T + epsilon*v_T; % theta_t+1 = theta_t + epsilon*v_t+1/2 (theta_L = theta_L-1 + epsilon*v_L-1/2 last step)
        if m~= L
            v_T = v_T - epsilon*grad_U(theta_T); % v_t+3/2 = v_t+1/2 - epsilon*grad_U(heta_t+1) 
        end
    end
    v_T = v_T - 0.5*epsilon*grad_U(theta_T); % v_L = v_L-1/2 - 0.5*epsilon*grad_U(theta_L)
    v_T = -v_T;
end


function ln_f = ln_posterior(theta, y, sigma, Sigmat)
    ln_f = -0.5/sigma^2 * sum((y - theta(1)^2 - theta(2)^2).^2) - 0.5 * theta * (Sigmat \ theta');
end


function [theta_vec, lambda_vec, t_vec, accept] = hybrid_MCMC(N, tau, tstart, tend, d, var_theta, rho)

    tmid = linspace(tstart, tend, d + 1);
    t = [tstart, tmid(2:d), tend];
    
    lambda = ones(1, d);
    theta_vec = zeros(N, 1);
    lambda_vec = zeros(N, d);
    t_vec = zeros(N, d - 1);
    
    accept = 0;
    
    for ii = 1:N
        % Gibbs where possible
        theta = gamrnd(2*d + 2, 1/(var_theta + sum(lambda)));
        % for each lambda i
        for jj = 1:d
            ni = sum(tau >= t(jj) & tau < t(jj + 1));
            lambda(jj) = gamrnd(ni + 2, 1/(theta + t(jj + 1) - t(jj)));
        end
        % Metropolis-Hastings
        for kk = 2:d
            R = rho*(t(kk + 1) - (t(kk - 1)));
            t_star = t(kk) + unifrnd(-R, R);
            if t_star > t(kk - 1) && t_star < t(kk + 1)
                t_proposal = t;
                t_proposal(kk) = t_star;
                n_old = arrayfun(@(x) sum(tau >= t(x) & tau < t(x + 1)), 1:d);
                n_new = arrayfun(@(x) sum(tau >= t_proposal(x) & tau < t_proposal(x + 1)), 1:d);
                
                num = exp(-sum(lambda.*diff(t_proposal)))*prod(lambda.^n_new)*prod(diff(t_proposal));
                den = exp(-sum(lambda.*diff(t)))*prod(lambda.^n_old)*prod(diff(t));
                alpha = min(1, num/den);
                if rand <= alpha
                    t(kk) = t_star;
                    accept = accept + 1;
                end
            end
        end
        theta_vec(ii) = theta;
        lambda_vec(ii,:) = lambda; 
        t_vec(ii,:) = t(2:d);
    end

end