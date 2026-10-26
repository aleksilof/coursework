%%  HA1 SF2955  %%
%%%%%%%%%%%%%%%%%%

%% Problem 1 & 2
clear all, close all

delta_t = 0.5;
alpha = 0.6;
sigma = 0.5;

Ztates = [0 0; 3.5 0; 0 3.5; 0 -3.5; -3.5 0]';
P = (1/20)*[16 1 1 1 1; 1 16 1 1 1; 1 1 16 1 1; 1 1 1 16 1; 1 1 1 1 16];

Phi_tilde = [1 delta_t 0.5*delta_t^2; 0 1 delta_t; 0 0 alpha];
Phi = [Phi_tilde zeros(3,3); zeros(3,3) Phi_tilde];

Psi_z_tilde = [0.5*delta_t^2; delta_t; 0];
Psi_w_tilde = [0.5*delta_t^2; delta_t; 1];
Psi_z = [Psi_z_tilde zeros(3,1); zeros(3,1) Psi_z_tilde];
Psi_w = [Psi_w_tilde zeros(3,1); zeros(3,1) Psi_w_tilde];

% Simulate trajectory of size m
m = 500;

X = zeros(6, m + 1);
Z_idx = zeros(1, m + 1);

X(:,1) = mvnrnd(zeros(6,1), diag([500,5,5,200,5,5]))';
Z_idx(1) = randi(5);

for n = 1:m
    Z_n = Ztates(:, Z_idx(n));
    W_n = sigma * randn(2, 1);
    X(:, n + 1) = Phi*X(:, n) + Psi_z*Z_n + Psi_w*W_n;
    Z_idx(n + 1) = randsample(5, 1, true, P(Z_idx(n), :));
end

figure;
plot(X(1,:), X(4,:), 'b-');   
hold on;
plot(X(1,1), X(4,1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold off;
xlabel('x (m)');
ylabel('y (m)');
title('Simulated target trajectory');

%% Problem 3
clear all, close all

delta_t = 0.5;
alpha = 0.6;
sigma = 0.5;

Ztates = [0 0; 3.5 0; 0 3.5; 0 -3.5; -3.5 0]';
P = (1/20)*[16 1 1 1 1; 1 16 1 1 1; 1 1 16 1 1; 1 1 1 16 1; 1 1 1 1 16];
Pcum = cumsum(P, 2);

Phi_tilde = [1 delta_t 0.5*delta_t^2; 0 1 delta_t; 0 0 alpha];
Phi = [Phi_tilde zeros(3,3); zeros(3,3) Phi_tilde];

Psi_z_tilde = [0.5*delta_t^2; delta_t; 0];
Psi_w_tilde = [0.5*delta_t^2; delta_t; 1];
Psi_z = [Psi_z_tilde zeros(3,1); zeros(3,1) Psi_z_tilde];
Psi_w = [Psi_w_tilde zeros(3,1); zeros(3,1) Psi_w_tilde];

stations = load("stations.mat");
stationx = stations.pos_vec(1,:);
stationy = stations.pos_vec(2,:);
%plot(stationx,stationy, "*r")
pi_mat = zeros(2,6);
for ii = 1:length(stationx)
    pi_mat(:,ii) = [stationx(ii); stationy(ii)];
end

load("RSSI-measurements.mat")
s = 6;
v = 90;
eta = 3;
var_sigma = 1.5;
N = 10000;
m = 500;

X = mvnrnd(zeros(N,6), diag([500,5,5,200,5,5]))';
Z_idx = randi(5, 1, N);

log_p = ln_obs_prob(X, Y(:,1), pi_mat, v, eta, var_sigma, N, s);
log_p = log_p - max(log_p);
w = exp(log_p);
w = w / sum(w);
selected = [0, 10, 20, 50];
W_stored = zeros(N, length(selected));
W_stored(:, 1) = w';

tau1 = zeros(1, m + 1);
tau2 = zeros(1, m + 1);
tau1(1) = w * X(1,:)';   
tau2(1) = w * X(4,:)';   

ESS = zeros(1, m+1);
ESS(1) = 1 / sum(w.^2);


% SIS
for ii = 1:m
    u = rand(N, 1);
  
    Z_n = Ztates(:, Z_idx);
    W_n = sigma * randn(2, N);  

    X = Phi*X + Psi_z*Z_n + Psi_w*W_n;
    Z_idx = sum(u > Pcum(Z_idx, :), 2)' + 1;
    
    log_p = ln_obs_prob(X, Y(:, ii + 1), pi_mat, v, eta, var_sigma, N, s);
    log_w = log(w) + log_p; %SIS
    log_w = log_w - max(log_w);
    w = exp(log_w); 
    w = w / sum(w);
    
    tau1(ii + 1) = w * X(1,:)';
    tau2(ii + 1) = w * X(4,:)';
    ESS(ii + 1) = 1 / sum(w.^2);

    if ismember(ii, selected)
        W_stored(:, selected == ii) = w';
    end
end

% estimated trajectory
figure;
plot(tau1, tau2, 'b-', 'LineWidth', 1.2); hold on;
plot(pi_mat(1,:), pi_mat(2,:), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(tau1(1), tau2(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
xlabel('x^1 (m)'); ylabel('x^2 (m)');
title('SIS estimated trajectory');
legend('Estimated path', 'Base stations', 'Start');

% ESS 
figure;
plot(0:m, ESS);
xlabel('n'); ylabel('ESS');
title('Effective sample size over time');

figure;
for k = 1:length(selected)
    subplot(2, 2, k);
    histogram(W_stored(:,k), 50);
    xlabel('Weight \omega^i');
    ylabel('Count');
    title(['Weights at n = ' num2str(selected(k))]);
end

%% Problem 4
clear all, close all

delta_t = 0.5;
alpha = 0.6;
sigma = 0.5;

Ztates = [0 0; 3.5 0; 0 3.5; 0 -3.5; -3.5 0]';
P = (1/20)*[16 1 1 1 1; 1 16 1 1 1; 1 1 16 1 1; 1 1 1 16 1; 1 1 1 1 16];
Pcum = cumsum(P, 2);

Phi_tilde = [1 delta_t 0.5*delta_t^2; 0 1 delta_t; 0 0 alpha];
Phi = [Phi_tilde zeros(3,3); zeros(3,3) Phi_tilde];

Psi_z_tilde = [0.5*delta_t^2; delta_t; 0];
Psi_w_tilde = [0.5*delta_t^2; delta_t; 1];
Psi_z = [Psi_z_tilde zeros(3,1); zeros(3,1) Psi_z_tilde];
Psi_w = [Psi_w_tilde zeros(3,1); zeros(3,1) Psi_w_tilde];

stations = load("stations.mat");
stationx = stations.pos_vec(1,:);
stationy = stations.pos_vec(2,:);
%plot(stationx,stationy, "*r")
pi_mat = zeros(2,6);
for ii = 1:length(stationx)
    pi_mat(:,ii) = [stationx(ii); stationy(ii)];
end

load("RSSI-measurements.mat")
s = 6;
v = 90;
eta = 3;
var_sigma = 1.5;
N = 10000;
m = 500;

X = mvnrnd(zeros(N,6), diag([500,5,5,200,5,5]))';
Z_idx = randi(5, 1, N);

log_p = ln_obs_prob(X, Y(:,1), pi_mat, v, eta, var_sigma, N, s);
log_p = log_p - max(log_p);
w = exp(log_p);
w = w / sum(w);

selected = [0, 50, 100, 500];
W_stored = zeros(N, length(selected));
W_stored(:, 1) = w';

probably_Z = zeros(1, m + 1);
probably_Z(1) = mode(Z_idx);  % after n=0 

tau1 = zeros(1, m + 1);
tau2 = zeros(1, m + 1);
ESS  = zeros(1, m + 1);
tau1(1) = w * X(1,:)';
tau2(1) = w * X(4,:)';
ESS(1)  = 1 / sum(w.^2);

idx = randsample(N, N, true, w);
X = X(:, idx);
Z_idx = Z_idx(idx);

for ii = 1:m

    u = rand(N, 1);

    Z_n = Ztates(:, Z_idx);
    W_n = sigma * randn(2, N);
    X = Phi*X + Psi_z*Z_n + Psi_w*W_n;
    Z_idx = sum(u > Pcum(Z_idx, :), 2)' + 1;

    log_p = ln_obs_prob(X, Y(:,ii + 1), pi_mat, v, eta, var_sigma, N, s);
    log_p = log_p - max(log_p);
    w = exp(log_p);
    w = w / sum(w);

    tau1(ii + 1) = w * X(1,:)'; 
    tau2(ii + 1) = w * X(4,:)';
    ESS(ii + 1) = 1 / sum(w.^2);

    if ismember(ii, selected)
        W_stored(:, selected == ii) = w';
    end

    idx = randsample(N, N, true, w);
    X = X(:, idx);
    Z_idx = Z_idx(idx);
    
    probably_Z(ii + 1) = mode(Z_idx);

end

figure;
plot(tau1, tau2, 'b-', 'LineWidth', 1.2); 
hold on;
grid on;
plot(pi_mat(1,:), pi_mat(2,:), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(tau1(1), tau2(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
xlabel('x^1 (m)'); ylabel('x^2 (m)');
title('SISR estimated trajectory');
legend('Estimated path', 'Base stations', 'Start');
% 
figure;
plot(0:m, ESS);
xlabel('n'); ylabel('ESS');
title('SISR effective sample size over time');
% 
Ztates_labels = {'No steering (0,0)', 'East (3.5,0)', 'North (0,3.5)', 'South (0,-3.5)', 'West (-3.5,0)'};
figure;
plot(0:m, probably_Z, 'b.');
yticks(1:5);
yticklabels(Ztates_labels);
xlabel('n');
title('Most probable driving command over time');
% 
figure;
for k = 1:length(selected)
     subplot(2, 2, k);
     histogram(W_stored(:,k), 50);
     xlabel('Weight \omega^i');
     ylabel('Count');
     title(['Weights at n = ' num2str(selected(k))]);
end


%% Problem 5
clear all, close all
delta_t = 0.5; 
alpha = 0.6; 
sigma = 0.5;
s = 6; 
v = 90; 
eta = 3; 
N = 10000; 
m = 500;
Ztates = [0 0; 3.5 0; 0 3.5; 0 -3.5; -3.5 0]';
P = (1/20)*[16 1 1 1 1; 1 16 1 1 1; 1 1 16 1 1; 1 1 1 16 1; 1 1 1 1 16];
Pcum = cumsum(P, 2);

Phi_tilde = [1 delta_t 0.5*delta_t^2; 0 1 delta_t; 0 0 alpha];
Phi = blkdiag(Phi_tilde, Phi_tilde);
Psi_z = blkdiag([0.5*delta_t^2; delta_t; 0], [0.5*delta_t^2; delta_t; 0]);
Psi_w = blkdiag([0.5*delta_t^2; delta_t; 1], [0.5*delta_t^2; delta_t; 1]);

stations = load("stations.mat");
pi_mat = stations.pos_vec;

load('RSSI-measurements-unknown-sigma.mat')

zeta_grid = 0.1:0.1:2.9;
ell_hat = zeros(1, length(zeta_grid));

tic;
for jj = 1:length(zeta_grid)
    var_sigma = zeta_grid(jj);
    
    X = mvnrnd(zeros(N,6), diag([500,5,5,200,5,5]))';
    Z_idx = randi(5, 1, N);
    log_normC = zeros(1, m + 1); 
    
    log_p = ln_obs_prob(X, Y(:,1), pi_mat, v, eta, var_sigma, N, s);
    log_normC(1) = log(mean(exp(log_p))); %
    log_p = log_p - max(log_p);
    w = exp(log_p);
    w = w/sum(w);

    idx = randsample(N, N, true, w);
    X = X(:,idx); Z_idx = Z_idx(idx); w = ones(1,N)/N;
    
    for ii = 1:m
        u = rand(N,1);

        Z_n = Ztates(:, Z_idx);
        W_n = sigma * randn(2, N);
        X = Phi*X + Psi_z*Z_n + Psi_w*W_n;
        Z_idx = sum(u > Pcum(Z_idx,:), 2)' + 1;
        
        log_p = ln_obs_prob(X, Y(:,ii + 1), pi_mat, v, eta, var_sigma, N, s);
        log_normC(ii + 1) = log(mean(exp(log_p))); %
        log_p = log_p - max(log_p);
        w = exp(log_p);
        w = w/sum(w);
        
        idx = randsample(N, N, true, w);
        X = X(:,idx); Z_idx = Z_idx(idx);
    end
    
    ell_hat(jj) = mean(log_normC);

end

time = toc;
disp(['Total time: ' num2str(time) ' seconds']);

[~, best_jj] = max(ell_hat);
disp(zeta_grid(best_jj));
%disp(ell_hat)

figure;
plot(zeta_grid, ell_hat, 'b-o');
xlabel('$\varsigma$', 'Interpreter', 'latex');
ylabel('$\ell_m^N(\varsigma)$', 'Interpreter', 'latex');
title('Log-likelihood over grid', 'Interpreter', 'latex');

%% 5 b) results, code from problem 4 with var_sigma = 2.2 (decided in problem5 above)
clear all, close all

delta_t = 0.5;
alpha = 0.6;
sigma = 0.5;

Ztates = [0 0; 3.5 0; 0 3.5; 0 -3.5; -3.5 0]';
P = (1/20)*[16 1 1 1 1; 1 16 1 1 1; 1 1 16 1 1; 1 1 1 16 1; 1 1 1 1 16];
Pcum = cumsum(P, 2);

Phi_tilde = [1 delta_t 0.5*delta_t^2; 0 1 delta_t; 0 0 alpha];
Phi = [Phi_tilde zeros(3,3); zeros(3,3) Phi_tilde];

Psi_z_tilde = [0.5*delta_t^2; delta_t; 0];
Psi_w_tilde = [0.5*delta_t^2; delta_t; 1];
Psi_z = [Psi_z_tilde zeros(3,1); zeros(3,1) Psi_z_tilde];
Psi_w = [Psi_w_tilde zeros(3,1); zeros(3,1) Psi_w_tilde];

stations = load("stations.mat");
stationx = stations.pos_vec(1,:);
stationy = stations.pos_vec(2,:);
%plot(stationx,stationy, "*r")
pi_mat = zeros(2,6);
for ii = 1:length(stationx)
    pi_mat(:,ii) = [stationx(ii); stationy(ii)];
end

load('RSSI-measurements-unknown-sigma.mat')
s = 6;
v = 90;
eta = 3;
var_sigma = 2.2; % from 5
N = 10000;
m = 500;

X = mvnrnd(zeros(N,6), diag([500,5,5,200,5,5]))';
Z_idx = randi(5, 1, N);

log_p = ln_obs_prob(X, Y(:,1), pi_mat, v, eta, var_sigma, N, s);
log_p = log_p - max(log_p);
w = exp(log_p);
w = w / sum(w);

selected = [0, 50, 100, 500];
W_stored = zeros(N, length(selected));
W_stored(:, 1) = w';

probably_Z = zeros(1, m + 1);
probably_Z(1) = mode(Z_idx);  % after n=0 

tau1 = zeros(1, m + 1);
tau2 = zeros(1, m + 1);
ESS  = zeros(1, m + 1);
tau1(1) = w * X(1,:)';
tau2(1) = w * X(4,:)';
ESS(1)  = 1 / sum(w.^2);

idx = randsample(N, N, true, w);
X = X(:, idx);
Z_idx = Z_idx(idx);

for ii = 1:m

    u = rand(N, 1);
    Z_n = Ztates(:, Z_idx);
    W_n = sigma * randn(2, N);
    X = Phi*X + Psi_z*Z_n + Psi_w*W_n;
    Z_idx = sum(u > Pcum(Z_idx, :), 2)' + 1;

    log_p = ln_obs_prob(X, Y(:,ii + 1), pi_mat, v, eta, var_sigma, N, s);
    log_p = log_p - max(log_p);
    w = exp(log_p);
    w = w / sum(w);

    tau1(ii + 1) = w * X(1,:)';
    tau2(ii + 1) = w * X(4,:)';
    ESS(ii + 1) = 1 / sum(w.^2);

    if ismember(ii, selected)
        W_stored(:, selected == ii) = w';
    end

    idx = randsample(N, N, true, w);
    X = X(:, idx);
    Z_idx = Z_idx(idx);
    
    probably_Z(ii + 1) = mode(Z_idx);

end

figure;
plot(tau1, tau2, 'b-', 'LineWidth', 1.2);
grid on;
hold on;
plot(pi_mat(1,:), pi_mat(2,:), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(tau1(1), tau2(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
xlabel('x^1 (m)'); ylabel('x^2 (m)');
title('SISR estimated trajectory');
legend('Estimated path', 'Base stations', 'Start');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Funktioner %%%%%%%%%%%%%%%%%%%%%


function log_p = ln_obs_prob(X, Y, pi_mat, v, eta, var_sigma, N, s)
    pos = X([1,4],:);
    log_p = zeros(1,N);
    for l = 1:s
        dist = sqrt(sum((pos - pi_mat(:,l)).^2 , 1));
        mu = v - 10*eta*log10(dist);
        log_p = log_p - 0.5*((Y(l) - mu)/var_sigma).^2 - log(sqrt(2*pi)*var_sigma);
    end
end 







