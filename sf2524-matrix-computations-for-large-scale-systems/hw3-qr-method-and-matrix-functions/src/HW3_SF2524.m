%%%%% HW3 %%%%%
% Aleksi Löf, Anton Blomgren

%% 1
clc, clear all, close all


errfun = @(A) max(max(abs(tril(A,-1))));

alpha_max = 100000;

tol = 1e-10;

alpha_values = 10:50:100000;
iterations = zeros(size(alpha_values));
predict_iter = zeros(size(alpha_values));

for jj = 1:length(alpha_values)
    alpha = alpha_values(jj);
    A = alpha_example(alpha);
    
    err = errfun(A);
    k = 0;
    while err > tol && k < 5000
        [Q,R] = qr(A);
        A = R*Q;
        k = k + 1;
        err = errfun(A);
    end
    eigvals = eig(A);
    iterations(jj) = k;
    predict_iter(jj) = log(tol)/log(abs(eigvals(2)/eigvals(1)));
end

figure;
semilogx(alpha_values, iterations, 'r-','LineWidth',1.5);
ylim([0 3000])
xlim([10 100000])
xlabel('\alpha');
ylabel('Number of iterations');
grid on
hold on
semilogx(alpha_values, predict_iter, 'k-.', 'LineWidth', 1.5);
legend('Number of iterations','Predicted number of iterations')

%% 2
clc, clear all, close all
mvec = [10 100 200 300 400];
tvecb = zeros(length(mvec), 1);
tvecalg2 = zeros(length(mvec), 1);
k = 1;

for m = mvec
    A = alpha_example(1,m);
    fprintf('For Naive Hessenberg m = %d\n', m);
    tic;
    H1 = naive_hessenberg_red(A);
    time_elapsed_1 = toc;
    tvecb(k) = time_elapsed_1;
    fprintf('Time: %.6f seconds\n\n', time_elapsed_1);

    fprintf('For Algorithm2 m = %d\n', m);
    tic;
    H2 = Algorithm2(A);
    time_elapsed_2 = toc;
    tvecalg2(k) = time_elapsed_2;
    fprintf('Time: %.6f seconds\n\n', time_elapsed_2);
    disp('------------------------------')
    fprintf('Condition number: %.6e\n', cond(A));
    disp('------------------------------')
    k = k + 1;
end

max(eig(A))
max(eig(H1))
max(eig(H2))
plot(mvec,tvecb,'r-o')
hold on
plot(mvec,tvecalg2,'b-o')
legend('Naive Hessenberg','Algorithm 2')
title('Time comparison')
xlabel('m values')
ylabel('Time')

X = hess(A);

%% 2 d)
clc, clear all, close all
format long

epsilon_vec =  [0.4 0.1 0.01 0.001 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 0];
sigma_0_vec = zeros(length(epsilon_vec), 1);
sigma_a22_vec = zeros(length(epsilon_vec), 1);
k = 1;
for epsilon = epsilon_vec
    A = [3 2; epsilon 1];

    sigma0 = 0;
    [Q,R] = qr(A - sigma0*eye(2));
    A_tilde0 = R*Q + sigma0*eye(2);
    sigma_0_vec(k) = abs(A_tilde0(2,1));

    sigma_a22 = A(2,2);
    [Q,R] = qr(A - sigma_a22*eye(2));
    A_tilde = R*Q + sigma_a22*eye(2);
    sigma_a22_vec(k) = abs(A_tilde(2,1));

    k = k + 1;
end

disp('\sigma = 0')
disp(sigma_0_vec)

disp('\sigma = a_2,2')
disp(sigma_a22_vec)

figure;
semilogy(1:12, sigma_0_vec, 'r-o', 'DisplayName', '\sigma = 0');
hold on;
semilogy(1:12, sigma_a22_vec, 'b-s', 'DisplayName', '\sigma = a_{22}');
legend('Sigma = 0', 'Sigma = a_{2,2}')
hold off

%% 3
clc, clear all, close all

A = [1 4 4; 3 -1 3; -1 4 4];

f = @(x) sin(x);
sinA = schur_parlett(A, f);

disp('sin(A) = ');
disp(sinA);

%% 3 b)
clc, clear all, close all

A = rand(100,100) + 1i*rand(100,100);
A = A/norm(A);
Nvec = [50 100 150 200 250 500 1000 1500 2000 4000 6000];

time_naive_vec = zeros(length(Nvec),1);
time_schur_vec = zeros(length(Nvec),1);

k = 1;
for N = Nvec
    B = A;
    tic;
    for ii = 1:N - 1
        B = B*A;
    end
    time_naive_vec(k) = toc;

    g = @(x) x^N;
    tic;
    schur_parlett(A, g);
    time_schur_vec(k) = toc;

    k = k + 1;
end

figure
plot(Nvec,time_naive_vec,'k*')
hold on
plot(Nvec,time_schur_vec,'ro')
xlabel('N')
ylabel('CPU-Time')
legend('Naive','Schur-Parlett',Location='northwest')
grid on
hold off

%% 4) 
clc, clear all, close all
eps_vec = logspace(-10, -1, 100);
errors = zeros(length(eps_vec), 1);
k = 1;
for epsilon = eps_vec
    A = [pi 1; 0 pi + epsilon];

    [V, D] = eig(A);

    F = diag(exp(diag(D)));
    FJ = V*F/V;

    A_exact = [exp(pi), exp(pi)*(exp(epsilon)-1)/epsilon; 
           0,       exp(pi + epsilon)];

    errors(k) = norm(FJ - A_exact);
    k = k + 1;
end

figure;
loglog(eps_vec, errors, 'b-', 'LineWidth', 1.5);
xlabel('\epsilon');
ylabel('||f(A) - F||');
grid on;
xlim([1e-10 1e-1])
xticks([10^-10 10^-9 10^-8 10^-7 10^-6 10^-5 10^-4 10^-3 10^-2 10^-1])
%ylim([1e-15 1])
title('Error for different \epsilon');
%% 5
clc, clear all, close all
n = 256;
randn('seed', 0)
A = randn(n, n);
A = A/norm(A);
A = A - 3*eye(n);

table1 = zeros(7, 8);
for m = 1:7
    for j = 0:7
        table1(m, j + 1) = (m - 1) + j;
    end
end

table2 = zeros(7, 8);
A_exp_exact = expm(A);
for m = 1:7
    for j = 0:7
        A_scaled = A / (2^j);
        Tm = taylor(A_scaled, m);
        for k = 1:j
            Tm = Tm * Tm;
        end
        table2(m, j + 1) = log10(norm(Tm - A_exp_exact));
    end
end

table3 = zeros(7, 8);
for m = 1:7
    for j = 0:7
        tic;
        A_scaled = A / (2^j);
        Tm = taylor(A_scaled, m);
        for k = 1:j
            Tm = Tm * Tm;
        end
        table3(m, j + 1) = toc;
    end
end

disp('Table 1: Number of matrix-matrix products');
disp(table1);
disp('Table 2: log10(||P_jm - exp(A)||)');
disp(table2);
disp('Table 3: CPU time');
disp(table3);

%% 7
clc, clear all, close all
load('mysterious_images.mat')
W = 68; H = 82;

%for ii = 1:288
%    x = X(:,ii);
%    myimage = reshape(x,W,H);
%    imagesc(myimage);
%    drawnow;
%end

A = create_adjacency(X, 0);

timevec_expm = zeros(100,1);
timevec_taylor = zeros(100,1);

m = 2;
j = 1;


for ii = 1:100
    tic;
    expA_exact = expm(A);
    estrada_exact = trace(expA_exact);
    centrality_exact = diag(expA_exact);
    timevec_expm(ii) = toc;

    tic;
    A_scaled = A / (2^j);
    Tm = taylor(A_scaled, m);
    for k = 1:j
        Tm = Tm * Tm;
    end

    estrada_taylor = trace(Tm);
    centrality_taylor = diag(Tm);
    timevec_taylor(ii) = toc;

end



fprintf('Avg Time with MATLAB expm: %.3e\n', mean(timevec_expm));
fprintf('Avg Time with Taylor method: %.3e\n', mean(timevec_taylor));


error_matrix = log10(norm(Tm - expA_exact));
error_estrada = log10(abs(estrada_taylor - estrada_exact));

fprintf('error matrix norm: %.3e\n', error_matrix);
fprintf('error Estrada: %.3e\n', error_estrada);

[~, rank_exact] = sort(centrality_exact, 'descend');
[~, rank_taylor] = sort(centrality_taylor, 'descend');

first_diff = find(rank_exact ~= rank_taylor, 1);

fprintf('\nFirst image where rankings differ: position %d\n', first_diff);
fprintf('  Exact method ranks image %d at position %d\n', rank_exact(first_diff), first_diff);
fprintf('  Taylor method ranks image %d at position %d\n', rank_taylor(first_diff), first_diff);

cut_off_vec = [0 0.5 0.6 0.7 0.9];
k = 1;

figure;
for ii = cut_off_vec
    A = create_adjacency(X, ii);
    subplot(2,3,k)
    plot(graph(A))
    title(sprintf('Cutoff: %.1f', ii))
    k = k + 1;
end

