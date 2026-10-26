%%% HW 2 %%%

%% 1)

%%% a) b) 
clear all, close all, clc
%alpha = 5;
%n = 100;
%rand('state',5);
%A = sprand(n,n,0.5);
%A = A + alpha*speye(n);
%A = A/norm(A,1);
%b = rand(n,1);
n = 100;

%x_tilda = GMRES(A,b,m,n);

alphavec = [1,5,10,100];
idx = 1;
for ii = alphavec
    x_rel = zeros(100, 1);
    rel_res = zeros(100, 1);
    rand('state',5);
    A = sprand(n,n,0.5);
    A = A + ii*speye(n);
    A = A/norm(A,1);
    b = rand(n,1);
    x_analytic = A\b;
    for m = 1:100
        x_tilda = GMRES(A,b,m);
        x_rel(m) = norm(x_analytic - x_tilda)/norm(x_analytic);
        rel_res(m) = norm(A*x_tilda - b)/norm(b);
    end

    figure(1)
    subplot(2, 2, idx);
    eigenvalues = eig(full(A));
    plot(eigenvalues,'k*')
    title(sprintf('\\alpha = %d',ii))
    xlabel('Im \lambda')
    ylabel('Re \lambda')
    grid on
    xlim([-0.2, 1.1])
    ylim([-0.2, 0.2])

    center = mean(eigenvalues(2:end));
    radius = max(abs(eigenvalues(2:end)-center));
    theta = linspace(0, 2*pi, 100);
    circle_x = real(center) + radius*cos(theta);
    circle_y = imag(center) + radius*sin(theta);
    hold on
    plot(circle_x,circle_y,'-')
    legend('Eigenvalues','Localization disk',location = 'northeast')
    hold off

    figure(2)
    subplot(2, 2, idx);
    %semilogy(x_tilda)
    semilogy(rel_res,'r')
    hold on
    semilogy(x_rel,'k--')
    xlabel('Iteration')
    ylabel('Error')
    title(sprintf('\\alpha = %d',ii))

    rho = radius / abs(center);
    m_val = 1:100;
    th_conv = rho.^m_val;

    semilogy(th_conv,'b--')

    legend('$\|Ax_m - b\| / \|b\|$','$\|x_\star - x_m\| / \|x_\star\|$','Theory prediction', 'Interpreter', 'latex')
    
    grid on
    ylim([1e-17, 10])
    hold off

    idx = idx + 1;
end

%% 1 c)
clear all, close all, clc
alpha = [1,100];
nvec = [200, 500, 1000];
mvec = [5, 10, 20, 50, 100];

for ii = alpha
    disp('a')
    disp(ii)
    for n = nvec
        disp('n')
        disp(n)
        rand('state',5);
        A = sprand(n,n,0.5);
        A = A + ii*speye(n);
        A = A/norm(A,1);
        b = rand(n,1);

        for m = mvec
            disp('m')
            disp(m)
            tic;
            x_tilde = GMRES(A,b,m);
            time_GMRES = toc;
            resnorm = norm(A*x_tilde-b)
        end
        tic;
        x_backslash = A\b;
        time_backs = toc;
        resnorm_backs = norm(A*x_backslash - b)
    end
end

%% 2 
clc, clear, close all
n = 100;
e = ones(n, 1);
A = -spdiags([e -6*e e], -1:1, n, n); % Symmetric
b = kron(ones(n/2, 1), [0;1]) + ones(n, 1);
m = 30;
% CG
[~, ~, res_norms_cg, ~] = cg(A, b, m, 1);
% GMRES
res_norms_GMRES = zeros(m, 1);
for k = 1:m
    x_approx_GMRES = GMRES(A, b, k);
    res_norms_GMRES(k) = norm(A*x_approx_GMRES - b);
end
% CG with LSQ
% Normal equation: z = ((Ax)'(Ax))^−1 * (Ax)'b <=> (x'A'Ax)^-1 * x'A'*b,
% A symmetric => (x'x)^-1 * x'A'*b <=> (x^-1)(x')^-1(x')A'b <=> (x^-1)A'b
[~, ~, res_norms_cg_LSQ, ~] = cg_LSQ(A, b, m, 1);
subplot(1, 2, 1)
semilogy(1:30, res_norms_cg, 'k*')
hold on
semilogy(1:30, res_norms_GMRES, 'bo')
hold on
semilogy(1:30, res_norms_cg_LSQ, 'r+')
grid on
ylim([1e-17 1e3])
xlim([0 30])
xlabel('Iteration $p$', 'Interpreter', 'latex')
ylabel('$\|Ax - b\|_2$', 'Interpreter', 'latex')
legend('Method 1: CG', 'Method 2: GMRES', 'Method 3: CG with LSQ')
subplot(1, 2, 2)
plot(1:30, res_norms_cg, 'k*')
hold on
plot(1:30, res_norms_GMRES, 'bo')
grid on
plot(1:30, res_norms_cg_LSQ, 'r+')
grid on
ylim([0.1528 0.1532])
xlim([1.6 2.4])
xlabel('Iteration $p$', 'Interpreter', 'latex')
ylabel('$\|Ax - b\|_2$', 'Interpreter', 'latex')
legend('Method 1: CG', 'Method 2: GMRES', 'Method 3: CG with LSQ')
set(findall(gcf,'-property','FontSize'), 'FontSize', 12)
set(findall(gcf,'Type','Line'),'LineWidth',1)
set(findall(gcf,'Type','Line'),'MarkerSize',8)


%% 5
clc, clear, close all
m = 80;

load('cgn_illustration.mat') %%% Ladda in nya sen!!!

[~, ~, res_norms_cgne, ~] = cgne(A, b, m, 1);
[~, ~, ~, t_vec_CGN] = cgne(A, b, m, 1);


res_norms_GMRES_5 = zeros(m, 1);

t_vec_GMRES = zeros(m,1);
t0 = tic;
for k = 1:m
    x_approx_GMRES_5 = GMRES(A, b, k);
    t_vec_GMRES(k) = toc(t0);
end

for k = 1:m
    x_approx_GMRES_5 = GMRES(A, b, k);
    res_norms_GMRES_5(k) = norm(A*x_approx_GMRES_5 - b);
end

figure
semilogy(1:m, res_norms_cgne, 'r-')
hold on
semilogy(1:m, res_norms_GMRES_5, 'b--')
grid on
ylim([1e-18 1e3])
xlim([0 80])
xlabel('Iteration', 'Interpreter', 'latex')
ylabel('$\|Ax - b\|_2$', 'Interpreter', 'latex')
legend('Method 1: CGN', 'Method 2: GMRES')
hold off

figure
semilogy(t_vec_GMRES,res_norms_GMRES_5,'b-')
hold on
semilogy(t_vec_CGN,res_norms_cgne,'r')
grid on
xlim([0 45])
ylim([1e-18 1e0])
xlabel('CPU time in seconds', 'Interpreter', 'latex')
ylabel('$\|Ax - b\|_2$', 'Interpreter', 'latex')
legend('Method 1: GMRES', 'Method 2: CGN')
hold off

%% 6 a)
clc, clear, close all

data1 = load('heatsink_exercise_ref1.mat');
data2 = load('heatsink_exercise_ref2.mat');

% a
t_vec1 = zeros(10,1);
for ii = 1:10
    [x1, ~, A1] = heatsink_build_matrices(data1);
    tic;
    b = A1*x1;
    t_vec1(ii) = toc;
end
avg_t_1 = mean(t_vec1);
disp('Average time for matrix vector product A, model X = 1 : ')
disp(avg_t_1)


t_vec2 = zeros(10,1);
for ii = 1:10
    [x2, ~, A2] = heatsink_build_matrices(data2);
    tic;
    b = A2*x2;
    t_vec2(ii) = toc;
end
avg_t_2 = mean(t_vec2);
disp('Average time for matric vector product A, model X = 2 : ')
disp(avg_t_2)

% b
n1 = length(x1);
figure
plot_heatsink(x1,'heatsink_exercise_ref1.mat')
figure
plot_heatsink(293*ones(n1,1),'heatsink_exercise_ref1.mat')



% c
t_vec1 = zeros(10,1);
for ii = 1:10
    [x1, t1, A1] = heatsink_build_matrices(data1);
    t_vec1(ii) = t1;
end
avg_t_1 = mean(t_vec1);
disp('Average time for backslash with model X = 1 : ')
disp(avg_t_1)

t_vec2 = zeros(10,1);
for ii = 1:10
    [x2, t2, A2] = heatsink_build_matrices(data2);
    t_vec2(ii) = t2;
end
avg_t_2 = mean(t_vec2);
disp('Average time for backslash with model X = 2 : ')
disp(avg_t_2)


%% 6 d)
clc, clear, close all
% d
N = 100;
data1_exercise1 =  load('heatsink_exercise_ref1.mat');
data2_exercise2 = load('heatsink_exercise_ref2.mat');
data1 = load('heatsink_precond_ref1.mat');
data2 = load('heatsink_precond_ref2.mat');

T1a = data1.T1;
T2a = data1.T2;

T1b = data2.T1;
T2b = data2.T2;

% data 1
n1 = length(T1a);
[x1, t1, A1, ~] = heatsink_build_matrices(data1_exercise1);
D = spdiags(diag(A1), 0, n1, n1);
r = zeros(N,1);
d = zeros(N,1);

t_vec1r = zeros(N,1);
for ii = 1:N
    z1 = rand(n1,1);
    tic;
    r(ii) = norm(T2a\(T1a\(A1*z1))-z1)/norm(A1*z1);
    t_vec1r(ii) = toc;
end
avg_t_1r = mean(t_vec1r);
disp('Average time for res, model X = 1 : ')
disp(avg_t_1r)

t_vec1d = zeros(N,1);
for ii = 1:N
    z1 = rand(n1,1);
    tic;
    d(ii) = norm(D\(A1*z1)-z1)/norm(A1*z1);
    t_vec1d(ii) = toc;
end
avg_t_1d = mean(t_vec1d);
disp('Average time for diag, model X = 1 : ')
disp(avg_t_1d)

figure
plot(1:N,r,'r')
hold on
plot(1:N,d,'b')
legend('norm(T2\(T1\(A*z))-z)/norm(A*z)', 'norm(D\(A*z)-z)/norm(A*z)')
title('X = 1')
xlabel('Iteration (Trial)')
ylabel('Relative error')

% data 2
n2 = length(T1b);
[x2, t2, A2, ~] = heatsink_build_matrices(data2_exercise2);
D2 = spdiags(diag(A2), 0, n2, n2);
r2 = zeros(N,1);
d2 = zeros(N,1);

t_vec2r = zeros(N,1);
for ii = 1:N
    z2 = rand(n2,1);
    tic;
    r2(ii) = norm(T2b\(T1b\(A2*z2))-z2)/norm(A2*z2);
    t_vec2r(ii) = toc;
end
avg_t_2r = mean(t_vec2r);
disp('Average time for res, model X = 2 : ')
disp(avg_t_2r)

t_vec2d = zeros(N,1);
for ii = 1:N
    z2 = rand(n2,1);
    tic;
    d2(ii) = norm(D2\(A2*z2)-z2)/norm(A2*z2);
    t_vec2d(ii) = toc;
end
avg_t_2d = mean(t_vec2d);
disp('Average time for diag, model X = 2 : ')
disp(avg_t_2d)

figure
plot(1:N,r,'r')
hold on
plot(1:N,d,'b')
legend('norm(T2\(T1\(A*z))-z)/norm(A*z)', 'norm(D\(A*z)-z)/norm(A*z)')
title('X = 2')
xlabel('Iteration (Trial)')
ylabel('Relative error')

%% 6 e) 
clc, clear, close all

data1 = load('heatsink_exercise_ref1.mat');
load('heatsink_precond_ref1.mat');
[x1, t1, A1, rhs] = heatsink_build_matrices(data1);

M_diag = @(r) r ./ diag(A1);
M_T = @(r) T2\(T1\r);

b = rhs;
m = 100;


[x_tildaT, rel_resT] = GMRES_leftpc(A1,b,m,M_T);
[x_tildaD, rel_resD] = GMRES_leftpc(A1,b,m,M_diag);

[x_pcgD, rvD, ~] = pcg_v1(A1, b, m, M_diag, 1);
[x_pcgT, rvT, ~] = pcg_v1(A1, b, m, M_T, 1);


[~, ~, tvecD] = pcg_v1(A1, b, m, M_diag, 0);
[~, ~, tvecT] = pcg_v1(A1, b, m, M_T, 0);

%t_vec_GMRES_MT = zeros(m,1);
%t0 = tic;
%for k = 1:m
%    [~,~] = GMRES_leftpc(A1,b,m,M_T);
%    t_vec_GMRES_MT(k) = toc(t0);
%end

%t_vec_GMRES_D = zeros(m,1);
%t0 = tic;
%for k = 1:m
%    [~,~] = GMRES_leftpc(A1,b,m,M_diag);
%    t_vec_GMRES_D(k) = toc(t0);
%end
tic;
x = A1\rhs;
t_backslash = toc;
res_back = norm(A1*x - rhs)/norm(rhs);
figure
semilogy(1:m, rel_resT)
hold on
semilogy(1:m, rel_resD)
hold on
semilogy(1:m, rvD)
hold on
semilogy(1:m, rvT)
grid on
%ylim([1e-18 1e3])
xlim([0 m])
xlabel('Iteration', 'Interpreter', 'latex')
ylabel('$\frac{\|Az - b\|_2}{\|b\|_2}$', 'Interpreter', 'latex')
legend('GMRES with precond M','GMRES with precond Diag','PCG with precond D','PCG with precond M')
hold off


figure
semilogy(tvecD,rvD)
hold on
semilogy(tvecT,rvT)
hold on
xline(t_backslash, 'r--', 'DisplayName', 'Backslash time')
grid on
xlim([0 0.09])
ylim([1e-7 1e0])
xlabel('CPU time in seconds', 'Interpreter', 'latex')
ylabel('$\frac{\|Az - b\|_2}{\|b\|_2}$', 'Interpreter', 'latex')
legend('CG-D', 'CG-M','Backslash')
hold off









%%% functions %%%

function [z,rel_res] = GMRES_leftpc(A,b,m,Minv)
    r0 = Minv(b);
    beta = norm(r0);

    rel_res = zeros(m, 1);

    [Q,H] = arnoldi_6e(A,r0,m,Minv);
    e1 = zeros(m + 1, 1);
    e1(1) = 1;

    for k = 1:m
        z_star_k = H(1:k + 1, 1:k) \ (beta * e1(1:k + 1));
        z_k = Q(:, 1:k) * z_star_k;
        
        rel_res(k) = norm(A*z_k - b) / norm(b);
    end
    z_star = H \ (beta * e1);
    z = Q(:,1:m)*z_star;
end

function [x, t, A, rhs] = heatsink_build_matrices(data) %FROM CODE GIVEN
    K = data.K;          % full stiffness
    F0 = data.F;          % base load vector (without Gaussian)
    coords = data.coords;     % N x 3
    tri = data.tri;        % M x 3
    baseNodes = data.baseNodes;  % indices of base nodes (y ≈ ymin)
    topNodes = data.topNodes;   % indices of top nodes  (y ≈ ymax)

    N = size(coords,1);
    Tinf = 293.15;   % K, ambient
    T0 = 363.15;   % K, peak at center (e.g. 90°C)
    sigma = 0.03;     % m, Gaussian width
    x_center = 0.0;     % m
    z_center = 0.0;     % m

    % Extract base node positio ns
    x_base = coords(baseNodes,1);
    z_base = coords(baseNodes,3);

    % Gaussian temperature profile at base
    dx2 = (x_base - x_center).^2;
    dz2 = (z_base - z_center).^2;
    Tbase = Tinf + (T0 - Tinf) .* exp(-(dx2 + dz2) / (2*sigma^2));

    % Initialize full temperature vector with ambient
    T = Tinf * ones(N,1);
    T(baseNodes) = Tbase;

    fixed = false(N,1);
    fixed(baseNodes) = true;
    free  = ~fixed;

    % Reduced system: A * Tfree = Ff - Kfc * Tfixed
    A = K(free,free);
    Kfc = K(free,fixed);
    Ff  = F0(free);

    rhs = Ff - Kfc * T(fixed);

    % Solve 
    tic;
    x= A\rhs;
    t = toc;
    % The vector x now represents the solution 

end
function [x, r_hist, res_norms, t_vec] = cg(A, b, m, compute_res)

    x = zeros(size(b));
    r = b;
    p = r;
    r_hist = zeros(m, 1);
    res_norms = zeros(m, 1);
    t_vec = zeros(m, 1);

    t0 = tic;

    for k = 1:m
        rr = r'*r;
        Ap = A*p;

        alpha = (rr) / (p'*Ap);
        x = x + (alpha*p);
        r = r - (alpha*Ap);
        beta = (r'*r)/(rr);
        p = r + beta*p;

        r_hist(k) = norm(r);
        
        if (compute_res)
            res_norm = norm(A*x - b);
            res_norms(k) = res_norm;
        end
        
        t_vec(k) = toc(t0);

        if norm(r) < 1e-15
            return;
        end
    end
end

function [x, r_hist, res_norms, t_vec] = cg_LSQ(A, b, m, compute_res)

    n = length(b);
    x = zeros(n, 1);
    r = b;
    p = r;
    r_hist = zeros(m, 1);
    res_norms = zeros(m, 1);
    x_vals = zeros(n, m);
    t_vec = zeros(m, 1);

    t0 = tic;

    for k = 1:m
        rr = r'*r;
        Ap = A*p;

        alpha = (rr) / (p'*Ap);
        x = x + (alpha*p);
        r = r - (alpha*Ap);
        beta = (r'*r)/(rr);
        p = r + beta*p;

        r_hist(k) = norm(r);
        x_vals(:, k) = x;
        
        % LSQ part
        if (compute_res)
            Ax = A*x_vals(:, 1:k);
            B = (Ax)'*(Ax);
            c = Ax'*b;
    
            z = B\c;
            Xz = x_vals(:, 1:k)*z;
            
            res_norm = norm(A*Xz - b);
            res_norms(k) = res_norm;
        end
        
        t_vec(k) = toc(t0);

        if norm(r) < 1e-15
            %return;
        end
    end
end

function [x, r_hist, res_norms, t_vec] = cgne(A, b, m, compute_res)

    x = zeros(size(b));
    r = A'*b; %%%
    p = r;
    r_hist = zeros(m, 1);
    res_norms = zeros(m, 1);
    t_vec = zeros(m, 1);

    t0 = tic;

    for k = 1:m
        rr = r'*r;
        Ap = A*p;

        alpha = (rr) / (Ap'*Ap); %%%
        x = x + (alpha*p);
        r = r - alpha*(A'*Ap); %%%
        beta = (r'*r)/(rr);
        p = r + beta*p;

        r_hist(k) = norm(r);
        
        if (compute_res)
            res_norm = norm(A*x - b);
            res_norms(k) = res_norm;
        end
        
        t_vec(k) = toc(t0);

        if norm(r) < 1e-15
            return;
        end
    end
end

function z = GMRES(A,b,m)
    beta = norm(b);
    [Q,H] = arnoldi(A,b,m);
    
    e1 = zeros(m + 1, 1);
    e1(1) = 1;
    z_star = H \ (beta * e1);
    z = Q(:,1:m)*z_star;
end