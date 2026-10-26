%% 2 a)
clear all
B1 = [4 1 1; 1 3 0; 1 0 2];
N = 50;
%PM
ev1 = eig(B1);
analytic_eig1 = max(ev1);
x0 = [1 1 1]';
A1 = B1;
v = x0;
lvec = zeros(N,1);

for i = 1:N
    w = A1*v;
    v = w/norm(w);
    lambda = v' * A1*v;
    lvec(i) = lambda;
end
ratio = abs(ev1(2)/ev1(3));
disp('Ratio of lambda 2 / lambda 1, Conv.Rate')
disp(ratio)
errvec = abs(lvec - analytic_eig1);

%semilogy(errvec,'o')
%title('Semilogy plot of Power Method eigenvalue error')
%xlabel('Iteration')
%ylabel('Error')
%grid on
%hold on

ratio = ratio^2;
iter_count = 0:24;
first_error = errvec(1);
t_line = first_error * ratio.^iter_count;

%semilogy(iter_count, t_line,'--');
%legend('Error', 'Theoretical Rate')


errvec_2a_mini = errvec(1:22);
p_estimate_2a = log(errvec_2a_mini(3:end) ./ errvec_2a_mini(2:end-1)) ./ log(errvec_2a_mini(2:end-1) ./ errvec_2a_mini(1:end-2));
disp('Order of convergence')
disp(p_estimate_2a)

%%%%%%%%

B2 = [4 1 1; 1 3 0; 1 0 2];
% Inverse iteration, fixed shift
N = 50;
mu = 3;
x0 = [1 1 1]';
I2 = eye(size(B2));
ev2 = eig(B2);

v = x0;
A2 = B2;
lvec_2b = zeros(N,1);

for ii = 1:N
    w = (A2 - mu*I2)\v;
    v = w/norm(w);
    lambda = v' * A2*v;
    lvec_2b(ii) = lambda;
end

analytic_eig2 = ev2(2);
errvec_2b = abs(lvec_2b - analytic_eig2);

lambda_target = analytic_eig2;
lambda_near = ev2(1);
iter_count_2b = 0:23;
rate_2b = abs(lambda_target - mu)/abs(lambda_near - mu);
first_err_2b = errvec_2b(1);
t_line_2b = first_err_2b*rate_2b.^iter_count_2b;

semilogy(errvec,'o')
title('Comparison of Error, PM & II')
xlabel('Iteration')
ylabel('Error')
grid on
hold on

lambda_target = ev2(2);
distances = abs(ev2 - mu);
[~, sorted_idx] = sort(distances);
lambda_next = ev2(sorted_idx(2));  
convergence_rate_2b = abs((lambda_target - mu) / (lambda_next - mu));
disp('Convergence rate of Inv. It.')
disp(convergence_rate_2b)

semilogy(errvec_2b,'o')
%title('Semilogy plot of Inverse Iteration eigenvalue error')
%xlabel('Iteration')
%ylabel('Error')
grid on

%semilogy(iter_count, t_line,'-');
%semilogy(iter_count_2b,t_line_2b,'-');
legend('Error of PM','Error of Inverse iteration')

errvec_2b_mini = errvec_2b(1:11);
p_estimate_2b = log(errvec_2b_mini(3:end) ./ errvec_2b_mini(2:end-1)) ./ log(errvec_2b_mini(2:end-1) ./ errvec_2b_mini(1:end-2));
disp('Order of convergence')
disp(p_estimate_2b)
%% 2 c)
clear all
B = [4 1 1; 1 3 0; 1 0 2];
% Rayleigh quotient iteration
N = 10;
mu = 3;
x0 = [1 1 1]';
I = eye(size(B));
ev = eig(B);
analytic_eig_2c = ev(3);
v = x0;
A = B;
lvec_2c = zeros(N,1);
for jj = 1:N
    w = (A - mu*I)\v;
    v = w/norm(w);
    mu = v' * A*v; % lambda
    lvec_2c(jj) = mu;
end

errvec_2c = abs(lvec_2c - analytic_eig_2c);

semilogy(errvec_2c,'o')
title('Semilogy plot of Rayleigh quotient iteration eigenvalue error')
xlabel('Iteration')
ylabel('Error')
grid on

errvec_2c_mini = errvec_2c(1:4);
p_estimate_2c = log(errvec_2c_mini(3:end) ./ errvec_2c_mini(2:end-1)) ./ log(errvec_2c_mini(2:end-1) ./ errvec_2c_mini(1:end-2));
disp('Order of convergence')
disp(p_estimate_2c)

%% 2 d)
clear all
B = [4 1 1; 1.5 3 0; 1 0 2];
% Perturbed Matrix, RQI
N = 10;
mu = 3;
x0 = [1 1 1]';
I = eye(size(B));
ev = eig(B);
analytic_eig_2d = ev(1);
v = x0;
A = B;
lvec_2d = zeros(N,1);
for jj = 1:N
    w = (A - mu*I)\v;
    v = w/norm(w);
    mu = v' * A*v; % lambda
    lvec_2d(jj) = mu;
end

errvec_2d = abs(lvec_2d - analytic_eig_2d);

semilogy(errvec_2d,'o')
title('Semilogy plot of Rayleigh quotient iteration eigenvalue error')
xlabel('Iteration')
ylabel('Error')
grid on

errvec_2d_mini = errvec_2d(1:4);
p_estimate_2d = log(errvec_2d_mini(3:end) ./ errvec_2d_mini(2:end-1)) ./ log(errvec_2d_mini(2:end-1) ./ errvec_2d_mini(1:end-2));
disp('Order of convergence')
disp(p_estimate_2d)

%% 3 a) & b)
nn = 30;
rand('seed',0); A = gallery('wathen',nn,nn);
size(A);
b = randn(size(A,1), 1);
m = 80;
arn_x = zeros(m,1);
Km = zeros(length(A),m);
eig_Km = zeros(m,1);

% b)
Km(:,1) = b/norm(b);
tic
for jj = 2:m
    w = A * Km(:, jj-1);
    w = w/norm(w);

    Km(:, jj) = w;

    Km_curr = Km(:, 1:jj);
    G = Km_curr' * Km_curr;
    B = Km_curr' * A * Km_curr;

    eigvals = eig(B, G);
    eig_Km(jj) = max(real(eigvals));
end

for ii = 1:m
    %a)
    [Q,H]=arnoldi(A,b,ii);
    should_be_zero1 = norm(Q*H-A*Q(:,1:ii));
    should_be_zero2 = norm(Q'*Q-eye(ii+1));

    % b)
    H = H(1:ii,:);
    arnoldi_eig_approx = real(eig(H));
    arn_x(ii) = max(arnoldi_eig_approx);

end
toc

figure
plot(1:m,arn_x)
hold on
plot(1:m,eig_Km,'o')
legend('Eigenvalue approx from AM','Eigenvalue approx from (2)')


%%
load('Bwedge.mat')
%plot(real(B_eigvals),imag(B_eigvals),'o')

m = [2 4 8 10 20 30 40];
b = randn(size(B,1), 1);
diffvec = [];
figure; plot(real(B_eigvals),imag(B_eigvals),'o'); hold on
for ii = m
    [Q,H]=arnoldi(B,b,ii);
    H = H(1:ii,:);
    x = eig(H);
    plot(x,'*')
    diff = abs(x(1)-B_eigvals(1));
    diffvec = [diffvec, diff];
end
legend('Eigenvalues','m=2','m=4','m=8','m=10','m=20','m=30','m=40',location='northwest')
grid on
xlabel('Re')
ylabel('Im')
title('Ritz value convergence over increasing iterations')
%% plot with circles and such
clear all, close all

load('Bwedge')

lambda1 = B_eigvals(1);
lambda2 = B_eigvals(2);
lambda3 = B_eigvals(3);

plot(real(B_eigvals), imag(B_eigvals), '.')
hold on
grid on 
xlabel('Re')
ylabel('Im')
plot(real(lambda1), imag(lambda1), '*r')
hold on
plot(real(lambda2), imag(lambda2), '*g')
hold on
plot(real(lambda3), imag(lambda3), '*b') 

c_cluster = mean(B_eigvals(4:end));
abs(c_cluster-lambda1)
abs(c_cluster-lambda2)
abs(c_cluster-lambda3)
legend('All eigenvalues', 'fastest convergence', '2nd fastest convergence', '3rd fastest convergence')

thetav = 0:0.01:2*pi;

B_eigvals_lambda1 = B_eigvals(2:length(B_eigvals))  ;
c1 = mean([(-17.8202+i*2.10471), (1.33199+i*12.642), B_eigvals(3)])
rho1 = max(abs(c1-B_eigvals_lambda1))
plot(c1+rho1*exp(1i*thetav), '--r')
hold on
alpha1 = rho1/norm(B_eigvals(1)-c1)

B_eigvals_lambda2 = [lambda1; B_eigvals(3:end)];
c2 = complex(-25,-16)
rho2 = max(abs(c2-lambda1))
plot(c2+rho2*exp(1i*thetav), '--g')
alpha2 = rho2/norm(lambda2-c2)

B_eigvals_lambda3 = [B_eigvals(1:2); B_eigvals(4:end)];
c3 = complex(-22,10)
rho3 = max(abs(c3-lambda1))
plot(c3+rho3*exp(1i*thetav), '--b')
alpha3 = rho3/norm(lambda3-c3)

legend('All eigenvalues', 'eigenvalue 1 (fastest)', 'eigenvalue 2 (2nd fastest)', 'eigenvalue 3 (3rd fastest)', 'circle 1', 'circle 2', 'circle 3')
axis equal
title('Eigenvalues with outliers + circles for each eigenvalue outlier')

%%
clear all, close all
load('Bwedge.mat')
lambda_star_ish = -9.8 + 2*i;
[min_dist, idx] = min(abs(B_eigvals - lambda_star_ish));
lambda_star = B_eigvals(idx);

sigma1 = -10;
sigma2 = -7 + 2*i;
sigma3 = -9.8 + 1.5*i;

figure
transformed_eig1 = 1./(B_eigvals - sigma1);
transformed_lambda_star1 = 1/(lambda_star - sigma1);
plot(transformed_eig1,'o')
hold on
plot(transformed_lambda_star1, 'r*')
hold off
legend('Eigenvalues','\lambda_{*}')
axis([-2 2 -2 2])
title('\sigma_1')
xlabel('Re')
ylabel('Im')
grid on

figure
transformed_eig2 = 1./(B_eigvals - sigma2);
transformed_lambda_star2 = 1/(lambda_star - sigma2);
plot(transformed_eig2,'o')
hold on
plot(transformed_lambda_star2, 'r*')
hold off
legend('Eigenvalues','\lambda_{*}')
axis([-2 2 -2 2])
title('\sigma_2')
xlabel('Re')
ylabel('Im')
grid on

figure
transformed_eig3 = 1./(B_eigvals - sigma3);
transformed_lambda_star3 = 1/(lambda_star - sigma3);
plot(transformed_eig3,'o')
hold on
plot(transformed_lambda_star3, 'r*')
hold off
legend('Eigenvalues','\lambda_{*}')
axis([-2 2 -2 2])
title('\sigma_3')
xlabel('Re') 
ylabel('Im')
grid on

b = ones(size(B,1), 1);
mvec = [10 20 30];
err1vec = [];
err2vec = [];
err3vec = [];

for ii = mvec
    [Q1,H1] = shift_inv_arnoldi(B,b,ii,sigma1);
    H1 = H1(1:ii,:);
    ritz1 = eig(H1);
    error1 = min(abs(ritz1 - transformed_lambda_star1));
    err1vec = [err1vec,error1];

    [Q2,H2] = shift_inv_arnoldi(B,b,ii,sigma2);
    H2 = H2(1:ii,:);
    ritz2 = eig(H2);
    error2 = min(abs(ritz2 - transformed_lambda_star2));
    err2vec = [err2vec,error2];

    [Q3,H3] = shift_inv_arnoldi(B,b,ii,sigma3);
    H3 = H3(1:ii,:);
    ritz3 = eig(H3);
    error3 = min(abs(ritz3 - transformed_lambda_star3));
    err3vec = [err3vec,error3];
end

%%
clear all, close all
disp('Faces')
load ("olivetti_faces.mat")
imagesc(faces(:,:,67));
colormap('gray');
axis image;
m = 400;
n = 64^2;
X = reshape(faces,n,m);

x_mean = mean(X,2);
Xc = X - x_mean;

disp('b)')
disp('(1):')
tic
C = (1/(m-1))*(Xc*Xc');
toc

disp('(2):')
tic 
for ii = 1:100
    b = rand(length(C),1);
    C_b = C*b;
end
toc
disp('(3):')
tic 
for ii = 1:100
    b = rand(length(C),1);
    X_prim = Xc*(Xc'*b);
end
toc

% 8) PM 
N = 30;
ev8 = eig(C);
[analytic_eig8,idx8] = max(ev8);

x0 = ones(length(C),1);
v = x0;
lvec82 = [];
lvec83 = [];

% b1.2
disp('PM with method b.1-2')
tic 
for ii = 1:N
    w = C*v;
    v = w/norm(w);
    lambda = v' * C*v;
    lvec82(ii) = lambda;
end
toc
figure
errvec82 = abs(lvec82 - analytic_eig8)/norm(analytic_eig8);
semilogy(1:N,errvec82,'o')
title('Semilogy plot of Power Method eigenvalue error') %, 8 c), PM with C*b'
xlabel('Iteration')
ylabel('Error')
grid on
hold on
% b1.3
v = x0;
disp('PM with method b.1-3')

tic 
for ii = 1:N
    w = Xc*(Xc'*v);
    v = w/norm(w);
    lambda = v' * w /(m-1);
    lvec83(ii) = lambda;
end
toc



errvec83 = abs(lvec83 - analytic_eig8)/norm(analytic_eig8);
semilogy(1:N,errvec83,'o')
%title('Semilogy plot of Power Method eigenvalue error, 8 c), PM with Xc(Xc^T)')
%xlabel('Iteration')
%ylabel('Error')
%grid on
legend('Part b1-2','Part b3')

% find face

[eigvecs,diag_abc] = eig(C);
ev8 = eig(C);
[analytic_eig8,idx8] = max(ev8);

figure
eig_vec_max = eigvecs(:,idx8);
eigface = reshape(eig_vec_max,64,64);
imagesc(eigface);
colormap('gray');
axis image;

% d)

b = ones(size(Xc, 1), 1);
m = 30;
errors8d = zeros(5, m);
tic
for ii = 1:m
    [Q, H] = arnoldi_8(Xc, b, ii);
    H = H(1:ii,1:ii);
    

    [V, D] = eigs(H, 5, 'largestabs');
    eigenvalues = diag(D);

    num_eig = min(5,ii);
    for jj = 1:num_eig
        lambda = eigenvalues(jj);
        
        x = Q(:,1:ii) * V(:, jj);
        
        Cx = Xc * (Xc' * x);

        errors8d(jj, ii) = norm(Cx - lambda*x) / norm(Cx);
    end

end
toc
figure
semilogy(1:m, errors8d');
legend('\lambda_1', '\lambda_2', '\lambda_3', '\lambda_4', '\lambda_5');
xlabel('Iteration');
ylabel('||Cx - \lambda x|| / ||Cx||');
grid on;

b = ones(size(Xc, 1), 1);
m = 50; 

[Q, H] = arnoldi_8(Xc, b, m);
H_k = H(1:m, 1:m);

[V5, D5] = eigs(H_k, 5, 'largestabs');
eigenvalues_5 = diag(D5);

eigenfaces_5 = Q(:, 1:m) * V5; 

figure
for i = 1:5
    subplot(4, 5, i);
    eigenface_img = reshape(eigenfaces_5(:, i), 64, 64);
    imagesc(eigenface_img);
    colormap('gray');
    axis image;
    title(sprintf('Eigenface %d', i));
end

b = ones(size(Xc, 1), 1);
m = 60;  
[Q, H] = arnoldi_8(Xc, b, m);
H_k = H(1:m, 1:m);

[V, D] = eig(H_k);
eigenvalues_all = diag(D);

eigenvalues_first20 = eigenvalues_all(1:20);
V_first20 = V(:, 1:20);


eigenfaces_first20 = Q(:, 1:m) * V_first20;

figure;
for i = 1:20
    subplot(4, 5, i);
    eigenface_img = reshape(eigenfaces_first20(:, i), 64, 64);
    imagesc(eigenface_img);
    colormap('gray');
    axis image;
    title(sprintf('Eigenface %d', i));
end


% Functions were originally spererate files, here they are pasted in:
function [Q,H]=arnoldi(A,b,m)
% [Q,H]=arnoldi(A,b,m)
% A simple implementation of the Arnoldi method.
% The algorithm will return an Arnoldi "factorization":
%   Q*H(1:m+1,1:m)-A*Q(:,1:m)=0
% where Q is an orthogonal basis of the Krylov subspace
% and H a Hessenberg matrix.
%
% Example:
%  A=randn(100); b=randn(100,1);
%  m=10;
%  [Q,H]=arnoldi(A,b,m);
%  should_be_zero1=norm(Q*H-A*Q(:,1:m))
%  should_be_zero2=norm(Q'*Q-eye(m+1))
    n = length(b);
    Q = zeros(n,m+1);
    Q(:,1) = b/norm(b);
    
    t0 = 0;
    t1 = 0;
    
    for k = 1:m
        w = A*Q(:,k); % Matrix-vector product 
                    % with last element        
        %%% Orthogonalize w against columns of Q
        %replace this with a orthogonalization   % COMMENT OUT TO CHANGE
        %FUNCTIONS BELOW:::
        %[h,beta,worth] = hw1_MV_gs(Q,w,k);
        [h,beta,worth] = hw1_MV2_gs(Q,w,k);
        %[h,beta,worth] = hw1_MV3_gs(Q,w,k);
        %[h,beta,worth] = hw1_for_gs(Q,w,k);
        %[h,beta,worth] = hw1_mod_gs(Q,w,k);
        %%% Put Gram-Schmidt coefficients into H
        H(1:(k+1),k) = [h;beta]; 
        %%% normalize
        Q(:,k+1) = worth/beta; 
    end
end

function [h,beta,worth] = hw1_MV_gs(Q,w,k)
        h = Q(:,1:k)'*w;
        worth = w - Q(:,1:k)*h;
        beta = norm(worth);
end

function [h,beta,worth] = hw1_for_gs(Q,w,k)
    h = zeros(k,1);
    for j = 1:k
        h(j) = Q(:,j)'*w;
    end
    worth = w - Q(:,1:k)*h;
    beta = norm(worth);
end

function [h,beta,worth] = hw1_mod_gs(Q,w,k)
    h = zeros(k,1);
    worth = w;
    for j = 1:k
        h(j) = Q(:,j)'*worth;     
        worth = worth - h(j)*Q(:,j);  
    end
    beta = norm(worth);
end

function [h,beta,worth] = hw1_MV2_gs(Q,w,k)
        h = Q(:,1:k)'*w;
        worth = w - Q(:,1:k)*h;
        g = Q(:,1:k)'*worth;
        worth = worth - Q(:,1:k)*g;
        h = h + g;
        beta = norm(worth);
end

function [h,beta,worth] = hw1_MV3_gs(Q,w,k)
        h = Q(:,1:k)'*w;
        worth = w - Q(:,1:k)*h;
        g = Q(:,1:k)'*worth;
        worth = worth - Q(:,1:k)*g;
        f = Q(:,1:k)'*worth;
        worth = worth - Q(:,1:k)*f;
        h = h + g + f;
        beta = norm(worth);
end

function [Q,H]=shift_inv_arnoldi(A,b,m,mu)
    n = length(b);
    Q = zeros(n,m+1);
    Q(:,1) = b/norm(b);
    
    t0 = 0;
    t1 = 0;

    I = eye(size(A));
    [L, U] = lu(A - mu*I); 
    
    for k = 1:m
        %w = L\(U\Q(:,k)); wrong from lecture?
        w = U\(L\Q(:,k));
        %%% Orthogonalize
        %[h,beta,worth] = hw1_MV_gs(Q,w,k);
        [h,beta,worth] = hw1_MV2_gs(Q,w,k);
        %[h,beta,worth] = hw1_MV3_gs(Q,w,k);
        %[h,beta,worth] = hw1_for_gs(Q,w,k);
        %[h,beta,worth] = hw1_mod_gs(Q,w,k);
        %%% Put Gram-Schmidt coefficients into H
        H(1:(k+1),k) = [h;beta]; 
        %%% normalize
        Q(:,k+1) = worth/beta; 
    end
end

function [Q,H]=arnoldi_8(Xc,b,m)
% [Q,H]=arnoldi(A,b,m)
% A simple implementation of the Arnoldi method.
% The algorithm will return an Arnoldi "factorization":
%   Q*H(1:m+1,1:m)-A*Q(:,1:m)=0
% where Q is an orthogonal basis of the Krylov subspace
% and H a Hessenberg matrix.
%
% Example:
%  A=randn(100); b=randn(100,1);
%  m=10;
%  [Q,H]=arnoldi(A,b,m);
%  should_be_zero1=norm(Q*H-A*Q(:,1:m))
%  should_be_zero2=norm(Q'*Q-eye(m+1))
    n = length(b);
    Q = zeros(n,m+1);
    Q(:,1) = b/norm(b);
    
    t0 = 0;
    t1 = 0;
    
    for k = 1:m
        w = Xc*(Xc'*Q(:,k)); % Matrix-vector product 
                    % with last element        
        %%% Orthogonalize w against columns of Q
        %replace this with a orthogonalization 
        %[h,beta,worth] = hw1_MV_gs(Q,w,k);
        [h,beta,worth] = hw1_MV2_gs(Q,w,k);
        %[h,beta,worth] = hw1_MV3_gs(Q,w,k);
        %[h,beta,worth] = hw1_for_gs(Q,w,k);
        %[h,beta,worth] = hw1_mod_gs(Q,w,k);
        %%% Put Gram-Schmidt coefficients into H
        H(1:(k+1),k) = [h;beta]; 
        %%% normalize
        Q(:,k+1) = worth/beta; 
    end
end