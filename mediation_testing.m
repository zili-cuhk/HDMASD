tic
format short
clc;clear;close;

%% === model parameters ===
n = 500; p = 10000;
N = 1; 
d = round(n/log(n));
% d = round(n/log(n)/3);

B = 500; % Bootstrap次数

rho = [0.25 0.5];
rho = rho(2);
mu = zeros(p,1); c = (1:p);
sigma = eye(p,p) + rho*(ones(p,p)-eye(p,p)); %% dependent structure

%% === 'control' parameter ===
% tau1 = 2.35;  tau2 = 0.50;  % rho=0.25  
tau1 = 3.25;  tau2 = 0.70; % rho=0.5   

%% === 'lambda' parameter 
% ===
lambda_mcp = 0.11*(0.001).^((1:100)/100);
lambda_lasso = 0.11*(0.001).^((1:100)/100);


%% === True Beta ===
Beta = [0.7;0.7;0.7;0.70;0.7;0;0;-3.5*rho;zeros(p-8,1)];
Alph = [0.35;0.35;0.30;0.25;0.30;0;0.30;0.30;zeros(p-8,1)];
% Alph = [0.30;0.30;0.30;0.30;0.30;0;0.30;0.30;zeros(p-8,1)];

gamma = 0.5;
theta = [0.5;0.5];
xi = [0.3;0.3];

index = find(Beta.*Alph~=0);
AB = Beta.*Alph;  % true value of beta*alpha


for iter = 1:N
    iter

    rng(iter)   
    [T,X,Z,M,Delta] = simulation_data(n,Beta,Alph,theta,xi,gamma,mu,sigma,tau1,iter);
    Censorrate(iter) = 1-mean(Delta);
    ini_value(:,iter) = initial_beta(n,[X Z M],1e-5,Delta);  % initial value

    %% Screening mediators
    Scre_beta(:,iter) = NonMargScr_Mediator(n,ini_value(:,iter),[X Z M],1e-3,Delta,d);

    % ini_value(:,iter) = initial_beta(n,[X Z M],1e-5,Delta);  % initial value
    % [opt_beta(:,iter),opt_alpha(:,iter),cov_beta(:,iter),variance_alpha(:,iter)] =...
    %     MIC(X,Z,M(:,index),Delta);

    for b = 1:B
        % 有放回地随机抽取n个样本
        indices = randi(n, 1, n); % 生成随机索引（允许重复）
        X = X(indices,:);
        Z = Z(indices,:);
        Delta = Delta(indices);
        M(:,index) = M(indices,index);
        [a_boot(:,b),alpha_boot(:,b),cov_beta1(:,b),variance_alpha1(:,b)] =...
        MIC(X,Z,M(:,index),Delta);
    end

    % 清理无效样本
    valid_idx = all(~isnan(alpha_boot),2) & all(~isnan(beta_boot),2);
    fprintf('有效Bootstrap样本: %d/%d\n', sum(valid_idx), B);
    alpha_boot = alpha_boot(valid_idx,:);
    beta_boot = beta_boot(valid_idx,:);

    %% 结果分析
    % alpha估计
    alpha_mean = mean(alpha_boot);
    alpha_se = std(alpha_boot);


    % beta估计
    beta_mean = mean(beta_boot);
    beta_se = std(beta_boot);

    opt_beat(:,iter) = mean(opt_beta1,2);

    %%non-parametric bootstrap estimate
    % data = [X,Z,M(:,index),Delta];
    % statistic = @(d) MIC(d(:,1), d(:,2), d(:,3), d(:,4));
    % bootstat = bootstrp(B, statistic, data);

end

mean(Censorrate)  % censoring rate

time = toc   % computing time


% %% Assessment criteria 
% beta_alpha = opt_beta.*opt_alpha;
% BIASE_beta_alpha = mean(beta_alpha,2)-AB(index) % Biase
  


%% ===================================================
%                 simulation_data()
% ============================================================
function  [T,X,Z,M,status] = simulation_data(n,Beta,Alph,theta,xi,gamma,mu,sigma,tau,iter)
%% Generating the survival data
p = length(Beta);
% d = round(n/log(n)/1.5);
% rng(iter);

% X = normrnd(0,1,[n,1]);
% Z1 = normrnd(0,2,[n,1]);
% Z2 = normrnd(0,2,[n,1]);

X = binornd(1,0.6,[n,1]);
Z1 = binornd(1,0.3,[n,1]);
Z2 = unifrnd(0,1,[n,1]);
Z = [Z1 Z2];

M = X.*repmat(Alph',n,1) + repmat(Z*xi,1,p) + mvnrnd(mu,sigma,n); % no intercept term

% M_1(:,1) = M(:,1);
% M_1(:,2:p) = M(:,2:p) + 0.2*M(:,1:p-1);
% M = M_1;

D = unifrnd(0,1,n,1);
Death_time = -log(1-D)./exp(gamma*X+M*Beta+Z*theta);  % % death time
C = unifrnd(0,tau,n,1);    % % censored time

status = (Death_time <= C);
T = min(Death_time,C);     % % survial time
[T,I] = sort(T,'descend'); % % sorting the time
% Y = bsxfun(@ge,T,T');    % % at risk process
X = X(I,:);
Z = Z(I,:);
M = M(I,:);
status = status(I);

end


%% ===================================================
%                 initial_beta()
% ============================================================
function initial_beta = initial_beta(n,Z,r,status)
[~,p] = size(Z);
beta=zeros(p,1);
k = 1; err = 0; 
% tk = 104;
tk = 240; 

while k<=1000&&err==0
    %k
    L_prime = -(sum(status.*(Z-cumsum((exp(Z*beta).*Z))./ cumsum(exp(Z*beta)))))'/n;
    beta1 =  beta - L_prime/tk;
    %beta1 =  beta - L_prime/(n*tk);
    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end

initial_beta = beta;

end



%% ===================================================
%                 NonMargScr_Mediator()
% ============================================================
function opt_beta = NonMargScr_Mediator(n,ini_beta,Z,r,status,d)
[~,p] = size(Z);
beta = ini_beta;
% beta = zeros(p,1);

k=1;err=0; tk = 104;
while k<=1000 && err==0
    %k
    L_prime = -(sum(status.*(Z-cumsum((exp(Z*beta).*Z))./ cumsum(exp(Z*beta)))))'/n;
    beta_tilde = beta - L_prime/tk;

    b1 = beta_tilde(1:3);
    b2 = sort(abs(beta_tilde(4:p)),'descend');
    b2 = beta_tilde(4:p).*(abs(beta_tilde(4:p)) >= b2(d));
    beta1 = [b1;b2];
    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end

opt_beta = beta1(4:p);

end


%% ===================================================
%                 ista_MIC()
% ============================================================
function [opt_beta,opt_alpha,cov_beta,variance_alpha] = MIC(X,Z,M,status)

P = [X,Z,M];
XZ = [X,Z];
[n,p1] = size(XZ);
[~,p2] = size(M);
p = p1+p2;
Q = zeros(p1+p2,1);

n0 = sum(status);
lambda0 = log(n0);
r = 1e-5;
%% === 'lambda' parameter ===
% a = 30:5:70;    % % tunning parameter for MIC-penalized variable slection method
a = 120;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  initial_beta
k1 = 1; err = 0; %tk1 = 24;
tk1 = 8;

while k1<=1000 && err==0
    %k1
    L_prime = -(sum(status.*(P-cumsum((exp(P*Q).*P))./ cumsum(exp(P*Q)))))'/n;
    Q1 =  Q-L_prime/tk1;
    w = Q1-Q;
    err = norm(w,2)^2 <= r*norm(Q,2)^2;
    Q = Q1;
    k1 = k1+1;

end

%%
Gamma_theta = Q(1:p1);
beta0 = Q((p1+1):p);
beta = beta0;

tk2 = 4; tk3 = 4;
k3 = 1; err2 = 0;
while k3<=1000 && err2==0

    %% calculate =%-%= beta
    k2 = 1; err1 = 0;
    while k2<=1000 && err1==0
        %k2

        % beta = beta.*(abs(beta)>1e-3);
        W1 = diag( (4*a*lambda0*exp(2*a*beta.^2))./((exp(2*a*beta.^2)+1).^2) );
        u = eye(p2) + 2/tk2*W1;

        if any(isnan(u))
            break;
        end

        L_prime1 = -(sum(status.*(M-cumsum((exp(M*beta+XZ*Gamma_theta).*M))./...
            cumsum(exp(M*beta+XZ*Gamma_theta)))))'/n;
        beta_tilde = beta - L_prime1/tk3;
        beta1 = u\beta_tilde;

        w1 = beta1-beta;
        err1 = norm(w1,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k2 = k2+1;

    end

    %% calculate =%-%= Gamma-theta
    L_prime_Gamma = -(sum(status.*(XZ-cumsum((exp(M*beta+XZ*Gamma_theta).*XZ))./...
        cumsum(exp(M*beta+XZ*Gamma_theta)))))'/n;
    Gamma_theta1 =  Gamma_theta - L_prime_Gamma/tk2;
    w2 = Gamma_theta1-Gamma_theta;
    err2 = norm(w2,2)^2 <= r*norm(Gamma_theta,2)^2;
    Gamma_theta = Gamma_theta1;
    k3 = k3+1;

end

beta = beta.*(abs(beta)>1e-3);
opt_beta = beta;


%%
[n,p] = size(M);
[~,q1] = size(X);
[~,q2] = size(Z);
Q = [X,Z];
[~,q] = size(Q);

for j = 1:p

    Y = M(:,j); % -unifrnd(0,1,[n,1]);
    beta = (Q'*Q+0.08*eye(q))\Q'*Y;  % by OLS formulaiton

    alpha = beta(1:q1);
    eta = beta((q1+1):(q1+q2));

    k = 1; err1 = 0; err2 = 0;  tk = 4; 

    while k<=1000 && (err1==0 || err2==0)
        %k

        L_prime1 = (X'*X*alpha-X'*Y+X'*(Z*eta))/n;
        alpha1 =  alpha - L_prime1/tk;

        w1 = alpha1-alpha;
        err1 = norm(w1,2)^2 <= r*norm(alpha,2)^2;
        alpha = alpha1;

        L_prime2 = (Z'*Z*eta-Z'*Y+Z'*X*alpha)/n;
        eta1 =  eta - L_prime2/tk;
        w2 = eta1-eta;
        err2 = norm(w2,2)^2 <= r*norm(eta,2)^2;
        eta = eta1;

        k = k+1;

    end

    %% variance for =%-%= alpha_eta
    Sigma_2 = 1/(n-q)*sum((Y - Q*[alpha; eta]).^2);
    covariance_alpha_eta = inv(Q'*Q)*Sigma_2;
    variance_alpha_eta = diag(covariance_alpha_eta);
    variance_alpha(j) = variance_alpha_eta(1);

    opt_alpha(j) = alpha;


    %% %% % * --- ASE:the average of estimated standard error; --- *% % %%
    beta = opt_beta;
    eta = XZ*Gamma_theta+M*beta;
    t = cumsum(exp(eta));  % % n*1
    v = cumsum( (exp(eta).*M) ); % % n*p

    for i=1:n
        Cel1(1,i) = {exp(eta(i))*M(i,:)'*M(i,:)};
        Cel2(1,i) = {v(i,:)'*v(i,:)};
    end

    f1 = cumsum(cat(3,Cel1{:}),3);   % % 二阶偏导的第一项
    f2 = cumsum(cat(3,Cel2{:}),3);   % % 二阶偏导的第二项

    for i=1:n
        Cel3(1,i) = {f1(:,:,i)};
        Cel4(1,i) = {f2(:,:,i)};
        Cel3(1,i) =  cellfun(@(x) (status(i)/t(i)).*x, Cel3(1, i),'UniformOutput',false);
        Cel4(1,i) =  cellfun(@(x) (status(i)/t(i)^2).*x, Cel4(1,i),'UniformOutput',false);
    end

    L_primeprime = -( sum(cat(3,Cel4{:}),3) + sum(cat(3,Cel3{:}),3) );
    cov_beta = diag( pinv( -L_primeprime + 1e-5) );  

end



end








