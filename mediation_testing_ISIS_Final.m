tic
format short
clc;clear;close;

%% === model parameters ===
n = 500; p = 10000;
N = 200; 
d = round(n/log(n));
% d = round(n/log(n)/3);

B = 500; %

rho = [0.25 0.5];
rho = rho(2);
mu = zeros(p,1); c = (1:p);
sigma = eye(p,p) + rho*(ones(p,p)-eye(p,p)); %% dependent structure

%% === 'control' parameter ===
% tau1 = 2.35;  tau2 = 0.50;  % rho=0.25  
tau1 = 3.25;  tau2 = 0.70; % rho=0.5   

%% === 'lambda' parameter 
lambda = 4;
lambda_n = lambda*sqrt(n)/log(n);


%% === 'lambda' parameter 
% ===
% lambda_mcp = 0.11*(0.001).^((1:100)/100);
% lambda_mcp = 0.08*(0.001).^((1:100)/100);
% lambda_mcp = 0.06*(0.001).^((1:100)/100);
lambda_mcp = 0.035*(0.001).^((1:100)/100);

%% === True Beta ===
Beta = [0.7;0.7;0.7;0.70;0.7;0;0;-3.5*rho;zeros(p-8,1)];
Alph = [0.35;0.35;0.30;0.25;0.30;0;0.30;0.30;zeros(p-8,1)];
% Alph = [0.30;0.30;0.30;0.30;0.30;0;0.30;0.30;zeros(p-8,1)];

gamma = 0.5;
theta = [0.5;0.5];
xi = [0.3;0.3];

index = find(Beta.*Alph~=0);
AB = Beta.*Alph;  % true value of beta*alpha

beta_alpha = zeros(d,N);
opt_beta = zeros(d,N);
opt_alpha = zeros(d,N);
% cov_beta = zeros(d,N);
variance_alpha = zeros(d,N);
T_beta = zeros(d,N);
T_alpha = zeros(d,N);

beta_alpha_boot = zeros(d,N);
beta_boot = zeros(d,N);
alpha_boot = zeros(d,N);
cov_beta_boot = zeros(d,N);
variance_alpha_boot = zeros(d,N);
T_beta_boot = zeros(d,N);
T_alpha_boot = zeros(d,N);

U_star = zeros(d,N);
q_up = zeros(d,N);
q_low = zeros(d,N);


%%
for iter = 1:N
    iter

    rng(iter)   
    [T,X,Z,M,Delta] = simulation_data(n,Beta,Alph,theta,xi,gamma,mu,sigma,tau2,iter);
    Censorrate(iter) = 1-mean(Delta);
    ini_value(:,iter) = initial_beta(n,[X Z M],1e-5,Delta);  % initial value

    %% Screening mediators
    % Scre_beta(:,iter) = NonMargScr_Mediator(n,ini_value(:,iter),[X Z M],1e-3,Delta,d);
    [ISIS_beta(:,iter),opt_alpha(:,iter),cov_beta(:,iter),variance_alpha(:,iter)] = ISIS_Mediator(n,ini_value(:,iter),X,Z,M,1e-3,Delta,d,lambda_mcp);  % ISIS mediation screening


    % if any(ISIS_beta(:,iter)==0)
    %     continue;
    % else
        
    
    % index_Screening(:,iter) = find(Scre_beta(:,iter)~=0);
    index_Screening = find(ISIS_beta(:,iter)~=0);

    index_Selection = 1:length(index_Screening);

    opt_beta(index_Selection,iter) =ISIS_beta(index_Screening,iter);

    % 
    % % ini_value(:,iter) = initial_beta(n,[X Z M],1e-5,Delta);  % initial value
    % [opt_beta(:,iter),opt_alpha(:,iter),cov_beta(:,iter),variance_alpha(:,iter)] =...
    %     MIC(X,Z,M(:,index_Screening(:,iter)),Delta);

    
    % index_Selection = find(opt_beta(:,iter)~=0);
    % 
    index_test = index_Screening;

    %%
    % beta_alpha(index_Selection,iter) = opt_beta(index_Selection,iter).*opt_alpha(index_Selection,iter);
    % T_beta(index_Selection,iter) = sqrt(n)*opt_beta(index_Selection,iter)./sqrt(cov_beta(index_Selection,iter));
    % T_alpha(index_Selection,iter) = sqrt(n)*opt_alpha(index_Selection,iter)./sqrt(variance_alpha(index_Selection,iter));

    beta_alpha(:,iter) = opt_beta(:,iter).*opt_alpha(:,iter);
    T_beta(index_Selection,iter) = sqrt(n)*opt_beta(index_Selection,iter)./sqrt(cov_beta(index_Selection,iter));
    T_alpha(index_Selection,iter) = sqrt(n)*opt_alpha(index_Selection,iter)./sqrt(variance_alpha(index_Selection,iter));

%% Adaptive bootstrap testing
    for b = 1:B
        indices = randi(n, 1, n);
        X_b = X(indices,:);
        Z_b = Z(indices,:);
        T_b = T(indices);
        Delta_b = Delta(indices);
        M_b = M(indices,index_test);
        [beta_boot(index_Selection,b),alpha_boot(index_Selection,b),cov_beta_boot(index_Selection,b),variance_alpha_boot(index_Selection,b)] =...
        bootstrap_Mediator(T_b,X_b,Z_b,M_b,Delta_b);


        %%
        beta_alpha_boot(index_Selection,b) = beta_boot(index_Selection,b).*alpha_boot(index_Selection,b);
        T_beta_boot(index_Selection,b) = sqrt(n)*beta_boot(index_Selection,b)./sqrt(cov_beta_boot(index_Selection,b));
        T_alpha_boot(index_Selection,b) = sqrt(n)*alpha_boot(index_Selection,b)./sqrt(variance_alpha_boot(index_Selection,b));

        U_star(index_Selection,b) = (beta_alpha_boot(index_Selection,b)-beta_alpha(index_Selection,iter)).*...
            (1-max([abs(T_beta(index_Selection,iter)),abs(T_alpha(index_Selection,iter)),...
            abs(T_beta_boot(index_Selection,b)),abs(T_alpha_boot(index_Selection,b))],[],2)<= lambda_n);

    end


    q_low(index_Selection,iter) = quantile(U_star(index_Selection,:),0.05/2,2);
    q_up(index_Selection,iter) = quantile(U_star(index_Selection,:),1-0.05/2,2);

    index_final = ...
        find((beta_alpha(index_Selection,iter)<=q_low(index_Selection,iter))|...
        (beta_alpha(index_Selection,iter)>=q_up(index_Selection,iter)));


    %% fdr
    TP = sum(ismember(index_test(index_final),index));
    FP = sum(~ismember(index_test(index_final),index));

    fdr(iter) = FP/(TP+FP);

    %% power
    power_j(:,iter) = mean(ismember(index,index_test(index_final)),2);
 
    % end

end

mean(Censorrate)  % censoring rate


%% 
FDR = mean(fdr)
Power_individual = mean(power_j,2)
Power_all = mean(Power_individual)


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
tk = 104;
% tk = 240; 

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
% a = 90;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  initial_beta
k1 = 1; err = 0; tk1 = 24;
% tk1 = 8;

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

        exp_arg = 2 * a * beta.^2;
        max_exp = 700; 
        exp_arg = min(exp_arg, max_exp);
        exp_term = exp(exp_arg);

        numerator = 4*a*lambda0*exp_term + 1e-6;
        denominator = (exp_term + 1).^2 + 1e-6;

        mask_large = exp_arg > 700;
        numerator(mask_large) = 0;
        denominator(mask_large) = 1;

        mask_small = exp_arg < -700;
        numerator(mask_small) = 1e-6;
        denominator(mask_small) = 1 + 1e-6;

        W1 = diag(numerator./ denominator);

        % beta = beta.*(abs(beta)>1e-3);
        % W1 = diag( (4*a*lambda0*exp(2*a*beta.^2))./((exp(2*a*beta.^2)+1).^2) );
        u = eye(p2) + 2/tk2*W1;
        % 
        % if any(isnan(u))
        %     break;
        % end

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
    index_beta = opt_beta~=0;
    beta = opt_beta(index_beta);
    M1 = M(:,index_beta);
    eta = XZ*Gamma_theta+M1*beta;
    t = cumsum(exp(eta));  % % n*1
    v = cumsum( (exp(eta).*M1) ); % % n*p

    for i=1:n
        Cel1(1,i) = {exp(eta(i))*M1(i,:)'*M1(i,:)};
        Cel2(1,i) = {v(i,:)'*v(i,:)};
    end

    f1 = cumsum(cat(3,Cel1{:}),3);   
    f2 = cumsum(cat(3,Cel2{:}),3);  

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



%% ===================================================
%                 bootstrap_Mediator()
% ============================================================
function [opt_beta,opt_alpha,cov_beta,variance_alpha] = bootstrap_Mediator(T,X,Z,M,status)
%%

[~,I] = sort(T,'descend'); % % sorting the time
% Y = bsxfun(@ge,T,T');    % % at risk process
X = X(I,:);
Z = Z(I,:);
M = M(I,:);
status = status(I);

%%
P = [X,Z,M];
XZ = [X,Z];
[n,p1] = size(XZ);
[~,p2] = size(M);
p = p1+p2;
Q = zeros(p1+p2,1);

r = 1e-5;


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

tk2 = 4;  tk3 = 4;
k3 = 1; err2 = 0;
while k3<=1000 && err2==0

    %% calculate =%-%= beta
    k2 = 1; err1 = 0;
    while k2<=1000 && err1==0
        %k2

        L_prime1 = -(sum(status.*(M-cumsum((exp(M*beta+XZ*Gamma_theta).*M))./...
            cumsum(exp(M*beta+XZ*Gamma_theta)))))'/n;
        beta1 = beta - L_prime1/tk3;

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

end


    %% %% % * --- ASE:the average of estimated standard error; --- *% % %%
    index_beta = opt_beta~=0;
    beta = opt_beta(index_beta);
    M1 = M(:,index_beta);
    % beta = opt_beta;
    eta = XZ*Gamma_theta+M1*beta;
    t = cumsum(exp(eta));  % % n*1
    v = cumsum( (exp(eta).*M1) ); % % n*p

    for i=1:n
        Cel1(1,i) = {exp(eta(i))*M1(i,:)'*M1(i,:)};
        Cel2(1,i) = {v(i,:)'*v(i,:)};
    end

    f1 = cumsum(cat(3,Cel1{:}),3);   
    f2 = cumsum(cat(3,Cel2{:}),3);   

    for i=1:n
        Cel3(1,i) = {f1(:,:,i)};
        Cel4(1,i) = {f2(:,:,i)};
        Cel3(1,i) =  cellfun(@(x) (status(i)/t(i)).*x, Cel3(1, i),'UniformOutput',false);
        Cel4(1,i) =  cellfun(@(x) (status(i)/t(i)^2).*x, Cel4(1,i),'UniformOutput',false);
    end

    L_primeprime = -( sum(cat(3,Cel4{:}),3) + sum(cat(3,Cel3{:}),3) );
    cov_beta = diag( pinv( -L_primeprime + 1e-5) ); 



end




%% ===================================================
%                 ISIS_Mediator()  Zhang et al.(2021)
% ============================================================
function [opt_beta,opt_alpha,cov_beta,variance_alpha] = ISIS_Mediator(n,ini_beta,X,Z,M,r,status,d,lambda_mcp)
[~,p] = size(M);
%beta = zeros(p,1);
% Q = [X,Z];
XZ = [X,Z];

%% First marginal screening 
for j = 1:p
beta = ini_beta([1 2 3 j]);
k=1;err=0; tk = 4;
while k<=1000 && err==0
    %k
    ZZ = [X,Z,M(:,j)];
    L_prime = -(sum(status.*(ZZ-cumsum((exp(ZZ*beta).*ZZ))./ cumsum(exp(ZZ*beta)))))'/n;
    beta1 = beta - L_prime/tk;

    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end

hat_beta(j) = beta1(4);

end

tem = sort(abs(hat_beta),'descend');
S1 = find(abs(hat_beta) >= tem(d));  % fisrt screening 
[mcp_beta,~] = ista_MCP(n,[X Z M(:,S1)],1e-5,status,lambda_mcp); % first selection
beta2 = mcp_beta;
beta2(1:3) = [];
% SS1 = find(beta2);
SS1 = beta2~=0;
[SS2,IA,~] = intersect(S1,S1(SS1));
S2 = setdiff(1:p, SS2);  % unselected variables index

% M_beta = mcp_beta(horzcat([1 2 3], SS2+3));
M_beta = mcp_beta(horzcat([1 2 3], IA'+3));
ZZ = [X,Z,M(:,SS2)];
hat_beta1 = zeros(p,1);

%% Second screening stage
for l=1:length(S2)
    beta = ini_beta(S2(l));
    k=1;err=0; tk = 4;
while k<=1000 && err==0
    %k
    ZZ1 = M(:,S2(l));
    L_prime = -(sum(status.*(ZZ1-cumsum((exp(ZZ*M_beta+ZZ1*beta).*ZZ1))./ cumsum(exp(ZZ*M_beta+ZZ1*beta)))))'/n;
    beta1 = beta - L_prime/tk;

    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end

hat_beta1(S2(l)) = beta1;
  
end

tem = sort(abs(hat_beta1),'descend');
SS3 = find(abs(hat_beta1) >= tem(d));
opt_S = union(SS2,SS3);

opt_beta = zeros(p,1);

[opt_beta1,~] = ista_MCP(n,[X Z M(:,opt_S)],1e-5,status,lambda_mcp);

Gamma_theta = opt_beta1(1:3);
opt_beta1(1:3) = [];
opt_beta(opt_S) = opt_beta1;

% opt_beta = opt_beta';
% 
% size(opt_beta)


index_beta = find(opt_beta~=0);
M1 = M(:,index_beta);

opt_alpha = zeros(d,1);
cov_beta = zeros(d,1);
variance_alpha = zeros(d,1);



%%
[n,p1] = size(M1);
[~,q1] = size(X);
[~,q2] = size(Z);
Q = [X,Z];
[~,q] = size(Q);
% opt_alpha = zeros(p1,1);

for j = 1:p1

    Y = M1(:,j); % -unifrnd(0,1,[n,1]);
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
    variance_alpha(j,1) = variance_alpha_eta(1);

    opt_alpha(j,1) = alpha;

end


%% %% % * --- ASE:the average of estimated standard error; --- *% % %%

% index_beta = opt_beta~=0;
% M1 = M(:,index_beta);
beta = opt_beta(index_beta);

eta = XZ*Gamma_theta+M1*beta;
t = cumsum(exp(eta));  % % n*1
v = cumsum( (exp(eta).*M1) ); % % n*p

for i=1:n
    Cel1(1,i) = {exp(eta(i))*M1(i,:)'*M1(i,:)};
    Cel2(1,i) = {v(i,:)'*v(i,:)};
end

f1 = cumsum(cat(3,Cel1{:}),3);   
f2 = cumsum(cat(3,Cel2{:}),3);   

for i=1:n
    Cel3(1,i) = {f1(:,:,i)};
    Cel4(1,i) = {f2(:,:,i)};
    Cel3(1,i) =  cellfun(@(x) (status(i)/t(i)).*x, Cel3(1, i),'UniformOutput',false);
    Cel4(1,i) =  cellfun(@(x) (status(i)/t(i)^2).*x, Cel4(1,i),'UniformOutput',false);
end

L_primeprime = -( sum(cat(3,Cel4{:}),3) + sum(cat(3,Cel3{:}),3) );
cov_beta(1:length(index_beta)) = diag( pinv( -L_primeprime + 1e-5) );


end





%% ===================================================
%                 ista_MCP()
% ============================================================
function [opt_beta,opt_lambda] = ista_MCP(n,Z,r,status,lambda)
[~,p] = size(Z);
beta = zeros(p,1);
tv = 2.7;    % according to Cao et al.(2017)
opt_BIC=1e+10;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  initial_beta
k = 1;err = 0; tk = 8;
while k<=1000&&err==0
    %k
    L_prime = -(sum(status.*(Z-cumsum((exp(Z*beta).*Z))./ cumsum(exp(Z*beta)))))'/n;
    beta1 =  beta - L_prime/tk;
    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end

for i=1:length(lambda)
    k=1;err=0; tk = 4; % beta = ini_beta;
    while k<=1000 && err==0
%        k
        W1=diag(max(0,(tv*lambda(i)-abs(beta))/(tv*lambda(i)) )./(abs(beta)+1e-6));
        u = eye(p) + W1/tk;

        L_prime = -(sum(status.*(Z-cumsum((exp(Z*beta).*Z))./ cumsum(exp(Z*beta)))))'/n;

        beta_tilde = beta - L_prime/tk;
        beta1 = u\beta_tilde;
        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k = k+1;
    end

    beta(4:p) = beta(4:p).*(abs(beta(4:p))>2*1e-4);

    ell = -sum(status.*(Z*beta-log(cumsum(exp(Z*beta)))))/n;
    sel = beta~=0;
    BIC = ell+sum(sel)*log(n)/n;
    if BIC<=opt_BIC
        opt_BIC = BIC;
        opt_lambda = lambda(i);
        opt_beta = beta;
    end
end

end





