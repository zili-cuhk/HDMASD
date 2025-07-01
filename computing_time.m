tic
format short
clc;clear;close;

%% === model parameters ===
n = 500; p = 20000;
N = 10; 
d = round(n/log(n));
% d = round(n/log(n)/3);

rho = [0.25 0.5 0.75];
rho = rho(2);
mu = zeros(p,1); c = (1:p);
% ama = bsxfun(@minus,c,c');
% sigma = rho.^(abs(ama));
sigma = eye(p,p) + rho*(ones(p,p)-eye(p,p)); %% dependent structure

%% === 'control' parameter ===
% tau1 = 2.35;  tau2 = 0.50; % rho=0.25  
tau1 = 3.25;  tau2 = 0.70; % rho=0.5

% tau1 = 2.35;  tau2 = 0.50;  % rho=0.25  
% % tau1 = 3.75;  tau2 = 0.735; % rho=0.5   

%% === 'lambda' parameter 
% ===
lambda_mcp = 0.11*(0.001).^((1:100)/100);
lambda_lasso = 0.11*(0.001).^((1:100)/100);

%% === 'lambda' parameter ===
% a = 100:10:130;    % % tunning parameter for MIC-penalized variable slection method
a = 90;

%% === True Beta ===
%% Zhang et al.(2021) Bioinformatics
% Beta = [0.55;0.6;0.65;0.7;0.5;0.5;zeros(p-6,1)];
% Alph = [0.5;0.45;0.5;0.4;0.45;0.45;zeros(p-6,1)];


Beta = [0.7;0.7;0.7;0.70;0.7;0;0;-3.5*rho;zeros(p-8,1)];
Alph = [0.35;0.35;0.30;0.25;0.30;0;0.30;0.30;zeros(p-8,1)];
% Alph = [0.30;0.30;0.30;0.30;0.30;0;0.30;0.30;zeros(p-8,1)];

gamma = 0.5;
theta = [0.5;0.5];
xi = [0.3;0.3];

index = find(Beta.*Alph~=0);
AB = Beta.*Alph;  % true value of beta*alpha


% t1=cputime;
%% Simulation Study
for iter = 1:N
    iter
    
    rng(iter)   
    [T,X,Z,M,Delta] = simulation_data(n,Beta,Alph,theta,xi,gamma,mu,sigma,tau1,iter);
    Censorrate(iter) = 1-mean(Delta);
    % ini_value(:,iter) = initial_beta(n,[X Z M],1e-5,Delta);  % initial value

    %% Screening mediators
    Scre_beta(:,iter) = NonMargScr_Mediator(n,[X Z M],1e-3,Delta,d);
    % Scre_beta(:,iter) = NonMargScr_Mediator(n,ini_value(:,iter),[X Z M],1e-3,Delta,d);
    % Scre_beta(:,iter) = NonMargScr_Mediator2(n,ini_value(:,iter),X,Z,M,1e-3,Delta,d);

    % SIS_beta(:,iter) = SIS_Mediator(n,ini_value(:,iter),X,Z,M,1e-3,Delta,d); % Luo et al. (2020)
    % opt_BetaAlpha(:,iter) = MargScr_Mediator(n,ini_value(:,iter),X,Z,M,1e-3,Delta,d);  % marginal mediation screening (Zhang et al.(2021))
    % ISIS_beta(:,iter) = ISIS_Mediator(n,ini_value(:,iter),X,Z,M,1e-3,Delta,d,lambda_mcp);  % ISIS mediation screening


end

% t2 = cputime;
% runtime = t2-t1;
% aver_runtime = vpa(runtime/N,3)
time = toc;   % computing time
time/N


%% Assessment criteria 
pj = mean(Scre_beta(index,:)~=0,2);  % proportion of individual active predictor is selected
pa = mean(all(Scre_beta(index,:)~=0));  % proportion of all active predictors are selected
TP = sum(sum(Scre_beta(index,:)~=0))/N;   % trues positives(TP)

% pj2 = mean(opt_BetaAlpha(index,:)~=0,2);  % proportion of individual active predictor is selected
% pa2 = mean(all(opt_BetaAlpha(index,:)~=0));  % proportion of all active predictors are selected
% TP2 = sum(sum(opt_BetaAlpha(index,:)~=0))/N;   % trues positives(TP)
% 
% pj3 = mean(ISIS_beta(index,:)~=0,2);  % proportion of individual active predictor is selected
% pa3 = mean(all(ISIS_beta(index,:)~=0));  % proportion of all active predictors are selected
% TP3 = sum(sum(ISIS_beta(index,:)~=0))/N;   % trues positives(TP)
% 
% pj4 = mean(SIS_beta(index,:)~=0,2);  % proportion of individual active predictor is selected
% pa4 = mean(all(SIS_beta(index,:)~=0));  % proportion of all active predictors are selected
% TP4 = sum(sum(SIS_beta(index,:)~=0))/N;   % trues positives(TP)


mean(Censorrate)  % censoring rate
% 
% time = toc   % computing time
% 
% [pj' pa TP; pj2' pa2 TP2; pj3' pa3 TP3; pj4' pa4 TP4 ]

[pj' pa TP]


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
% function opt_beta = NonMargScr_Mediator(n,ini_beta,Z,r,status,d)
function opt_beta = NonMargScr_Mediator(n,Z,r,status,d)
[~,p] = size(Z);
% beta = ini_beta;
beta = zeros(p,1);

k = 1; err = 0; 
% tk = 104;
tk = 400; 

while k<=1000&&err==0
    %k
    L_prime = -(sum(status.*(Z-cumsum((exp(Z*beta).*Z))./ cumsum(exp(Z*beta)))))'/n;
    beta1 =  beta - L_prime/tk;
    %beta1 =  beta - L_prime/(n*tk);
    w = beta1-beta;
    err = norm(w,2)^2 <= 1e-3*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end


k=1;err=0; %tk = 104;
tk = 240;

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
%                 NonMargScr_Mediator()
% ============================================================
function opt_beta = NonMargScr_Mediator2(n,ini_value,X,Z,M,r,status,d)
[~,p] = size(M);
Q = [X,Z];
[~,q] = size(Q);
% ini_value = zeros(p+q,1);
beta = ini_value(4:(p+q));
gt = ini_value(1:q);

k=1; err1=0; err2=0; 
% tk = 64;
tk = 104;

while k<=1000 && (err1==0 || err2==0)
    %k

    L_prime2 = -(sum(status.*(Q-cumsum((exp(M*beta+Q*gt).*Q))./cumsum(exp(M*beta+Q*gt)))))'/n;
    gt1 = gt - L_prime2/tk;
    w1 = gt1-gt;
    err1 = norm(w1,2)^2 <= r*norm(gt,2)^2;
    gt = gt1;

    %%
    L_prime2 = -(sum(status.*(M-cumsum((exp(M*beta+Q*gt).*M))./cumsum(exp(M*beta+Q*gt)))))'/n;
    beta_tilde = beta - L_prime2/tk;
    tem = sort(abs(beta_tilde),'descend');
    beta1 = beta_tilde.*(abs(beta_tilde) >= tem(d));
    w2 = beta1-beta;
    err2 = norm(w2,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;

    k = k+1;

end

opt_beta = beta1;

end


%% ===================================================
%                 MargScr_Mediator()  Zhang et al.(2021)
% ============================================================
function opt_BetaAlpha = MargScr_Mediator(n,ini_beta,X,Z,M,r,status,d)
Q = [X,Z];
[~,p] = size(M);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = M(:,j); 
eta = regress(Y,Q);  % regression
hat_alpha(j) = eta(1);

end


tem = sort(abs(hat_beta.*hat_alpha),'descend');
opt_BetaAlpha  = (hat_beta.*hat_alpha).*(abs(hat_beta.*hat_alpha) >= tem(d));

end


%% ===================================================
%                 SIS_Mediator()  Luo et al.(2020)
% ============================================================
function opt_beta = SIS_Mediator(n,ini_beta,X,Z,M,r,status,d)
[~,p] = size(M);
Q = [X,Z];

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


% Y = M(:,j); 
% %     beta = (Q'*Q)\Q'*Y;
% alpha = zeros(size(Q,2),1);
% 
% k1 = 1; err1 = 0; tk1 = 4;
% 
% while k1<=1000 && err1==0
%     %k1
%     L_prime = (Q'*Q*alpha-Q'*Y)/n;
%     alpha1 =  alpha - L_prime/tk1;
%     w = alpha1-alpha;
%     err1 = norm(w,2)^2 <= r*norm(alpha,2)^2;
%     alpha = alpha1;
%     k1 = k1+1;
% end

% hat_alpha(j) = alpha(1);

end

tem = sort(abs(hat_beta),'descend');
opt_beta = hat_beta.*(abs(hat_beta) >= tem(d));  % fisrt screening

% tem = sort(abs(hat_alpha),'descend');
% opt_beta = hat_alpha.*(abs(hat_alpha) >= tem(d));  % fisrt screening

end



%% ===================================================
%                 ISIS_Mediator()  Zhang et al.(2021)
% ============================================================
function opt_beta = ISIS_Mediator(n,ini_beta,X,Z,M,r,status,d,lambda_mcp)
[~,p] = size(M);
%beta = zeros(p,1);
% Q = [X,Z];

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
opt_beta1(1:3) = [];
opt_beta(opt_S) = opt_beta1;


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



