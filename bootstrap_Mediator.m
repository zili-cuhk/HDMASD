%% ===================================================
%                 bootstrap_Mediator()
% ============================================================
function [opt_beta,opt_alpha,cov_beta,variance_alpha,opt_Gamma_theta] = bootstrap_Mediator(T,X,Z,M,status)
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

opt_Gamma_theta = Gamma_theta;


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
    cov_beta = diag( inv( -L_primeprime + 1e-5) ); 



end