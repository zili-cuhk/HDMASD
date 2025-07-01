%% ===================================================
%                 NonMargScr_Mediator()
% ============================================================
function opt_beta = NonMargScr_Mediator(n,ini_beta,Z,r,status,d)
[~,p] = size(Z);
beta = ini_beta;
% beta = zeros(p,1);

Z(:,8:365313) = normalize(Z(:,8:365313)); 

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