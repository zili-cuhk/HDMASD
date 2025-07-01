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
    k
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
