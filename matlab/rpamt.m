function [p,pex,e] = rpamt(beta, epsilon, tau, mu, sigma, r,n,k,pc)
% Ranking probability algorithm for multiple team games
% Usage: [p] =
% rpamt(4.166666666666667,0.13,1,[25.000,25.000,25.000,25.000],[8.333,8.333,8.333,8.333],[1,2],[2,2],4);
% disp([p]);

lower = 1;
upper = 0;

for j = 1 : k,
    upper = upper + n(j);
    for i = lower : upper, 
        sigma(i) = sqrt(sigma(i)^2 + tau^2);
    end
    lower = upper+1;
end

% Set values of TeamPlayer matrix and integration limits
A = zeros(pc,k-1);
ai = zeros(k-1);
bi = zeros(k-1);

lower1 = 1; lower2 = 0; 
upper1 = 0; upper2 = 0;

for j = 1 : k-1,
    upper1 = upper1 + n(j);
    
    for i = lower1 : upper1, 
        A(i,j) = (2./(n(j)+n(j+1)));
    end
    
    lower1 = upper1+1;
    upper2 = upper1 + n(j+1);
    lower2 = upper1 +1; 
    
    for i = lower2 : upper2, 
        A(i,j) = ((-2.)/(n(j)+n(j+1)));
    end
    
    % Set integration limits
    if r(j) == r(j+1), 
        ai(j) = -epsilon;
        bi(j) = epsilon;
    else
        ai(j) = epsilon;
        bi(j) = inf;
    end
end

diag = zeros(pc,pc);

for k = 1 : pc,
    diag(k,k) = sigma(k)^2;
end

% Set up mean and covariance to approximate truncated gaussian
disp(transpose(A));
u = mu*A; C = transpose(A)*(beta^2*eye(pc)+diag)*A;
disp(u);
pex = 0; 
% Compute value of constant function of a truncated Gaussian
[p,pex,e] = qsimvn( 5000, u, C, ai, bi );
return
