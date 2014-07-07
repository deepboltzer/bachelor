function [mu] = rpamt(beta, epsilon, tau, mu, sigma, r,n,k,pc)
% Ranking probability algorithm for multiple team games
% Usage: 
% [v] = rpamt(4.166666666666667,4.943921970341054,1,[25.00,25.00,25.0000,25.000],[8.333,8.333,8.333,8.333,8.333],[2,1],[2,2],2,4); disp(v)
% disp([v]);
% Shows updated means for a two team setup with two players in each team
% for given starting means and varainces. 
%
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

for j = 1 : pc,
    diag(j,j) = sigma(j)^2;
end

% Set up mean and covariance to approximate truncated gaussian
disp(transpose(A));
u = mu*A; C = transpose(A)*(beta^2*eye(pc)+diag)*A;
disp(u); disp(1/C);
pex = 0; psig = 0; 
%
% Compute value of constant function of a truncated Gaussian and mean z of
% a truncated Gaussian with input parameters [u,C]
%
[p,z,Z,e] = qsimvn( 5000, u, C, ai, bi );
v = A * (1/C) * (1/p); v = v * (z-u);
%
% Update players means.
%
lower = 1;
upper = 0;
disp(k);
for j = 1 : k,
    upper = upper + n(j);
    for i = lower : upper, 
        mu(i) = mu(i) + 10*sigma(i)^2*v(i);
    end
    lower = upper+1;
end
return