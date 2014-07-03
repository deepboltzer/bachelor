function [If] = sgcc(l,d)
 If = 0; m = l+d-1;
 for l = d : m,
     while k ~= 0,
        p = 1; k1= ones(d);k(1) = m; k2 = zeros(d); k2(:) = m;
        product = 1;
        If = If + product;
        k = iess(k1,k2,2,d);
     end
 end
return

function [x,w] = clencur(nl)
% Clenshaw Curtis Quadrature rule. 
x = 1:(nl -1); w = 1:(nl -1); 
w(1) = 1./(nl*(nl - 2)); x(1) = -cos(0);
for i=1 : (nl)
    x(i) = (1./2.)*(1-cos(pi*i/(nl+1)));
    sum = 0;
    for j = 1 : ((nl+1)/2)
        sum = sum + (1./(2*j-1))*sin(pi*i*(2*j-1)/(nl+1)); 
    end
    w(i) = (2/(nl +1))*sin(pi*i/(nl+1))*sum;
end
return

function [ind] = iepe(p,i,n,d)
% Drop algorithm for the iterative enumeration of all product indices for
% product rule method
while (i~=0)
    i(p) = i(p) +1;
    if i(p) > n(p), 
        if p == d, 
            ind = 0; 
            return; 
        end, 
        i(p) = 1; 
        p = p+1; 
    else
        p = 1; 
        ind = i; 
        return; 
    end
end
return

function [value] = f(~) 
    value = 1; 
return

function[ind] = iess(k1,k2,p,d)
% Drop algorithm for the iterative enumeration of all product indices for 
% sparse grid method
while(k1 ~= 0)
    k1(p) = k1(p)+1;
    if k1(p) > k2(p),
        if p == d, ind = 0; return; end
        k1(p) = 1;
        p = p+1;
    else
        for j = 1 : p-1,
            k2(j) = k2(p) - k1(p) +1;
        end
        k1(1) = k2(1);
        p = 1;
        ind = k1;
        return;
    end
end
return
% Values of the nk array must have values of the form nk = 2^(k-1) + 1; 
function [xev] = pqr(d,nk)
% Computation of product quadrature rule
n = 1;
for k = 1 : d, 
    n = n*nk(k); 
end
n = 1; i = ones(d); xev = zeros(d,n); wev = 1;
If = 0; p = 1; dim = max(nk); x = zeros(dim,d); w = zeros(dim,d);
% Set entries of the x and w matrices, which hold integration points and
% weights for Clenshaw Curtis rule.
for k = 1 : d, 
    h1 = zeros(nk(k)); h2 = zeros(nk(k)); [h1,h2] = clencur(nk(k));
    for j = 1: nk(k), 
        x(j,k) = h1(j); 
        w(j,k) = h2(j); 
    end
end
ind = 1;
while (i ~= 0)
        wev = 1;
        for j=1 : d, 
            xev(j,ind) = x(i(j),j);  
            wev = wev*w(i(j),j);  
        end
        If = If + wev*f(xev(:,ind));
        [i] = iepe(1,i,nk,d);
        ind = ind+1;
end
return
