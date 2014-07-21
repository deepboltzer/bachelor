function[v1,v2,v3] = gcfnew()
e = linspace(-6,6,600);
for i = 1 : 600
 v1(i) = w(e(i),-0.5,0.5);
end
for i = 1 : 600
 v2(i) = w(e(i),-1,1);
end
for i = 1 : 600
 v3(i) = w(e(i),-4,4);
end
return

function [c] = v(t,l,u)
 a = normpdf(l-t);
 b = normpdf(u-t);
 nominator = a - b;
 denominator = phi(u-t) - phi(l-t);
 c = nominator/denominator;
return

function [c] = w(t,l,u)
 c = v(t,l,u) * v(t,l,u);
 nominator = (u-t) * normpdf(u-t) - (l-t)*normpdf(l-t);
 denominator = phi(u-t) - phi(l-t);
 c = c + nominator/denominator;
return

function [e] = normpdf(x)
 e = (1./sqrt(2*pi)) * exp(- (x*x)/(2));
return

function p =   phi(z), p =  erfc( -z/sqrt(2) )/2;
