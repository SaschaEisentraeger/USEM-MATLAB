function [p,q,r,C] = TimeIntgCoeff(M,a)
%

% M: order
% a: a weight factor. It equals to the asymptotic of the spectrum at high frequency 

L = M;
ii=0:L; 
p1=(factorial(M+L-ii))./(factorial(ii).*factorial(L-ii)); 
ii=0:M; 
q1=(factorial(M+L-ii)).*((-1).^ii)./(factorial(ii).*factorial(M-ii))...
    *factorial(M)/factorial(L); 
L = M-1;
ii=0:L; 
p2=(factorial(M+L-ii))./(factorial(ii).*factorial(L-ii)); 
ii=0:M; 
q2=(factorial(M+L-ii)).*((-1).^ii)./(factorial(ii).*factorial(M-ii))...
    *factorial(M)/factorial(L); 
p = a*p1 + (1-a)*[p2 zeros(M-L)];
q = a*q1 + (1-a)*q2;
r = roots(q(end:-1:1));
r(imag(r)<-1.d-6) = [];
[~, idx] = sort(imag(r));
r = r(idx);

tmp = p - q;
C = zeros(M+1,M);
C(1,:) = tmp(2:end);
f1 = 1;
f2 = 1;
for ii = 1:M
    f1 = -f1/2;
    f2 = -f2;
    tmp = ii*[C(ii,:) zeros(M-L)] + f1*(p-f2*q);
    C(ii+1,:) = tmp(2:end);
end

end