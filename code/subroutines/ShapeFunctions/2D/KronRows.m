function C = KronRows(A,B)
%KRONROWS Row-wise Kronecker product of two basis matrices.
%
% Column-wise Kronecker product: C(i*nB+j+1,:) = A(i+1,:).*B(j+1,:)

nA = size(A,1);
nB = size(B,1);
C = reshape(reshape(B,nB,1,[]).*reshape(A,1,nA,[]),nA*nB,[]);
end
