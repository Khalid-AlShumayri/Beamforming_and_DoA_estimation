function [E,L] = nnzmat(Q,sig)
% This function consider the nonzero eigenvalues of a matrix and return
% the nonzero block and the corresponding columns of the matrix

L = diag(sig);
flag = find(L>0.0001);
L = L(flag);
L = diag(L);
E = Q(:,flag);