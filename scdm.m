function [U, piv] = scdm(U,ortho)
% Computes the orthogonalized or non-orthogonal SCDM from

% Anil Damle, Lin Lin, Lexing Ying, 
% "Compressed representation of Kohn?Sham orbitals 
% via selected columns of the density matrix,"
% J. Chem. Theory Comput., 2015, 11 (4), pp 1463?1469.


% Input:
%
%  U:      N x n_e sub-unitary matrix where each column is a kohn-sham
%          orbital
%  ortho:  1 if the orthognal orbitals are desired, 0 if the non-orthogonal
%          ones are desired

% Output:
%
%  U:      N x n_e localized basis
%  piv:    1 x n_e vector of selected columns

[Q, ~, piv] = qr(U',0);
piv = piv(1:size(U,2));
if ortho
    U = U*Q;
else    
    U = U*(U(piv,:)');
end