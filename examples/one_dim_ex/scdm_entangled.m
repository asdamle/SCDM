function [U, piv] = scdm_entangled(U,f,r)

[~, ~, piv] = qr((U*f)',0);
piv = piv(1:r);

T = (U(piv,:)*f)';
U = U * (T / sqrtm(T'*T));
