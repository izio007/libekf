function info=nav_innovation(P,obs)
% Dimensionless normalized innovation using Cholesky solves, never inv(S).
% Diagonal equilibration makes the result invariant to physical unit changes.
H=obs.H; R=obs.R; r=obs.residual(:);
assert(isequal(size(P),[size(H,2) size(H,2)]) && all(isfinite(P(:))), ...
    'libekf:Covariance','Invalid prior covariance.');
assert(norm(P-P.','fro')<=1e-12*max(norm(P,'fro'),realmin),'libekf:Covariance','Prior is not symmetric.');
S=H*P*H.'+R;
assert(all(isfinite(S(:))) && all(diag(S)>0),'libekf:Innovation','Nonfinite/nonpositive innovation covariance.');
d=sqrt(diag(S));
% Sequential division avoids forming d*d', which may overflow/underflow.
A=(S./d)./d.';
assert(norm(A-A.','fro')<1e-10,'libekf:Innovation','Innovation covariance lost symmetry.');
[L,flag]=chol((A+A.')/2,'lower');
assert(flag==0,'libekf:Innovation','Innovation covariance is not positive definite.');
u=L\(r./d);
assert(all(isfinite(u)),'libekf:Innovation','Whitened innovation is not finite.');
info=struct('S',S,'scale',d,'factor',L,'whitened',u,'nis',u.'*u);
end
