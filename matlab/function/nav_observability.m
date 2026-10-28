function out=nav_observability(H,R,scales,targets,relative_tolerance)
% Локальная определимость заданных величин из измерений без априорной ковариации.
% H и targets относятся к одним ошибкам состояния; scales задаёт их физические масштабы.
% R включает зависимости всех строк наблюдений. Это не проверка целостности.
[m,n]=size(H); scales=scales(:);
assert(m>0&&n>0&&isreal(H)&&all(isfinite(H(:)))&&numel(scales)==n&& ...
    isreal(scales)&&all(isfinite(scales)&scales>0)&&size(targets,2)==n&& ...
    isreal(targets)&&all(isfinite(targets(:)))&& ...
    isscalar(relative_tolerance)&&isfinite(relative_tolerance)&& ...
    relative_tolerance>0&&relative_tolerance<1, ...
    'libekf:Observability','Invalid observation geometry, scales or tolerance.');
assert(isequal(size(R),[m m])&&isreal(R)&&all(isfinite(R(:)))&& ...
    norm(R-R.','fro')<=1e-12*max(norm(R,'fro'),realmin)&&all(diag(R)>0), ...
    'libekf:Covariance','Expected full measurement covariance.');
d=sqrt(diag(R)); normalized=(R./d)./d.';
[L,flag]=chol((normalized+normalized.')/2,'lower');
assert(flag==0,'libekf:Covariance','Measurement covariance must be positive definite.');
A=L\((H.*scales.')./d); C=targets.*scales.';
assert(all(isfinite(A(:)))&&all(isfinite(C(:))), ...
    'libekf:Observability','Scaling overflow.');
[~,S,V]=svd(A,'econ'); singular=diag(S);
tol=max(relative_tolerance,100*eps*max(m,n));
threshold=tol*max(singular); retained=singular>threshold;
V=V(:,retained); values=singular(retained);
projection=C*V; unresolved=C-projection*V.';
lengths=sqrt(sum(C.^2,2)); fraction=sqrt(sum(unresolved.^2,2));
nonzero=lengths>0; fraction(nonzero)=fraction(nonzero)./lengths(nonzero);
identified=fraction<=tol;
variance=sum((projection./values.').^2,2);
variance(~identified)=Inf;
out=struct('rank',sum(retained),'identified',identified,'unresolved_fraction',fraction, ...
    'variance',variance,'singular_values',singular,'threshold',threshold);
end
