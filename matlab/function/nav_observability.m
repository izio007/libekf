function out=nav_observability(H,R,scales,targets,relative_tolerance,target_noise)
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
[left,S,V]=svd(A,'econ'); singular=diag(S);
tol=max(relative_tolerance,100*eps*max(m,n));
threshold=tol*max(singular); retained=singular>threshold;
V=V(:,retained); values=singular(retained);
projection=C*V; unresolved=C-projection*V.';
lengths=sqrt(sum(C.^2,2)); fraction=sqrt(sum(unresolved.^2,2));
nonzero=lengths>0; fraction(nonzero)=fraction(nonzero)./lengths(nonzero);
identified=fraction<=tol;
variance=sum((projection./values.').^2,2);
if nargin>=6
    % Целевая величина может содержать процессную ошибку, общую с измерениями.
    q=size(targets,1);
    assert(isstruct(target_noise)&&isscalar(target_noise)&& ...
        all(isfield(target_noise,{'covariance','cross_cov'}))&& ...
        isequal(size(target_noise.covariance),[q q])&&isequal(size(target_noise.cross_cov),[q m]), ...
        'libekf:Covariance','Invalid target noise dimensions.');
    W=target_noise.cross_cov; J=[R W.';W target_noise.covariance];
    assert(isreal(J)&&all(isfinite(J(:)))&& ...
        norm(J-J.','fro')<=1e-12*max(norm(J,'fro'),realmin)&&all(diag(J)>=0), ...
        'libekf:Covariance','Invalid joint measurement/target noise.');
    ds=sqrt(diag(J)); zero=ds==0;
    assert(~any(any(J(zero,:)~=0)),'libekf:Covariance','Impossible zero-variance dependence.');
    ds(zero)=1; normalized=(J./ds)./ds.';
    [vectors,eigenvalues]=eig((normalized+normalized.')/2,'vector');
    assert(min(eigenvalues)>=-100*eps*(m+q)*max(1,max(abs(eigenvalues))), ...
        'libekf:Covariance','Indefinite measurement/target noise.');
    G=(((W./d.')/L.')/L)./d.';
    effective=(targets-G*H).*scales.';
    weights=G+(((effective*V)./values.')*left(:,retained).'/L)./d.';
    % Норма через корень совместной ковариации избегает вычитания близких дисперсий.
    errors=([weights -eye(q)].*ds.')*vectors;
    variance=sum(errors.^2.*max(eigenvalues.',0),2);
end

variance(~identified)=Inf;
out=struct('rank',sum(retained),'identified',identified,'unresolved_fraction',fraction, ...
    'variance',variance,'singular_values',singular,'threshold',threshold);
end
