function out=nav_save_reference(filter)
% Сохранить текущее состояние как опору; это не новое измерение.
% Старая опора заменяется, история использованных исходных данных сохраняется.
P=filter.P;
assert(isequal(size(P),[21 21])&&all(isfinite(P(:)))&& ...
    norm(P-P.','fro')<=1e-12*max(norm(P,'fro'),realmin), ...
    'libekf:Covariance','Invalid current covariance.');
d=sqrt(diag(P)); assert(isreal(d)&&all(isfinite(d)&d>0),'libekf:Covariance','Invalid variances.');
scaled=(P./d)./d.'; [~,flag]=chol((scaled+scaled.')/2);
assert(flag==0,'libekf:Covariance','Current covariance must be positive definite.');
assert(isscalar(filter.state.t)&&isfinite(filter.state.t)&& ...
    ischar(filter.state.frame_id)&&~isempty(filter.state.frame_id), ...
    'libekf:State','Invalid current epoch or frame.');
nav_rotation(filter.state.q);
assert(iscell(filter.used)&&all(cellfun(@(x)ischar(x)&&~isempty(x),filter.used)), ...
    'libekf:Lineage','Invalid source history.');
out=filter; out.reference=filter.state; out.jointP=[P P;P P];
if isfield(filter,'reference_bank')
    cross=filter.reference_bank.cross_cov;
    assert(any(size(cross,1)==[21 42]),'libekf:Correlation','Invalid reference-set dependence.');
    cross=cross(end-20:end,:); out.reference_bank.cross_cov=[cross;cross];
end
end
