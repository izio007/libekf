function out=nav_observability_sequence(history,scales,targets,tolerance)
% Определимость текущих величин по временной последовательности без априорной P.
% x(k+1)=T{k}*x(k)+w(k); y(k)=H{k}*x(k)+n(k).
% Размер состояния может меняться; пустое наблюдение не добавляет информации.
if isstruct(history)&&isscalar(history)&&isfield(history,'steps'), history=from_steps(history); end
assert(isstruct(history)&&isscalar(history)&&all(isfield(history,{'H','transitions','times', ...
    'process_covariance','measurement_covariance','process_measurement_covariance','frame_id'})), ...
    'libekf:Observability','Incomplete temporal observation model.');
H=history.H; transitions=history.transitions; times=history.times(:);
assert(iscell(H)&&isvector(H)&&~isempty(H)&&iscell(transitions)&&numel(transitions)==numel(H)-1&& ...
    isnumeric(times)&&isreal(times)&&numel(times)==numel(H)&&all(isfinite(times))&&all(diff(times)>=0)&& ...
    ischar(history.frame_id)&&isrow(history.frame_id)&&~isempty(history.frame_id), ...
    'libekf:Observability','Invalid history epochs, transitions or frame.');
assert(all(cellfun(@(h)isnumeric(h)&&isreal(h)&&ismatrix(h)&&all(isfinite(h(:))),H)), ...
    'libekf:Observability','Invalid observation Jacobian in history.');
dimensions=cellfun(@(h)size(h,2),H); n=dimensions(1); count=numel(H);
rows=cellfun(@(h)size(h,1),H); m=sum(rows); nw=sum(dimensions(2:end));
assert(all(dimensions>0)&&m>0&&size(targets,2)==dimensions(end), ...
    'libekf:Observability','Invalid state or observation dimensions.');
Q=history.process_covariance; R=history.measurement_covariance; U=history.process_measurement_covariance;
assert(isequal(size(Q),[nw nw])&&isequal(size(R),[m m])&&isequal(size(U),[nw m]), ...
    'libekf:Covariance','Invalid full temporal noise dimensions.');
joint=[Q U;U.' R];
assert(isreal(joint)&&all(isfinite(joint(:)))&& ...
    norm(joint-joint.','fro')<=1e-12*max(norm(joint,'fro'),realmin)&&all(diag(joint)>=0), ...
    'libekf:Covariance','Invalid process/measurement covariance.');
d=sqrt(diag(joint)); zero=d==0;
assert(~any(any(joint(zero,:)~=0)),'libekf:Covariance','Impossible zero-variance dependence.');
d(zero)=1; normalized=(joint./d)./d.'; eigenvalues=eig((normalized+normalized.')/2);
assert(min(eigenvalues)>=-100*eps*(nw+m)*max(1,max(abs(eigenvalues))), ...
    'libekf:Covariance','Indefinite process/measurement covariance.');
Phi=eye(n); B=zeros(n,nw); O=zeros(m,n); F=zeros(m,nw); first=1; noise_first=1;
for k=1:count
    hk=H{k};
    if k>1
        T=transitions{k-1};
        assert(isnumeric(T)&&isreal(T)&&isequal(size(T),[dimensions(k) dimensions(k-1)])&&all(isfinite(T(:))), ...
            'libekf:Observability','Invalid error transition.');
        Phi=T*Phi; B=T*B;
        ix=noise_first:noise_first+dimensions(k)-1;
        B(:,ix)=eye(dimensions(k)); noise_first=noise_first+dimensions(k);
    end
    ix=first:first+rows(k)-1; first=first+rows(k);
    O(ix,:)=hk*Phi; F(ix,:)=hk*B;
end
effective=F*Q*F.'+R+F*U+U.'*F.'; effective=(effective+effective.')/2;
target_process=targets*B;
noise=struct('covariance',target_process*Q*target_process.', ...
    'cross_cov',target_process*(Q*F.'+U));
out=nav_observability(O,effective,scales,targets*Phi,tolerance,noise);
out.t=times(end); out.frame_id=history.frame_id;
end


function h=from_steps(input)
% Независимость шумов всей последовательности должна быть задана явно.
assert(all(isfield(input,{'initial_times','frame_id','noise_independent'}))&& ...
    isequal(input.noise_independent,true),'libekf:Correlation', ...
    'Step compilation requires explicitly independent process and measurement noises.');
epochs=input.initial_times(:);
assert(isnumeric(epochs)&&isreal(epochs)&&any(numel(epochs)==[1 2])&& ...
    all(isfinite(epochs))&&all(diff(epochs)>=0)&&iscell(input.steps)&& ...
    isvector(input.steps)&&~isempty(input.steps)&&~any(isfield(input, ...
    {'H','transitions','times','process_covariance','measurement_covariance','process_measurement_covariance'})), ...
    'libekf:History','Invalid initial epochs or ambiguous history representation.');
dimension=21*numel(epochs);
h=struct('H',{{zeros(0,dimension)}},'transitions',{{}},'times',epochs(end), ...
    'frame_id',input.frame_id,'process_covariance',zeros(0), ...
    'measurement_covariance',zeros(0));
for k=1:numel(input.steps)
    step=input.steps{k};
    assert(isstruct(step)&&isscalar(step)&&all(isfield(step,{'kind','data'})), ...
        'libekf:History','Incomplete accepted step.');
    d=step.data; T=[]; Q=[];
    if strcmp(step.kind,'propagation')
        assert(all(isfield(d,{'transition','process_covariance','t_start','t_end','frame_id'}))&& ...
            strcmp(d.frame_id,input.frame_id)&&isscalar(d.t_start)&&isscalar(d.t_end)&& ...
            d.t_start==epochs(end)&&isfinite(d.t_end)&&d.t_end>d.t_start, ...
            'libekf:History','Gap or frame mismatch in accepted propagation.');
        T=d.transition; Q=d.process_covariance; epochs(end)=d.t_end;
    else
        assert(isfield(d,'accepted')&&isequal(d.accepted,true)&&isfield(step,'time')&& ...
            isscalar(step.time)&&step.time==epochs(end), ...
            'libekf:History','Only accepted actions at the current epoch can enter the history.');
        if strcmp(step.kind,'save_reference')
            assert(isfield(d,'reference_change'),'libekf:History','Missing reference change.');
            c=d.reference_change;
            assert(all(isfield(c,{'transition','input_times','output_times','frame_id'}))&& ...
                isequal(c.input_times(:),epochs)&&isequal(c.output_times(:),[epochs(end);epochs(end)])&& ...
                strcmp(c.frame_id,input.frame_id),'libekf:History','Invalid reference change epochs or axes.');
            T=c.transition; Q=zeros(42); epochs=c.output_times(:);
        elseif isfield(d,'observation')
            o=d.observation;
            assert(all(isfield(d,{'state_times','error_reset'}))&&isequal(d.state_times(:),epochs)&& ...
                all(isfield(o,{'H','R','frame_id','noise_model'}))&&strcmp(o.frame_id,input.frame_id), ...
                'libekf:History','Observation does not match the current history state.');
            assert(strcmp(o.noise_model,'independent')&&~isfield(o,'peer_covariance')&& ...
                (~isfield(o,'state_noise_cov')||all(o.state_noise_cov(:)==0))&& ...
                (~isfield(o,'state_input_cov')||all(o.state_input_cov(:)==0)), ...
                'libekf:Correlation','Conditional or shared-reference errors require a full temporal model.');
            assert(size(o.H,2)==dimension&&isempty(h.H{end}),'libekf:History','Invalid observation dimensions.');
            h.H{end}=o.H;
            h.measurement_covariance=blkdiag(h.measurement_covariance,o.R);
            T=d.error_reset; Q=zeros(dimension);
        elseif any(strcmp(step.kind,{'source_health','reference_health'}))
            continue; % Исправность управляет допуском, но не добавляет наблюдения состояния.
        else
            error('libekf:History','Unsupported accepted action.');
        end
    end
    next_dimension=21*numel(epochs);
    assert(isequal(size(T),[next_dimension dimension])&&isequal(size(Q),[next_dimension next_dimension]), ...
        'libekf:History','Invalid transition or noise dimension.');
    h.transitions{end+1}=T; h.H{end+1}=zeros(0,next_dimension); h.times(end+1)=epochs(end);
    h.process_covariance=blkdiag(h.process_covariance,Q); dimension=next_dimension;
end
if isfield(input,'current_only')
    assert(isequal(input.current_only,true)||isequal(input.current_only,false), ...
        'libekf:History','Current-state selection must be logical.');
    if input.current_only&&dimension==42
        h.transitions{end+1}=[zeros(21) eye(21)]; h.H{end+1}=zeros(0,21);
        h.times(end+1)=epochs(end); h.process_covariance=blkdiag(h.process_covariance,zeros(21));
    end
end
h.process_measurement_covariance=zeros(size(h.process_covariance,1),size(h.measurement_covariance,1));
end
