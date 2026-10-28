function out = nav_apply_error(state, dx)
% Right/local error injection. Error order: p,v,theta,ba,bg,wind,bh,clock,drift.
dx=dx(:);
assert(numel(dx)==21 && all(isfinite(dx)), 'libekf:ErrorState', 'Expected 21 finite errors.');
nav_rotation(state.q); % Validate before constructing candidate.
out=state;
out.p=state.p(:)+dx(1:3); out.v=state.v(:)+dx(4:6);
a=norm(dx(7:9));
if a<1e-8
    scale=0.5-a*a/48;
else
    scale=sin(a/2)/a;
end
dq=[cos(a/2); scale*dx(7:9)];
q=state.q(:); w=q(1); v=q(2:4);
out.q=[w*dq(1)-v.'*dq(2:4); w*dq(2:4)+dq(1)*v+cross(v,dq(2:4))];
out.q=out.q/norm(out.q);
out.ba=state.ba(:)+dx(10:12); out.bg=state.bg(:)+dx(13:15);
out.wind=state.wind(:)+dx(16:18);
out.bh=state.bh+dx(19); out.clock=state.clock+dx(20); out.drift=state.drift+dx(21);
end
