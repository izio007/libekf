function C = nav_rotation(q)
% Hamilton scalar-first unit quaternion, active rotation (body-to-local in state).
q = q(:);
assert(numel(q)==4 && all(isfinite(q)) && abs(norm(q)-1)<1e-10, ...
    'libekf:Quaternion', 'Expected unit Hamilton quaternion [w;x;y;z].');
w=q(1); v=q(2:4);
C=(w*w-v.'*v)*eye(3)+2*(v*v.')+2*w*nav_skew(v);
end
