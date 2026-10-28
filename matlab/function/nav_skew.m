function S = nav_skew(v)
% Cross-product matrix: S*w = cross(v,w).
v = v(:);
assert(numel(v)==3 && all(isfinite(v)), 'libekf:Vector', 'Expected finite 3-vector.');
S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end
