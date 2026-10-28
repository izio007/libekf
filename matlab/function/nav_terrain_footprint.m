function hit=nav_terrain_footprint(origin,beam,halfAngle,vertices,faces,edgeMargin)
% Minimum Euclidean distance to a triangle mesh inside a circular cone.
% Continuous geometry, not a sampling of beam rays or just triangle vertices.
% Echo strength, sidelobes, DEM uncertainty/coverage and multiple returns are
% source/model obligations; this function alone establishes no integrity.
o=origin(:); b=beam(:);
assert(numel(o)==3 && numel(b)==3 && all(isfinite([o;b])) && abs(norm(b)-1)<1e-10 && ...
    isscalar(halfAngle) && isfinite(halfAngle) && halfAngle>0 && halfAngle<pi/2, ...
    'libekf:Terrain','Invalid circular footprint.');
assert(size(vertices,2)==3 && all(isfinite(vertices(:))) && size(faces,2)==3 && ~isempty(faces) && ...
    all(isfinite(faces(:))) && all(faces(:)==fix(faces(:))) && all(faces(:)>=1) && ...
    all(faces(:)<=size(vertices,1)) && isscalar(edgeMargin) && isfinite(edgeMargin) && edgeMargin>=0, ...
    'libekf:Terrain','Invalid local terrain mesh.');
c=cos(halfAngle); sn=sin(halfAngle);
hit=struct('range',Inf,'point',nan(3,1),'normal',nan(3,1),'smooth',false, ...
    'cone_active',false,'direction_cosine',NaN,'plane_distance',NaN);
for j=1:size(faces,1)
    triangle=vertices(faces(j,:),:).'; a=triangle(:,1); e1=triangle(:,2)-a; e2=triangle(:,3)-a;
    n=cross(e1,e2); area=norm(n);
    assert(area>0,'libekf:Terrain','Degenerate triangle in footprint.'); n=n/area;
    signedDistance=n.'*(a-o);
    if signedDistance~=0
        n=sign(signedDistance)*n; planeDistance=abs(signedDistance); k=max(-1,min(1,n.'*b));
        if k>=c
            u=n; cosine=1; coneActive=false;
        else
            transverse=sqrt(max(0,1-k*k));
            if transverse>1e-14
                u=c*b+sn*(n-k*b)/transverse; cosine=k*c+transverse*sn;
            else
                u=b; cosine=k;
            end
            coneActive=true;
        end
        if cosine>0
            distance=planeDistance/cosine; point=o+distance*u;
            gram=[e1.'*e1 e1.'*e2;e1.'*e2 e2.'*e2]; uv=gram\[e1.'*(point-a);e2.'*(point-a)];
            bary=[1-sum(uv);uv];
            if min(bary)>=-1e-12
                edgeDistance=min(bary.*area./[norm(e2-e1);norm(e2);norm(e1)]);
                hit=accept(hit,distance,point,n,edgeDistance>edgeMargin,coneActive,cosine,planeDistance);
            end
        end
    end
    % If the plane optimum misses the triangle, a constrained edge optimum
    % or vertex completes the convex triangle/cone minimum.
    for edge=[1 2;2 3;3 1].'
        aa=triangle(:,edge(1))-o; dd=triangle(:,edge(2))-triangle(:,edge(1));
        axial=b.'*aa; axialStep=b.'*dd;
        polynomial=[axialStep^2-c*c*(dd.'*dd), ...
            2*(axial*axialStep-c*c*(aa.'*dd)),axial^2-c*c*(aa.'*aa)];
        scale=max(abs(polynomial)); breaks=[0;1];
        if scale>0
            polynomial=polynomial/scale;
            if abs(polynomial(1))<1e-14
                if abs(polynomial(2))>1e-14, rr=-polynomial(3)/polynomial(2); else, rr=[]; end
            else, rr=roots(polynomial); end
            rr=real(rr(abs(imag(rr))<1e-10)); breaks=[breaks;rr(rr>0 & rr<1)];
        end
        if axialStep~=0
            r=-axial/axialStep; if r>0 && r<1, breaks=[breaks;r]; end
        end
        breaks=unique(sort(breaks));
        for k=1:numel(breaks)-1
            lo=breaks(k); hi=breaks(k+1); mid=aa+((lo+hi)/2)*dd;
            if b.'*mid>=c*norm(mid)-1e-12*max(1,norm(mid))
                t=min(hi,max(lo,-aa.'*dd/(dd.'*dd)));
                relative=aa+t*dd; distance=norm(relative);
                if distance>0
                    hit=accept(hit,distance,o+relative,n,false,false,NaN,NaN);
                end
            end
        end
        % A cone can touch an edge at a single isolated feasible point.
        for t=breaks.'
            relative=aa+t*dd; distance=norm(relative);
            if distance>0 && b.'*relative>=c*distance-1e-12*max(1,distance)
                hit=accept(hit,distance,o+relative,n,false,false,NaN,NaN);
            end
        end
    end
end
assert(isfinite(hit.range),'libekf:TerrainCoverage','Footprint has no surface within supplied terrain.');
end

function h=accept(h,d,p,n,smooth,active,cosine,planeDistance)
if ~isfinite(h.range) || d<h.range-1e-10*max(1,h.range)
    h=struct('range',d,'point',p,'normal',n,'smooth',smooth,'cone_active',active, ...
        'direction_cosine',cosine,'plane_distance',planeDistance);
elseif abs(d-h.range)<=1e-10*max(1,h.range)
    % Repeated coplanar mesh representation is harmless. Distinct closest
    % points/surface normals do not support a single differentiable return.
    if norm(p-h.point)>1e-8 || abs(abs(n.'*h.normal)-1)>1e-10
        h.smooth=false;
    end
end
end
