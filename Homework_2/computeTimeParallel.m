function [a1, a2] = computeTimeParallel(a1_old, q0, q)

ne = (length(q) + 1) / 4 - 1; % Number of edges
tangent0 = computeTangent( q0 ); % Old tangents : ne x 3
tangent = computeTangent( q ); % New tangents : ne x 3
a1 = zeros(ne, 3); % First reference frame director
a2 = zeros(ne, 3); % Second reference frame director

for c=1:ne % Loop over edges
    t0 = tangent0(c,:); % Old tangent on c-th edge
    t = tangent(c,:); % New tangent on c-th edge
    a1_local = parallel_transport( a1_old(c,:), t0, t );
    
    % Just to be careful
    a1_local = a1_local - dot(a1_local, t) * t; % Enforcing a1 is perp. to t
    a1_local = a1_local / norm(a1_local); % Enforcing unit vector
    % Store
    a1(c,:) = a1_local;
    a2(c,:) = cross(t, a1_local);
end

end