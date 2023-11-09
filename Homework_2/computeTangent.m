function tangent = computeTangent( q )

% Input: q is a vector of size 4*nv - 1
% Output: tangent is a matrix of size (ne, 3) where ne = nv - 1
nv = (length(q) + 1) / 4; % length(q) = 4*nv - 1
ne = nv - 1;
tangent = zeros(ne, 3);

for c=1:ne % Loop over each edge
    xc = q(4*c-3 : 4*c-1);
    xcp1 = q(4*c+1 : 4*c+3);
    dx = xcp1 - xc;
    tangent(c,:) = dx / norm(dx);
end

end