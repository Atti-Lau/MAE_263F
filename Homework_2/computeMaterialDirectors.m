function [m1, m2] = computeMaterialDirectors(a1, a2, theta)

ne = length(theta);
m1 = zeros(ne, 3); % Material frame director 1 for all the edges
m2 = zeros(ne, 3); % Material frame director 2 for all the edges

for c=1:ne % loop over the edges
    cs = cos(theta(c));
    ss = sin(theta(c));
    m1(c,:) = cs * a1(c,:) + ss * a2(c,:);
    m2(c,:) = - ss * a1(c,:) + cs * a2(c,:);
end

end