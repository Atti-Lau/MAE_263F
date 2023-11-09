function [Fs, Js] = getFs(q)

global EA refLen

nv = (length(q)+1) / 4;
ne = nv - 1;
Fs = zeros(size(q));
Js = zeros(length(q), length(q));

for c = 1:ne
    node0 = [q(4*c-3), q(4*c-2), q(4*c-1)];
    node1 = [q(4*c+1), q(4*c+2), q(4*c+3)];
    
    [dF, dJ] = ...
        gradEs_hessEs(node0, node1, ...
        refLen(c), EA);
    
    ind = [4*c-3, 4*c-2, 4*c-1, 4*c+1, 4*c+2, 4*c+3]; % Size 6
    
    Fs(ind) = Fs(ind) - dF;
    Js(ind, ind) = Js(ind, ind) - dJ;
end

end