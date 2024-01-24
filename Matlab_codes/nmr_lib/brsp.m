function [spec] = brsp(struc)
%UNTITLED2 Summary of this function goes here
%   criar spec bruker

for k = 1:1:length(struc.sim.fid{1})
spec.fid(k) = struc.sim.fid{1}(k);
end

for k = 1:1:length(struc.sim.esp{1})
spec.esp(k) = struc.sim.esp{1}(k);
end

spec.sw = struc.sw;

end

