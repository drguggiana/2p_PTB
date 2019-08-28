function [long, lat]=patches_deg(n_stim, size_fov, offset_fov, size_patch_rel)
    % [lo,la]=patches_deg([3,4], [100,100], [-50,-50], 0.5)
    
    %convert to radians
    size_fov = (size_fov ./ 180) .* pi;
    offset_fov = (offset_fov ./ 180) .* pi;
    
    %calculate the center x,y points of the patches
    start_cent = offset_fov + size_fov./n_stim./2;  % center first patch
    inc_cent = size_fov./n_stim;                    % distance patch centers
    end_cent = offset_fov + size_fov - size_fov./n_stim./2;  % center last patch
    cent_lo = start_cent(1):inc_cent(1):end_cent(1);  % centers along x axis
    cent_la = start_cent(2):inc_cent(2):end_cent(2);  % centers along y axis
    
    % vectors with all corner points:
    long = zeros(prod(n_stim),4);
    lat = zeros(prod(n_stim),4);
    i=1;
    for x=1:n_stim(1) 
        for y=1:n_stim(2)
            long(i, [1,4]) = cent_lo(x)-inc_cent(1)*size_patch_rel/2; 
            long(i, [2,3]) = cent_lo(x)+inc_cent(1)*size_patch_rel/2; 
            lat(i, [1,2]) = cent_la(y)-inc_cent(2)*size_patch_rel/2; 
            lat(i, [3,4]) = cent_la(y)+inc_cent(2)*size_patch_rel/2; 
            i=i+1;
        end
    end
        
end