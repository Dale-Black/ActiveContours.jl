"""
    region_seg(img, init_mask, max_iterations, alpha=0.1, display)

This is adapted from https://www.mathworks.com/matlabcentral/fileexchange/19567-active-contour-segmentation.
The code implements the paper "Active Contours Without Edges" By Chan Vese. 
"""
function region_seg(img, init_mask, max_iterations, α=0.1, display)
    # Compute the signed distance map from the init_mask
    ϕ = compute_dtm(init_mask)

    ϵ = 1e-3
    for its=1:max_iterations
        idx = findall(ϕ <= 1.2 && ϕ >= -1.2);  #get the curve's narrow band
        
        # find interior and exterior mean
        upts = findall(ϕ <= 0);                     # interior points
        vpts = findall(ϕ > 0);                      # exterior points
        u = sum(img[upts]) / (length(upts) + ϵ);    # interior mean
        v = sum(img[vpts]) / (length(vpts) + ϵ);    # exterior mean

        F = (img[idx]-u).^2 - (img[idx]-v).^2;      # force from image information
        curvature = get_curvature(ϕ, idx);          # force from curvature penalty
    
        dϕ = F ./ max(abs(F)) + α * curvature;      # gradient descent to minimize energy
    end
end