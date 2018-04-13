export tod2map, baseline2map, applyz, applyz_and_sum, applya, destriped_map, destripe

using LinearMaps

import Healpix
import IterativeSolvers


#from TOD to binned map: it computes the mean value of the samples taken in each pixel
function tod2map(pix_idx, tod, num_of_pixels; unseen = NaN)
    binned_map = zeros(Float64,num_of_pixels)
    hits = zeros(Int, num_of_pixels)
    
    for i in eachindex(pix_idx)
        binned_map[pix_idx[i]] = binned_map[pix_idx[i]] + tod[i]
        hits[pix_idx[i]] = hits[pix_idx[i]] + 1
    end

    for i in eachindex(binned_map)
        if (hits[i] != 0)
            binned_map[i] = binned_map[i] / hits[i]
        else
            binned_map[i] = unseen
        end
    end 
    
    binned_map
end


function baseline2map(pix_idx, baselines, baseline_dim, num_of_pixels; unseen=NaN)
    noise_map = zeros(Float64, num_of_pixels)
    hits = zeros(Int, num_of_pixels)
    startidx = 1
    
    for i in eachindex(baseline_dim)
        endidx = baseline_dim[i] + startidx - 1        
        
        for j in startidx:endidx
        noise_map[pix_idx[j]] = noise_map[pix_idx[j]] +  baselines[i]
        hits[pix_idx[j]] = hits[pix_idx[j]] + 1
        end
        startidx += baseline_dim[i]
    end 
    
    for i in eachindex(noise_map)
        if (hits[i] != 0)
            noise_map[i] = noise_map[i] / hits[i]
        else
            noise_map[i] = unseen
        end
    end 
    
    noise_map
end



function applyz_and_sum(pix_idx, tod, baseline_dim, num_of_pixels; unseen=NaN)
    @assert length(tod) == length(pix_idx)

    binned_map = tod2map(pix_idx, tod, num_of_pixels, unseen=unseen)

    startidx = 1
    baselines_sum = zeros(Float64, length(baseline_dim))

    for i in eachindex(baseline_dim)
        endidx = baseline_dim[i] + startidx - 1

        # The inner for is equivalent to
        #
        #   baselines_sum[i] += sum(tod[startidx:endidx] - binned_map[pix_idx[startidx:endidx]])
        #
        # but roundoff errors are reduced
        for j in startidx:endidx
            baselines_sum[i] += tod[j] - binned_map[pix_idx[j]]
        end

        startidx += baseline_dim[i]
    end
    
    baselines_sum
end



function applya(baselines, pix_idx, tod, baseline_dim, num_of_pixels; unseen=NaN)
    @assert length(tod) == length(pix_idx)

    binned_map = baseline2map(pix_idx, baselines, baseline_dim, num_of_pixels; unseen=NaN)
    startidx = 1
    baselines_sum = zeros(Float64, length(baseline_dim))

    for i in eachindex(baseline_dim)
        endidx = baseline_dim[i] + startidx - 1
        
        for j in startidx:endidx
            baselines_sum[i] += baselines[i] - binned_map[pix_idx[j]]
        end

        startidx += baseline_dim[i]
    end
    
    baselines_sum
end



function destriped_map(baselines, pix_idx, tod, baseline_dim, num_of_pixels; unseen=NaN)
    @assert length(tod) == length(pix_idx)
    tod2map(pix_idx, tod, num_of_pixels) - baseline2map(pix_idx, baselines, baseline_dim, num_of_pixels)
end



function destripe(tod, pix_idx, baseline_dim, num_of_pixels; unseen=NaN)
    @assert sum(baseline_dim) == length(tod)

    A = LinearMap(x -> applya(x, pix_idx, tod, baseline_dim, num_of_pixels),
                  length(baseline_dim), length(baseline_dim),
                  issymmetric = true, ishermitian = true, isposdef = true)

    baselines_sum = applyz_and_sum(pix_idx, tod, baseline_dim, num_of_pixels, unseen=unseen)
    baselines = IterativeSolvers.cg(A, baselines_sum)

    # once we have an estimate of the baselines, we can build the destriped map
    pixels = destriped_map(baselines, pix_idx, tod, baseline_dim, num_of_pixels, unseen=unseen)

    (pixels, baselines)
end