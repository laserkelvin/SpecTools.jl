

using Peaks, DataFrames

"""
Function to grab transitions that fit into the range of observations, and then grab
chunks of spectra that overlap with a transition, given a window size.

This function takes a `Spectrum` type, nominally an `Experiment` object, and works
with a `transitions` DataFrame that contains ν and intensity keys, corresponding to
the PGopher output.

"""
function match_spectrum_transitions(
        obs::Spectrum, transitions::DataFrame; velocity::AbstractFloat=0., window_size::AbstractFloat=10., threshold::AbstractFloat=1e-3, interloper_threshold::AbstractFloat=0.1
    )
    rest_ν = transitions.ν
    @show min_obs, max_obs = minimum(obs.ν), maximum(obs.ν)
    @show max_int = maximum(transitions.intensity)
    # isolate transitions that are actually inband and match our intensity cutoff
    trans_mask = (min_obs .< transitions.ν .< max_obs) .& (transitions.intensity .>= max_int * threshold)
    @show sum(trans_mask)
    # calculate the transition frequencies shifted to the nominal source velocity
    offset_ν = rest_ν[trans_mask] .+ doppler_shift_ν.(rest_ν[trans_mask], velocity)
    obs_mask = BitArray(zero.(obs.ν))
    noise = zero.(obs.intensity)
    bar = Progress(length(offset_ν), 1, "Running through catalog...", 50)
    actual_trans = []
    for (index, ν) in enumerate(offset_ν)
        next!(bar)
        velocity_window = doppler_shift_ν(ν, window_size)
        # find regions of spectra that fall into our velocity windows
        lower_ν, upper_ν = ν - velocity_window, ν + velocity_window
        lower_idx, upper_idx = searchsortednearest(obs.ν, lower_ν, 1.), searchsortednearest(obs.ν, upper_ν, 1.)
        if ~isnothing(lower_idx) && ~isnothing(upper_idx)
            # look for peaks; if any are found then we immediately ignore
            peak_idx, peak_vals = peakprom(Maxima(), obs.intensity[lower_idx:upper_idx], minprom=interloper_threshold)
            if length(peak_idx) == 0
                obs_mask[lower_idx:upper_idx] .= 1
                chunk_int = obs.intensity[lower_idx:upper_idx]
                # update the noise array with the estimated values
                update_noise!(noise, chunk_int, lower_idx:upper_idx)
                # add this to the list of transitions we actually used
                push!(actual_trans, index)
            else
                @info "Ignoring transition at $ν due to interloper."
            end
        end
    end
    @show sum(obs_mask)
    new_spectrum = Experiment(
        obs.ν[obs_mask],
        obs.intensity[obs_mask],
        noise[obs_mask],
        ""
    )
    return new_spectrum, transitions[trans_mask,:][actual_trans,:]
end


function update_noise!(noise, chunk_int, mask)
    noise_std = std(chunk_int)
    noise[mask] .= sqrt.(noise_std^2 .+ (chunk_int .* 0.1).^2)
    nothing
end