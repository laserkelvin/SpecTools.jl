

function search_frequency(rest_ν, ν)
  min_i = 1
  min_Δ = Inf
  for (i, value) in enumerate(ν)
    current_Δ = abs(value - rest_ν)
    (min_Δ, min_i) = current_Δ < min_Δ ? (current_Δ, i) : (min_Δ, min_i)
  end
  return min_Δ, min_i
end


"""Loops through two vectors A and B, searching for the best matches
for between them.

Matches are only considered if the difference of the values meet
a preliminary threshold, as a fraction of the value of A.

This is the preferred method of matching the catalog and obs. frequencies,
as it means you only run the loop once.

Returns a list of index pairs of A and B.
"""
function intersect_frequencies(A::Vector{T}, B::Vector{S}; thres=1f-4) where {T<:Real,S<:Real}
  matches = []
  # preallocate up to the length of the rest frequencies
  sizehint!(matches, length(A))
  for (i, rest_value) in enumerate(A)
    to_match = rest_value * thres
    best_diff = Inf
    best_match = nothing
    for (j, obs_value) in enumerate(B)
      diff = abs(rest_value - obs_value)
      if diff <= to_match && diff < best_diff
        best_diff = diff
        best_match = [i, j]
      end
      # if the diff increases from the best diff, we've passed it
      if diff > best_diff
        push!(matches, best_match)
        break
      end
    end
  end
  return matches
end

#function extract_windows(rest_ν, obs::Observation, vlsr; thres=1f-4,)
#  matches = intersect_frequencies(rest_ν, obs.frequency; thres=thres)
#end


