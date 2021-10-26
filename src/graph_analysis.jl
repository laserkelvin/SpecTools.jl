
"""a is a vector of length M, b is a matrix shape M x N
"""
function encoding_distance(a, b)
    @tullio distances[n] := a[m] - b[m,n]
end


function link_nodes_by_Δ!(sg::SpectroscopicGraph; Δ=1)
    encodings = reduce(hcat, sg.levels.encoding)
    for node ∈ 1:length(sg.levels)
        distances = abs.(encoding_distance(encodings[:,node], encodings))
        for (i, d) in enumerate(distances)
            d != Δ ? continue : MetaGraphs.add_edge!(sg.graph, node, i, Dict(:Δ => d))
        end
    end
end