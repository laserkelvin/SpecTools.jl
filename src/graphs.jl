
"""Provides an abstract type for spectroscopic graphs,
allowing for other representations of spectroscopic graphs
to be developed.
"""
abstract type SpectroscopicGraph end

"""A projected spectroscopic monograph, where nodes represent
energy levels, and edges between nodes represent transitions
that connect them.
"""
struct ProjectedSG <: SpectroscopicGraph
  levels
  transitions
  graph
end

"""A bipartite spectroscopic graph, with two classes of nodes:
energy levels and transitions. Transition nodes are connected to
two energy levels nodes.
"""
struct BipartiteSG <: SpectroscopicGraph
  levels
  transitions
  graph
end

Base.length(sg::SpectroscopicGraph) = length(sg.levels) + length(sg.transitions)
offset(sg::SpectroscopicGraph) = length(sg.levels)
indices(sg::SpectroscopicGraph) = 1:length(sg)

function ProjectedSG(ls, ts)
	sg = MetaGraph(length(ls))
	for t ∈ ts
        properties = Dict(:I => t.I, :ν => t.ν)
		lq, uq = lower_state_qnos(t), upper_state_qnos(t)
		l_index, u_index = match_level(ls, lq), match_level(ls, uq)
        MetaGraphs.add_edge!(sg, l_index, u_index, properties)
	end
	return ProjectedSG(ls, ts, sg)	
end

function BipartiteSG(ls, ts)
	sg = MetaGraph(length(ls) + length(ts))
	for (t_index, t) ∈ enumerate(ts)
		t_offset = length(ls) + t_index + 1
        # get the transition strength and frequency, potentially for plotting
        properties = Dict(:I => t.I, :ν => t.ν)
		lq, uq = lower_state_qnos(t), upper_state_qnos(t)
		l_index, u_index = match_level(ls, lq), match_level(ls, uq)
		for index in [l_index, u_index]
          MetaGraphs.add_edge!(sg, t_offset, index, properties)
		end
	end
	return BipartiteSG(ls, ts, sg)
end