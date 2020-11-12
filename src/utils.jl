
using YAML: load_file
using Logging

export
    searchsortednearest, setup_logging, get_settings

function searchsortednearest(a, x, threshold::AbstractFloat=1.)
    idx = searchsortedfirst(a,x)
    # if the returned index is actually nowhere near our target
    # then we just return nothing
    if abs(a[idx] - x) > threshold
        return nothing
    end
    if (idx==1); return idx; end
    if (idx>length(a)); return length(a); end
    if (a[idx]==x); return idx; end
    if (abs(a[idx]-x) < abs(a[idx-1]-x))
      return idx
    else
      return idx-1
    end
end

# Grab a YAML file with settings
get_settings(yaml_file::String="settings.yml") = load_file(yaml_file)

function setup_logging(name::String)
    log_file = open("$name.log", "a+")
    global_logger(SimpleLogger(log_file))
    nothing
end

setup_logging(yaml_dict::Dict) = setup_logging(yaml_dict["name"])