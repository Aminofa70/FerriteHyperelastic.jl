# src/InputStruc.jl
"""
Create a dynamic structure for FEM parameters
"""
mutable struct InputStruct
    parameters::Dict{Symbol, Any}
    function InputStruct()
        new(Dict{Symbol, Any}())
    end
end

function Base.getproperty(obj::InputStruct, name::Symbol)
    if name == :parameters || haskey(obj.parameters, name)
        return name == :parameters ? getfield(obj, :parameters) : obj.parameters[name]
    else
        throw(ArgumentError("Property '$name' does not exist."))
    end
end

function Base.setproperty!(obj::InputStruct, name::Symbol, value)
    if name == :parameters
        throw(ArgumentError("Cannot modify the `parameters` field directly."))
    else
        obj.parameters[name] = value
    end
end

