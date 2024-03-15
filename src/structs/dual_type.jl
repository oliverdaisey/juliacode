"""
Model for the dual pictures that have been implemented so far.
"""
abstract type DualType end

struct Hypersurface <: DualType end
struct Linear <: DualType end
struct InvertedLinear <: DualType end