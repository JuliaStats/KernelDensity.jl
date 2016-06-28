type InterpKDE{K, G}
    kde::K
    grid::G
end

function InterpKDE{BC<:BoundaryCondition, IT<:InterpType}(k::UnivariateKDE, bc::Type{BC}=BCnan, it::Type{IT}=InterpQuadratic)
    g = CoordInterpGrid(k.x, k.density, bc, it)
    InterpKDE(k, g)
end
function InterpKDE{BC<:BoundaryCondition, IT<:InterpType}(k::BivariateKDE, bc::Type{BC}=BCnan, it::Type{IT}=InterpQuadratic)
    g = CoordInterpGrid((k.x, k.y), k.density, bc, it)
    InterpKDE{typeof(k),typeof(g)}(k, g)
end

# Support BCfill
function InterpKDE{IT<:InterpType}(k::UnivariateKDE, bc::Number, it::Type{IT}=InterpQuadratic)
    g = CoordInterpGrid(k.x, k.density, bc, it)
    InterpKDE(k, g)
end
function InterpKDE{IT<:InterpType}(k::BivariateKDE, bc::Number, it::Type{IT}=InterpQuadratic)
    g = CoordInterpGrid((k.x, k.y), k.density, bc, it)
    InterpKDE{typeof(k),typeof(g)}(k, g)
end


pdf(ik::InterpKDE, x...) = ik.grid[x...]

pdf(k::UnivariateKDE, x) = pdf(InterpKDE(k), x)
pdf(k::BivariateKDE, x, y) = pdf(InterpKDE(k), x, y)
