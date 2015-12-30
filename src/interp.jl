type InterpKDE{K,G}
    kde::K
    grid::G
end

function InterpKDE{BC<:BoundaryCondition, IT<:InterpType}(k::UnivariateKDE, bc::Type{BC}=BCnan,it::Type{IT}=InterpQuadratic)
    g = CoordInterpGrid(k.x, k.density, bc, it)
    InterpKDE(k,g)
end
function InterpKDE{BC<:BoundaryCondition, IT<:InterpType}(k::BivariateKDE, bc::Type{BC}=BCnan,it::Type{IT}=InterpQuadratic)
    g = CoordInterpGrid((k.x,k.y),k.density,bc,it)
    InterpKDE{typeof(k),typeof(g)}(k,g)
end


pdf(ik::InterpKDE,x...) = ik.grid[x...]

pdf(k::UnivariateKDE,x) = pdf(InterpKDE(k),x)
pdf(k::BivariateKDE,x,y) = pdf(InterpKDE(k),x,y)

function rand(k::UnivariateKDE, range, min)
  while true
    x = rand() * range + min
    b = rand()
    if b <= pdf(k, x)
      return x
    end
  end
end

function rand(k::BivariateKDE, xrange, xmin, yrange, ymin)
  while true
    x = rand() * xrange + xmin
    y = rand() * yrange + ymin
    b = rand()
    if b <= pdf(k, x, y)
      return (x, y)
    end
  end
end
