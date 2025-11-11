
mutable struct KernelEstimate{VF<:VariateForm, VS<:ValueSupport, KernelType<:Distribution, PriorDist<:UnivariateDistribution{Discrete}, PDF<:DiscretisedPDF} <: AbstractMixtureModel{VF, VS, KernelType}
    const data::Matrix{Float64}
    const prior::PriorDist
    const kernel::KernelType
    const bandwidth::Float64
    precomputedPDF::Union{Nothing, PDF}

    # these constructors guarantee type agreement
    function KernelEstimate(data::Matrix{Float64}, prior::PriorDist, kernel::KernelType, bandwidth::Float64, precomputedPDF::PDF) where {
            PriorDist<:UnivariateDistribution{Discrete},
            KernelType <: Distribution,
            PDF <: DiscretisedPDF}
        VF, VS = supertype(KernelType).parameters
        new{VF,VS,KernelType,PriorDist,PDF}(data, prior, kernel, bandwidth, precomputedPDF)
    end
    function KernelEstimate(data::Matrix{Float64}, prior::PriorDist, kernel::KernelType, bandwidth::Float64, ::Nothing) where {
            PriorDist<:UnivariateDistribution{Discrete},
            KernelType <: Distribution}
        VF, VS = supertype(KernelType).parameters

        # the default PDF type is based on Float64 and Int
        R = Base.return_types(range,(Float64,Float64,Int))[1] 
        PDF = DiscretisedPDF{size(data)[1],R,Float64}
        new{VF,VS,KernelType,PriorDist,PDF}(data, prior, kernel, bandwidth, nothing)
    end

end