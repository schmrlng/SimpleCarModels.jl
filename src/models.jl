import DifferentialDynamicsModels: SteeringBVP, propagate, state_dim, control_dim

export SimpleCarDynamics
export SE2State, SE2vState, SE2κState, SE2vκState
export VelocityCurvatureControl, AccelerationCurvatureControl, VelocityCurvRateControl, AccelerationCurvRateControl
export VelocityCurvatureStep, VelocityCurvRateStep

### Simple Car Dynamics (with integrators in speed v and curvature κ)
struct SimpleCarDynamics{Dv,Dκ} <: DifferentialDynamics end

state_dim(::SimpleCarDynamics{Dv,Dκ}) where {Dv,Dκ} = Dv + Dκ + 3
control_dim(::SimpleCarDynamics{Dv,Dκ}) where {Dv,Dκ} = 2

struct SE2State{T<:AbstractFloat} <: FieldVector{3,T} x::T; y::T; θ::T end
struct SE2vState{T<:AbstractFloat} <: FieldVector{4,T} x::T; y::T; θ::T; v::T end
struct SE2κState{T<:AbstractFloat} <: FieldVector{4,T} x::T; y::T; θ::T; κ::T end
struct SE2vκState{T<:AbstractFloat} <: FieldVector{5,T} x::T; y::T; θ::T; v::T; κ::T end

struct VelocityCurvatureControl{T<:AbstractFloat} <: FieldVector{2,T}; v::T; κ::T end
struct AccelerationCurvatureControl{T<:AbstractFloat} <: FieldVector{2,T}; a::T; κ::T end
struct VelocityCurvRateControl{T<:AbstractFloat} <: FieldVector{2,T}; v::T; σ::T end
struct AccelerationCurvRateControl{T<:AbstractFloat} <: FieldVector{2,T}; a::T; σ::T end

for V in (SE2State, SE2vState, SE2κState, SE2vκState,    # required for StaticArrays to keep type through +,*,... operations
          VelocityCurvatureControl, AccelerationCurvatureControl, VelocityCurvRateControl, AccelerationCurvRateControl)
    @eval StaticArrays.similar_type(::Type{A}, ::Type{T}, ::Size{$(size(V))}) where {A <: $V, T} = $V{T}
end
const VelocityCurvatureStep{T} = StepControl{2,T,VelocityCurvatureControl{T}}    # alias used for dubins, reedsshepp
const VelocityCurvRateStep{T} = StepControl{2,T,VelocityCurvRateControl{T}}      # alias used for dubinsCC, elementary

### Dv = 0, Dκ = 0 (Dubins Car, Reeds-Shepp Car)
function propagate(f::SimpleCarDynamics{0,0}, q::StaticVector{3,T}, c::StepControl{2,T}) where {T}
    t, v, κ = c.t, c.u[1], c.u[2]
    x, y, θ = q[1], q[2], q[3]
    s, c = sincos(θ)
    if abs(κ) > sqrt(eps(T))
        st, ct = sincos(θ + v*κ*t)
        typeof(q)(x + (st - s)/κ,
                  y + (c - ct)/κ,
                  mod2piF(θ + v*κ*t))
    else
        typeof(q)(x + v*c*t - v*v*s*κ*t*t/2,
                  y + v*s*t + v*v*c*κ*t*t/2,
                  mod2piF(θ + v*κ*t))
    end
end

### Dv = 1, Dκ = 0 (Dubins Car, Reeds-Shepp Car with Acceleration)
function propagate(f::SimpleCarDynamics{1,0}, q::StaticVector{4,T}, c::StepControl{2,T}) where {T}
    t, a, κ = c.t, c.u[1], c.u[2]
    x, y, θ, v = q[1], q[2], q[3], q[4]
    s = v*t + a*t*t/2
    xyθ = propagate(SimpleCarDynamics{0,0}(), SE2State(x, y, θ), StepControl(s, VelocityCurvatureControl(T(1), κ)))
    typeof(q)(xyθ[1], xyθ[2], xyθ[3], v + a*t)
end

### Dv = 0, Dκ = 1 (Continuous Curvature Car)
function propagate(f::SimpleCarDynamics{0,1}, q::StaticVector{4,T}, c::StepControl{2,T}) where {T}
    t, v, σ = c.t, c.u[1], c.u[2]
    x, y, θ, κ = q[1], q[2], q[3], q[4]
    if abs(σ) > sqrt(eps(T))
        spi = sqrt(T(pi))
        θ_ = θ - v*κ*κ/σ/2
        s, c = sincos(θ_)
        sv_σ = sqrt(abs(v/σ))
        PK = sv_σ*κ/spi
        PT = sv_σ*σ*t/spi
        FCK = flipsign(flipsign(fresnelC(PK), σ), v)
        FCT = flipsign(flipsign(fresnelC(PK + PT), σ), v)
        FSK = fresnelS(PK)
        FST = fresnelS(PK + PT)
        typeof(q)(x + spi*sv_σ*(c*(FCT-FCK) + s*(FSK-FST)),
                  y + spi*sv_σ*(c*(FST-FSK) + s*(FCT-FCK)),
                  mod2piF(θ + v*κ*t + v*σ*t*t/2),
                  κ + σ*t)
    elseif abs(κ) > sqrt(eps(T))
        s, c = sincos(θ)
        st, ct = sincos(θ + v*κ*t)
        typeof(q)(x + (st - s)/κ,    # perhaps include higher order σ terms
                  y + (c - ct)/κ,
                  mod2piF(θ + v*κ*t + v*σ*t*t/2),
                  κ + σ*t)
    else
        s, c = sincos(θ)
        typeof(q)(x + v*c*t,    # perhaps include higher order κ, σ terms
                  y + v*s*t,
                  mod2piF(θ + v*κ*t + v*σ*t*t/2),
                  κ + σ*t)
    end
end

### General Simple Car Dynamics (ẋ = f(x, u))

function (::SimpleCarDynamics{0,0})(q::StaticVector{3}, u::StaticVector{2})
    v, κ = u[1], u[2]
    x, y, θ = q[1], q[2], q[3]
    SVector(v*cos(θ), v*sin(θ), v*κ)
end

function (::SimpleCarDynamics{1,0})(q::StaticVector{4}, u::StaticVector{2})
    a, κ = u[1], u[2]
    x, y, θ, v = q[1], q[2], q[3], q[4]
    SVector(v*cos(θ), v*sin(θ), v*κ, a)
end

function (::SimpleCarDynamics{0,1})(q::StaticVector{4}, u::StaticVector{2})
    v, σ = u[1], u[2]
    x, y, θ, κ = q[1], q[2], q[3], q[4]
    SVector(v*cos(θ), v*sin(θ), v*κ, σ)
end

function (::SimpleCarDynamics{1,1})(q::StaticVector{5}, u::StaticVector{2})
    a, σ = u[1], u[2]
    x, y, θ, v, κ = q[1], q[2], q[3], q[4], q[5]
    SVector(v*cos(θ), v*sin(θ), v*κ, a, σ)
end

@generated function (::SimpleCarDynamics{Dv,Dκ})(q::StaticVector{N}, u::StaticVector{2}, x::Int) where {Dv,Dκ,N}
    @assert Dv + Dκ + 3 == N
    θ = :(q[3])
    v = :(q[4])
    κ = :(q[$(4+Dv)])
    qdot = [:($v*cos($θ)); :($v*sin($θ)); :($v*$κ);
            [:(q[$(4+i)]) for i in 1:Dv-1]; :(u[1]);
            [:(q[$(4+Dv+i)]) for i in 1:Dκ-1]; :(u[2])]
    return quote
        SVector{N}(tuple($(qdot...)))
    end
end

# ### Steering
# abstract type SimpleCarDirectionality end
# struct ForwardsOnly <: SimpleCarDirectionality end
# struct Bidirectional <: SimpleCarDirectionality end
# struct TurningRadius{D<:SimpleCarDirectionality,T<:AbstractFloat} <: SteeringParams
#     r::T
# end
# struct MaxCurvature{D<:SimpleCarDirectionality,T<:AbstractFloat} <: SteeringParams
#     κ_max::T
# end
# TurningRadius(p::MaxCurvature{D,T}) where {D,T} = TurningRadius{D,T}(1/p.κ_max)
# MaxCurvature(p::TurningRadius{D,T}) where {D,T} = MaxCurvature{D,T}(1/p.κ_max)
#
# struct MaxCurvatureAndCurvatureRate{D<:SimpleCarDirectionality,T<:AbstractFloat} <: SteeringParams
#     κ_max::T
#     σ_max::T
# end
#
# function SteeringBVP(f::SimpleCarDynamics{0,0,T}, j::Time,
#                      r::T, ::Type{D} = ForwardsOnly) where {T,D<:SimpleCarDirectionality}
#     SteeringBVP(f, j, TurningRadius{D,T}(r))
# end
# function SteeringBVP(f::SimpleCarDynamics{0,1,T}, j::Time,
#                      κ_max::T, σ_max::T, ::Type{D} = ForwardsOnly) where {T,D<:SimpleCarDirectionality}
#     SteeringBVP(f, j, TurningRadius{D,T}(r))
# end
