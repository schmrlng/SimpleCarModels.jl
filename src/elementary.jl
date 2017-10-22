export elementary, bi_elementary
export elementary_waypoints, bi_elementary_waypoints

function elementary(q1::SE2State{T}, q2::SE2State{T}, v::T = T(1), λ::T = T(.5),
                    ψ::T = atan2(q2.y - q1.y, q2.x - q1.x)) where {T<:AbstractFloat}
    adiff(ψ, q1.θ) ≈ adiff(q2.θ, ψ) || error("$q1 and $q2 are not symmetric for elementary path")
    abs(adiff(ψ, q1.θ)) < T(π/2) || error("$q2 is not in front of $q1 for elementary path")
    α = adiff(q2.θ, q1.θ)
    d = sqrt((q2.x - q1.x)^2 + (q2.y - q1.y)^2)
    α ≈ 0 && return SVector(VelocityCurvRateStep(d/v, VelocityCurvRateControl(v, T(0))))
    L = d / D(α, λ)
    σ = 4*α / (L*L*(1 - λ*λ))
    ctrl = SVector(VelocityCurvRateStep(L*(1-λ)/2, VelocityCurvRateControl(T(1), σ)),
                   VelocityCurvRateStep(L*λ, VelocityCurvRateControl(T(1), T(0))),
                   VelocityCurvRateStep(L*(1-λ)/2, VelocityCurvRateControl(T(1), -σ)))
    scalespeed.(ctrl, v)
end

function bi_elementary(q1::SE2State{T}, q2::SE2State{T}, v::T = T(1), γ::T = T(.5), λ::T = T(.5)) where {T<:AbstractFloat}
    adiff(q1.θ, q2.θ) ≈ 0 || error("angles of $q1 and $q2 must be equal for bi-elementary path")
    q1 ≈ q2 && return SVector(zero(VelocityCurvRateStep{T}))
    ψ = atan2(q2.y - q1.y, q2.x - q1.x)
    θi = mod2piF(2*ψ - q1.θ)
    qi = SE2State(q1.x + γ*(q2.x - q1.x), q1.y + γ*(q2.y - q1.y), θi)
    [elementary(q1, qi, v, λ, ψ); elementary(qi, q2, v, λ, ψ)]
end

function elementary_waypoints(q1::SE2State{T}, q2::SE2State{T}, dt_or_N::Union{T,Int},
                              v::T = T(1), λ::T = T(.5)) where {T}
    ctrl = elementary(q1, q2, v, λ)
    waypoints(SimpleCarDynamics{0,1,T}(), SE2κState(q1), ctrl, dt_or_N)
end

function bi_elementary_waypoints(q1::SE2State{T}, q2::SE2State{T}, dt_or_N::Union{T,Int},
                                 v::T = T(1), γ::T = T(.5), λ::T = T(.5)) where {T}
    ctrl = bi_elementary(q1, q2, v, γ, λ)
    waypoints(SimpleCarDynamics{0,1,T}(), SE2κState(q1), ctrl, dt_or_N)
end

@inline function D(α::T, λ::T) where {T}
    α = abs(α)
    c = α/(1 + λ)
    s = sqrt(α/(T(π)*(1 - λ*λ)))
    sin(c*λ)/c + (cos(α/2)*fresnelC(s*(1-λ)) + sin(α/2)*fresnelS(s*(1-λ)))/s
end