export reedsshepp, reedsshepp_length, reedsshepp_waypoints

reedsshepp_length(q1::SE2State{T}, q2::SE2State{T}, r::T = T(1)) where {T} = reedsshepp(q1, q2, r)[1]
function reedsshepp_waypoints(q1::SE2State{T}, q2::SE2State{T}, dt_or_N::Union{T,Int}, r::T = T(1), v::T = T(1)) where {T}
    ctrl = reedsshepp(q1, q2, r, v)[2]
    waypoints(SimpleCarDynamics{0,0,T}(), q1, ctrl, dt_or_N)
end

const POST, POST_T, POST_R, POST_B, POST_R_T, POST_B_T, POST_B_R, POST_B_R_T = 0, 1, 2, 3, 4, 5, 6, 7
@inline timeflip(c::VelocityCurvatureStep) = StepControl(c.t, VelocityCurvatureControl(-c.u.v, c.u.κ))
@inline reflect(c::VelocityCurvatureStep) = StepControl(c.t, VelocityCurvatureControl(c.u.v, -c.u.κ))
@inline reverse(ctrl::SVector{5}) = ctrl[SVector(5,4,3,2,1)]

function reedsshepp(q1::SE2State{T},
                    q2::SE2State{T},
                    r::T = T(1),
                    v::T = T(1)) where {T<:AbstractFloat}
    dx = (q2.x - q1.x) / r
    dy = (q2.y - q1.y) / r
    st, ct = sincos(q1.θ)
    target = SE2State(dx*ct + dy*st, -dx*st + dy*ct, mod2piF(q2.θ - q1.θ))

    tTarget = timeflip(target)
    rTarget = reflect(target)
    trTarget = reflect(tTarget)
    bTarget = backwards(target)
    btTarget = timeflip(bTarget)
    brTarget = reflect(bTarget)
    btrTarget = reflect(btTarget)

    c, ctrl, post = T(Inf), zeros(SVector{5,VelocityCurvatureStep{T}}), POST
    # (8.1) C S C
    b, c, ctrl = LpSpLp(target, c, ctrl); b && (post = POST)
    b, c, ctrl = LpSpLp(tTarget, c, ctrl); b && (post = POST_T)
    b, c, ctrl = LpSpLp(rTarget, c, ctrl); b && (post = POST_R)
    b, c, ctrl = LpSpLp(trTarget, c, ctrl); b && (post = POST_R_T)

    # (8.2) C S C
    b, c, ctrl = LpSpRp(target, c, ctrl); b && (post = POST)
    b, c, ctrl = LpSpRp(tTarget, c, ctrl); b && (post = POST_T)
    b, c, ctrl = LpSpRp(rTarget, c, ctrl); b && (post = POST_R)
    b, c, ctrl = LpSpRp(trTarget, c, ctrl); b && (post = POST_R_T)

    # (8.3) C|C|C
    b, c, ctrl = LpRmLp(target, c, ctrl); b && (post = POST)
    # b, c, ctrl = LpRmLp(tTarget, c, ctrl); b && (post = POST_T) # (redundant)
    b, c, ctrl = LpRmLp(rTarget, c, ctrl); b && (post = POST_R)
    # b, c, ctrl = LpRmLp(trTarget, c, ctrl); b && (post = POST_R_T) # (redundant)

    # (8.4) C|C C
    b, c, ctrl = LpRmLm(target, c, ctrl); b && (post = POST)
    b, c, ctrl = LpRmLm(tTarget, c, ctrl); b && (post = POST_T)
    b, c, ctrl = LpRmLm(rTarget, c, ctrl); b && (post = POST_R)
    b, c, ctrl = LpRmLm(trTarget, c, ctrl); b && (post = POST_R_T)
    b, c, ctrl = LpRmLm(bTarget, c, ctrl); b && (post = POST_B)
    b, c, ctrl = LpRmLm(btTarget, c, ctrl); b && (post = POST_B_T)
    b, c, ctrl = LpRmLm(brTarget, c, ctrl); b && (post = POST_B_R)
    b, c, ctrl = LpRmLm(btrTarget, c, ctrl); b && (post = POST_B_R_T)

    # (8.7) C Cu|Cu C
    b, c, ctrl = LpRpuLmuRm(target, c, ctrl); b && (post = POST)
    b, c, ctrl = LpRpuLmuRm(tTarget, c, ctrl); b && (post = POST_T)
    b, c, ctrl = LpRpuLmuRm(rTarget, c, ctrl); b && (post = POST_R)
    b, c, ctrl = LpRpuLmuRm(trTarget, c, ctrl); b && (post = POST_R_T)

    # (8.8) C|Cu Cu|C
    b, c, ctrl = LpRmuLmuRp(target, c, ctrl); b && (post = POST)
    b, c, ctrl = LpRmuLmuRp(tTarget, c, ctrl); b && (post = POST_T)
    b, c, ctrl = LpRmuLmuRp(rTarget, c, ctrl); b && (post = POST_R)
    b, c, ctrl = LpRmuLmuRp(trTarget, c, ctrl); b && (post = POST_R_T)

    # (8.9)
    b, c, ctrl = LpRmSmLm(target, c, ctrl); b && (post = POST)
    b, c, ctrl = LpRmSmLm(tTarget, c, ctrl); b && (post = POST_T)
    b, c, ctrl = LpRmSmLm(rTarget, c, ctrl); b && (post = POST_R)
    b, c, ctrl = LpRmSmLm(trTarget, c, ctrl); b && (post = POST_R_T)
    b, c, ctrl = LpRmSmLm(bTarget, c, ctrl); b && (post = POST_B)
    b, c, ctrl = LpRmSmLm(btTarget, c, ctrl); b && (post = POST_B_T)
    b, c, ctrl = LpRmSmLm(brTarget, c, ctrl); b && (post = POST_B_R)
    b, c, ctrl = LpRmSmLm(btrTarget, c, ctrl); b && (post = POST_B_R_T)

    # (8.10)
    b, c, ctrl = LpRmSmRm(target, c, ctrl); b && (post = POST)
    b, c, ctrl = LpRmSmRm(tTarget, c, ctrl); b && (post = POST_T)
    b, c, ctrl = LpRmSmRm(rTarget, c, ctrl); b && (post = POST_R)
    b, c, ctrl = LpRmSmRm(trTarget, c, ctrl); b && (post = POST_R_T)
    b, c, ctrl = LpRmSmRm(bTarget, c, ctrl); b && (post = POST_B)
    b, c, ctrl = LpRmSmRm(btTarget, c, ctrl); b && (post = POST_B_T)
    b, c, ctrl = LpRmSmRm(brTarget, c, ctrl); b && (post = POST_B_R)
    b, c, ctrl = LpRmSmRm(btrTarget, c, ctrl); b && (post = POST_B_R_T)

    # (8.11) C|Cpi/2 S Cpi/2|C
    b, c, ctrl = LpRmSmLmRp(target, c, ctrl); b && (post = POST)
    b, c, ctrl = LpRmSmLmRp(tTarget, c, ctrl); b && (post = POST_T)
    b, c, ctrl = LpRmSmLmRp(rTarget, c, ctrl); b && (post = POST_R)
    b, c, ctrl = LpRmSmLmRp(trTarget, c, ctrl); b && (post = POST_R_T)

    ctrl = scalespeed.(scaleradius.(ctrl, r), v)
    if post == POST_T
        ctrl = timeflip.(ctrl)
    elseif post == POST_R
        ctrl = reflect.(ctrl)
    elseif post == POST_B
        ctrl = reverse(ctrl)
    elseif post == POST_R_T
        ctrl = reflect.(timeflip.(ctrl))
    elseif post == POST_B_T
        ctrl = timeflip.(ctrl)
        ctrl = reverse(ctrl)
    elseif post == POST_B_R
        ctrl = reflect.(ctrl)
        ctrl = reverse(ctrl)
    elseif post == POST_B_R_T
        ctrl = reflect.(timeflip.(ctrl))
        ctrl = reverse(ctrl)
    end
    c*r, ctrl
end

# Utilities (pedantic about typing to guard against problems)
@inline R(x::T, y::T) where {T} = sqrt(x*x + y*y), atan2(y, x)
@inline function M(t::T) where {T}
    m = mod2piF(t)
    m > pi ? m - 2*T(pi) : m
end
@inline function Tau(u::T, v::T, E::T, N::T) where {T}
    delta = u - v
    A = sin(u) - sin(delta)
    B = cos(u) - cos(delta) - 1
    r, θ = R(E*A + N*B, N*A - E*B)
    t = 2*cos(delta) - 2*cos(v) - 2*cos(u) + 3
    t < 0 ? M(θ + pi) : M(θ)
end
@inline Omega(u::T, v::T, E::T, N::T, t::T) where {T} = M(Tau(u, v, E, N) - u + v - t)
@inline timeflip(q::SE2State{T}) where {T} = SE2State(-q.x, q.y, -q.θ)
@inline reflect(q::SE2State{T}) where {T} = SE2State(q.x, -q.y, -q.θ)
@inline backwards(q::SE2State{T}) where {T} = SE2State(q.x*cos(q.θ) + q.y*sin(q.θ), q.x*sin(q.θ) - q.y*cos(q.θ), q.θ)

@inline function LpSpLp(tgt::SE2State{T}, c::T, ctrl::SVector{5,VelocityCurvatureStep{T}}) where {T}
    tx, ty, tt = tgt.x, tgt.y, tgt.θ
    r, θ = R(tx - sin(tt), ty - 1 + cos(tt))
    u = r
    t = mod2piF(θ)
    v = mod2piF(tt - t)
    cnew = t + u + v
    c <= cnew && return false, c, ctrl
    ctrl = SVector(
        carsegment2stepcontrol(1, t),
        carsegment2stepcontrol(0, u),
        carsegment2stepcontrol(1, v),
        zero(VelocityCurvatureStep{T}),
        zero(VelocityCurvatureStep{T})
    )
    true, cnew, ctrl
end

@inline function LpSpRp(tgt::SE2State{T}, c::T, ctrl::SVector{5,VelocityCurvatureStep{T}}) where {T}
    tx, ty, tt = tgt.x, tgt.y, tgt.θ
    r, θ = R(tx + sin(tt), ty - 1 - cos(tt))
    r*r < 4 && return false, c, ctrl
    u = sqrt(r*r - 4)
    r1, θ1 = R(u, T(2))
    t = mod2piF(θ + θ1)
    v = mod2piF(t - tt)
    cnew = t + u + v
    c <= cnew && return false, c, ctrl
    ctrl = SVector(
        carsegment2stepcontrol(1, t),
        carsegment2stepcontrol(0, u),
        carsegment2stepcontrol(-1, v),
        zero(VelocityCurvatureStep{T}),
        zero(VelocityCurvatureStep{T})
    )
    true, cnew, ctrl
end

@inline function LpRmLp(tgt::SE2State{T}, c::T, ctrl::SVector{5,VelocityCurvatureStep{T}}) where {T}
    tx, ty, tt = tgt.x, tgt.y, tgt.θ
    E = tx - sin(tt)
    N = ty + cos(tt) - 1
    E*E + N*N > 16 && return false, c, ctrl
    r, θ = R(E, N)
    u = acos(1 - r*r/8)
    t = mod2piF(θ - u/2 + pi)
    v = mod2piF(pi - u/2 - θ + tt)
    u = -u
    cnew = t - u + v
    c <= cnew && return false, c, ctrl
    ctrl = SVector(
        carsegment2stepcontrol(1, t),
        carsegment2stepcontrol(-1, u),
        carsegment2stepcontrol(1, v),
        zero(VelocityCurvatureStep{T}),
        zero(VelocityCurvatureStep{T})
    )
    true, cnew, ctrl
end

@inline function LpRmLm(tgt::SE2State{T}, c::T, ctrl::SVector{5,VelocityCurvatureStep{T}}) where {T}
    tx, ty, tt = tgt.x, tgt.y, tgt.θ
    E = tx - sin(tt)
    N = ty + cos(tt) - 1
    E*E + N*N > 16 && return false, c, ctrl
    r, θ = R(E, N)
    u = acos(1 - r*r/8)
    t = mod2piF(θ - u/2 + pi)
    v = mod2piF(pi - u/2 - θ + tt) - 2*T(pi)
    u = -u
    cnew = t - u - v
    c <= cnew && return false, c, ctrl
    ctrl = SVector(
        carsegment2stepcontrol(1, t),
        carsegment2stepcontrol(-1, u),
        carsegment2stepcontrol(1, v),
        zero(VelocityCurvatureStep{T}),
        zero(VelocityCurvatureStep{T})
    )
    true, cnew, ctrl
end

@inline function LpRpuLmuRm(tgt::SE2State{T}, c::T, ctrl::SVector{5,VelocityCurvatureStep{T}}) where {T}
    tx, ty, tt = tgt.x, tgt.y, tgt.θ
    E = tx + sin(tt)
    N = ty - cos(tt) - 1
    p = (2 + sqrt(E*E + N*N)) / 4
    (p < 0 || p > 1) && return false, c, ctrl
    u = acos(p)
    t = mod2piF(Tau(u, -u, E, N))
    v = mod2piF(Omega(u, -u, E, N, tt)) - 2*T(pi)
    cnew = t + 2*u - v
    c <= cnew && return false, c, ctrl
    ctrl = SVector(
        carsegment2stepcontrol(1, t),
        carsegment2stepcontrol(-1, u),
        carsegment2stepcontrol(1, -u),
        carsegment2stepcontrol(-1, v),
        zero(VelocityCurvatureStep{T})
    )
    true, cnew, ctrl
end

@inline function LpRmuLmuRp(tgt::SE2State{T}, c::T, ctrl::SVector{5,VelocityCurvatureStep{T}}) where {T}
    tx, ty, tt = tgt.x, tgt.y, tgt.θ
    E = tx + sin(tt)
    N = ty - cos(tt) - 1
    p = (20 - E*E - N*N) / 16
    (p < 0 || p > 1) && return false, c, ctrl
    u = -acos(p)
    t = mod2piF(Tau(u, u, E, N))
    v = mod2piF(Omega(u, u, E, N, tt))
    cnew = t - 2*u + v
    c <= cnew && return false, c, ctrl
    ctrl = SVector(
        carsegment2stepcontrol(1, t),
        carsegment2stepcontrol(-1, u),
        carsegment2stepcontrol(1, u),
        carsegment2stepcontrol(-1, v),
        zero(VelocityCurvatureStep{T})
    )
    true, cnew, ctrl
end

@inline function LpRmSmLm(tgt::SE2State{T}, c::T, ctrl::SVector{5,VelocityCurvatureStep{T}}) where {T}
    tx, ty, tt = tgt.x, tgt.y, tgt.θ
    E = tx - sin(tt)
    N = ty + cos(tt) - 1
    D, β = R(E, N)
    D < 2 && return false, c, ctrl
    γ = acos(2/D)
    F = sqrt(D*D/4 - 1)
    t = mod2piF(pi + β - γ)
    u = 2 - 2*F
    u > 0 && return false, c, ctrl
    v = mod2piF(-3*T(pi)/2 + γ + tt - β) - 2*T(pi)
    cnew = t + T(pi)/2 - u - v
    c <= cnew && return false, c, ctrl
    ctrl = SVector(
        carsegment2stepcontrol(1, t),
        carsegment2stepcontrol(-1, -T(pi)/2),
        carsegment2stepcontrol(0, u),
        carsegment2stepcontrol(1, v),
        zero(VelocityCurvatureStep{T})
    )
    true, cnew, ctrl
end

@inline function LpRmSmRm(tgt::SE2State{T}, c::T, ctrl::SVector{5,VelocityCurvatureStep{T}}) where {T}
    tx, ty, tt = tgt.x, tgt.y, tgt.θ
    E = tx + sin(tt)
    N = ty - cos(tt) - 1
    D, β = R(E, N)
    D < 2 && return false, c, ctrl
    t = mod2piF(β + T(pi)/2)
    u = 2 - D
    u > 0 && return false, c, ctrl
    v = mod2piF(-pi - tt + β) - 2*T(pi)
    cnew = t + T(pi)/2 - u - v
    c <= cnew && return false, c, ctrl
    ctrl = SVector(
        carsegment2stepcontrol(1, t),
        carsegment2stepcontrol(-1, -T(pi)/2),
        carsegment2stepcontrol(0, u),
        carsegment2stepcontrol(-1, v),
        zero(VelocityCurvatureStep{T})
    )
    true, cnew, ctrl
end

@inline function LpRmSmLmRp(tgt::SE2State{T}, c::T, ctrl::SVector{5,VelocityCurvatureStep{T}}) where {T}
    tx, ty, tt = tgt.x, tgt.y, tgt.θ
    E = tx + sin(tt)
    N = ty - cos(tt) - 1
    D, β = R(E, N)
    D < 2 && return false, c, ctrl
    γ = acos(2/D)
    F = sqrt(D*D/4 - 1)
    t = mod2piF(pi + β - γ)
    u = 4 - 2*F
    u > 0 && return false, c, ctrl
    v = mod2piF(pi + β - tt - γ)
    cnew = t + pi - u + v
    c <= cnew && return false, c, ctrl
    ctrl = SVector(
        carsegment2stepcontrol(1, t),
        carsegment2stepcontrol(-1, -T(pi)/2),
        carsegment2stepcontrol(0, u),
        carsegment2stepcontrol(1, -T(pi)/2),
        carsegment2stepcontrol(-1, v)
    )
    true, cnew, ctrl
end