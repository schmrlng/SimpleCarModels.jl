export dubins, dubins_length, dubins_waypoints

dubins_length(q1::SE2State{T}, q2::SE2State{T}, r::T = T(1)) where {T} = dubins(q1, q2, r)[1]
function dubins_waypoints(q1::SE2State{T}, q2::SE2State{T}, dt_or_N::Union{T,Int}, r::T = T(1), v::T = T(1)) where {T}
    ctrl = dubins(q1, q2, r, v)[2]
    waypoints(SimpleCarDynamics{0,0,T}(), q1, ctrl, dt_or_N)
end

@inline scaleradius(c::VelocityCurvatureStep{T}, α::T) where {T} = StepControl(c.t*α, VelocityCurvatureControl(c.u.v, c.u.κ/α))
@inline scalespeed(c::VelocityCurvatureStep{T}, λ::T) where {T} = StepControl(c.t/λ, VelocityCurvatureControl(c.u.v*λ, c.u.κ))
@inline carsegment2stepcontrol(t::Int, d::T) where {T} = StepControl(abs(d), VelocityCurvatureControl(T(sign(d)), T(t)))

function dubins(q1::SE2State{T},
                q2::SE2State{T},
                r::T = T(1),
                v::T = T(1)) where {T<:AbstractFloat}
    dx = (q2.x - q1.x) / r
    dy = (q2.y - q1.y) / r
    d = sqrt(abs2(dx) + abs2(dy))
    th = atan2(dy, dx)
    a = q1.θ - th
    b = q2.θ - th
    sa, ca = sincos(a)
    sb, cb = sincos(b)
    cmin = T(Inf)
    ctrl = zeros(SVector{3,VelocityCurvatureStep{T}})

    # LSL
    tmp = 2 + d*d - 2*(ca*cb + sa*sb - d*(sa - sb))
    if tmp >= 0
        th = atan2(cb - ca, d + sa - sb)
        t = mod2piF(-a + th)
        p = sqrt(max(tmp, T(0)))
        q = mod2piF(b - th)
        c = t + p + q
        if c < cmin
            cmin = c
            ctrl = SVector(
                carsegment2stepcontrol(1, t),
                carsegment2stepcontrol(0, p),
                carsegment2stepcontrol(1, q)
            )
        end
    end

    # RSR
    tmp = 2 + d*d - 2*(ca*cb + sa*sb - d*(sb - sa))
    if tmp >= 0
        th = atan2(ca - cb, d - sa + sb)
        t = mod2piF(a - th)
        p = sqrt(max(tmp, T(0)))
        q = mod2piF(-b + th)
        c = t + p + q
        if c < cmin
            cmin = c
            ctrl = SVector(
                carsegment2stepcontrol(-1, t),
                carsegment2stepcontrol(0, p),
                carsegment2stepcontrol(-1, q)
            )
        end
    end

    # RSL
    tmp = d*d - 2 + 2*(ca*cb + sa*sb - d*(sa + sb))
    if tmp >= 0
        p = sqrt(max(tmp, T(0)))
        th = atan2(ca + cb, d - sa - sb) - atan2(T(2), p)
        t = mod2piF(a - th)
        q = mod2piF(b - th)
        c = t + p + q
        if c < cmin
            cmin = c
            ctrl = SVector(
                carsegment2stepcontrol(-1, t),
                carsegment2stepcontrol(0, p),
                carsegment2stepcontrol(1, q)
            )
        end
    end

    # LSR
    tmp = -2 + d*d + 2*(ca*cb + sa*sb + d*(sa + sb))
    if tmp >= 0
        p = sqrt(max(tmp, T(0)))
        th = atan2(-ca - cb, d + sa + sb) - atan2(-T(2), p)
        t = mod2piF(-a + th)
        q = mod2piF(-b + th)
        c = t + p + q
        if c < cmin
            cmin = c
            ctrl = SVector(
                carsegment2stepcontrol(1, t),
                carsegment2stepcontrol(0, p),
                carsegment2stepcontrol(-1, q)
            )
        end
    end

    # RLR
    tmp = (6 - d*d  + 2*(ca*cb + sa*sb + d*(sa - sb))) / 8
    if abs(tmp) < 1
        p = 2*T(pi) - acos(tmp)
        th = atan2(ca - cb, d - sa + sb)
        t = mod2piF(a - th + p/2)
        q = mod2piF(a - b - t + p)
        c = t + p + q
        if c < cmin
            cmin = c
            ctrl = SVector(
                carsegment2stepcontrol(-1, t),
                carsegment2stepcontrol(1, p),
                carsegment2stepcontrol(-1, q)
            )
        end
    end

    # LRL
    tmp = (6 - d*d  + 2*(ca*cb + sa*sb - d*(sa - sb))) / 8
    if abs(tmp) < 1
        p = 2*T(pi) - acos(tmp)
        th = atan2(-ca + cb, d + sa - sb)
        t = mod2piF(-a + th + p/2)
        q = mod2piF(b - a - t + p)
        c = t + p + q
        if c < cmin
            cmin = c
            ctrl = SVector(
                carsegment2stepcontrol(1, t),
                carsegment2stepcontrol(-1, p),
                carsegment2stepcontrol(1, q)
            )
        end
    end

    ctrl = scalespeed.(scaleradius.(ctrl, r), v)
    cmin*r, ctrl
end