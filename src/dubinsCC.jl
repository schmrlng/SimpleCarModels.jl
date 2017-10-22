export dubinsCC, dubinsCC_length, dubinsCC_waypoints

dubinsCC_length(q1::SE2κState{T}, q2::SE2κState{T}, κ_max::T = T(1), σ_max::T = T(1)) where {T} = dubinsCC(q1, q2, κ_max, σ_max)[1]
function dubinsCC_waypoints(q1::SE2κState{T}, q2::SE2κState{T}, dt_or_N::Union{T,Int},
                            κ_max::T = T(1), σ_max::T = T(1), v::T = T(1)) where {T}
    ctrl = dubinsCC(q1, q2, κ_max, σ_max, v)[2]
    waypoints(SimpleCarDynamics{0,1,T}(), q1, ctrl, dt_or_N)
end

dubinsCC_length(q1::SE2State{T}, q2::SE2State{T}, κ_max::T = T(1), σ_max::T = T(1)) where {T} = dubinsCC(q1, q2, κ_max, σ_max)[1]
function dubinsCC_waypoints(q1::SE2State{T}, q2::SE2State{T}, dt_or_N::Union{T,Int},
                            κ_max::T = T(1), σ_max::T = T(1), v::T = T(1)) where {T}
    ctrl = dubinsCC(q1, q2, κ_max, σ_max, v)[2]
    waypoints(SimpleCarDynamics{0,1,T}(), SE2κState(q1), ctrl, dt_or_N)
end

@inline scalespeed(c::VelocityCurvRateStep{T}, λ::T) where {T} = StepControl(c.t/λ, c.u*λ)
@inline function θrγ(κ_max::T, σ_max::T) where {T}
    xi = sqrt(T(pi)/σ_max)*fresnelC(κ_max/sqrt(T(pi)*σ_max))
    yi = sqrt(T(pi)/σ_max)*fresnelS(κ_max/sqrt(T(pi)*σ_max))
    θi = κ_max*κ_max/(2*σ_max)
    s, c = sincos(θi)
    xΩ = xi - s/κ_max
    yΩ = yi + c/κ_max
    2*θi, sqrt(xΩ*xΩ + yΩ*yΩ), atan(xΩ/yΩ)
end
SE2State(q::SE2κState) = SE2State(q.x, q.y, q.θ)
SE2κState(q::SE2State{T}) where {T} = SE2κState(q.x, q.y, q.θ, T(0))

# using Plots

function dubinsCC(q1::SE2κState{T},
                  q2::SE2κState{T},
                  κ_max::T = T(1),
                  σ_max::T = T(1),
                  v::T = T(1),
                  θrγ = θrγ(κ_max, σ_max/v)) where {T<:AbstractFloat}
    if abs(q1.κ) > κ_max || abs(q2.κ) > κ_max
        error("$q1, $q2 endpoint curvatures exceed bound of $κ_max")
    end

    D = SimpleCarDynamics{0,1,T}()
    cmin = T(Inf)
    ctrl = zeros(SVector{11,VelocityCurvRateStep{T}})

    t1 = abs(q1.κ)/σ_max
    t2 = abs(q2.κ)/σ_max
    c1_towards = StepControl(t1, VelocityCurvRateControl(v, -flipsign(σ_max, q1.κ)))
    c1_away    = StepControl(-t1, VelocityCurvRateControl(v, flipsign(σ_max, q1.κ)))
    c2_towards = StepControl(t2, VelocityCurvRateControl(v, -flipsign(σ_max, q2.κ)))
    c2_away    = StepControl(-t2, VelocityCurvRateControl(v, flipsign(σ_max, q2.κ)))

    # waypoints(D, q1, c1_towards, 10) |> wps -> plot!([q.x for q in wps], [q.y for q in wps], color=:blue)
    # q1x = propagate(D, q1, c1_away)
    # waypoints(D, q1x, StepControl(t1, c1_away.u), 10)    |> wps -> plot!([q.x for q in wps], [q.y for q in wps], color=:red)
    #
    # waypoints(D, q2, c2_towards, 10) |> wps -> plot!([q.x for q in wps], [q.y for q in wps], color=:blue)
    # q2x = propagate(D, q2, c2_away)
    # waypoints(D, q2x, StepControl(t2, c2_away.u), 10)    |> wps -> plot!([q.x for q in wps], [q.y for q in wps], color=:red)

    # towards 0, towards 0
    q1_0 = SE2State(propagate(D, q1, c1_towards))
    q2_0 = SE2State(propagate(D, q2, c2_towards))
    cnew, ctrl_0 = dubinsCC(q1_0, q2_0, κ_max, σ_max, v, θrγ, true, false, true, true, q2.κ >= 0, q2.κ <= 0)
    cnew = cnew + v*(t1 - t2)
    if cnew < cmin
        cmin = cnew
        ctrl_0 = setindex(ctrl_0, StepControl(ctrl_0[9].t - t2, ctrl_0[9].u), 9)
        ctrl = [SVector(c1_towards); ctrl_0; SVector(zero(VelocityCurvRateStep{T}))]
        # q1_best, q2_best = q1_0, q2_0
    end

    # towards 0, away from 0
    q1_0 = SE2State(propagate(D, q1, c1_towards))
    q2_0 = SE2State(propagate(D, q2, c2_away))
    cnew, ctrl_0 = dubinsCC(q1_0, q2_0, κ_max, σ_max, v, θrγ, true, true, true, true, true, true)
    cnew = cnew + v*(t1 + t2)
    if cnew < cmin
        cmin = cnew
        ctrl = [SVector(c1_towards); ctrl_0; SVector(StepControl(t2, c2_away.u))]
        # q1_best, q2_best = q1_0, q2_0
    end

    # away from 0, towards 0
    q1_0 = SE2State(propagate(D, q1, c1_away))
    q2_0 = SE2State(propagate(D, q2, c2_towards))
    cnew, ctrl_0 = dubinsCC(q1_0, q2_0, κ_max, σ_max, v, θrγ, false, false, q1.κ >= 0, q1.κ <= 0, q2.κ >= 0, q2.κ <= 0)
    cnew = cnew + v*(-t1 - t2)
    if cnew < cmin
        cmin = cnew
        ctrl_0 = setindex(ctrl_0, StepControl(ctrl_0[1].t - t1, ctrl_0[1].u), 1)
        ctrl_0 = setindex(ctrl_0, StepControl(ctrl_0[9].t - t2, ctrl_0[9].u), 9)
        ctrl = [SVector(zero(VelocityCurvRateStep{T})); ctrl_0; SVector(zero(VelocityCurvRateStep{T}))]
        # q1_best, q2_best = q1_0, q2_0
    end

    # away from 0, away from 0
    q1_0 = SE2State(propagate(D, q1, c1_away))
    q2_0 = SE2State(propagate(D, q2, c2_away))
    cnew, ctrl_0 = dubinsCC(q1_0, q2_0, κ_max, σ_max, v, θrγ, false, true, q1.κ >= 0, q1.κ <= 0, true, true)
    cnew = cnew + v*(-t1 + t2)
    if cnew < cmin
        cmin = cnew
        ctrl_0 = setindex(ctrl_0, StepControl(ctrl_0[1].t - t1, ctrl_0[1].u), 1)
        ctrl = [SVector(zero(VelocityCurvRateStep{T})); ctrl_0; SVector(StepControl(t2, c2_away.u))]
        # q1_best, q2_best = q1_0, q2_0
    end

    # scatter!([q1_best.x, q2_best.x], [q1_best.y, q2_best.y], color=:green)

    cmin, ctrl
end

function dubinsCC(q1::SE2State{T},
                  q2::SE2State{T},
                  κ_max::T = T(1),
                  σ_max::T = T(1),
                  v::T = T(1),
                  θrγ = θrγ(κ_max, σ_max/v),
                  allow_short_turn_1 = true,
                  allow_short_turn_3 = true,
                  LXX = true,
                  RXX = true,
                  XXL = true,
                  XXR = true) where {T<:AbstractFloat}
    θ_lim, r, γ = θrγ
    @inline turn_length(β) = CC_turn_length(β, κ_max, σ_max/v, θ_lim, r, γ)
    @inline turn_control(β) = CC_turn_control(β, κ_max, σ_max/v, θ_lim, r, γ)

    dx = (q2.x - q1.x) / r
    dy = (q2.y - q1.y) / r
    d = sqrt(abs2(dx) + abs2(dy))
    th = atan2(dy, dx)
    a = q1.θ - th
    b = q2.θ - th

    sap, cap = sincos(a+γ)
    sbp, cbp = sincos(b+γ)
    sam, cam = sincos(a-γ)
    sbm, cbm = sincos(b-γ)
    sγ, cγ = sincos(γ)
    cmin = T(Inf)
    ctrl = zeros(SVector{9,VelocityCurvRateStep{T}})

    ### 1. LSL: a-γ, b+γ
    if LXX && XXL
        ca, sa, cb, sb = cam, sam, cbp, sbp
        tmp = 2 + d*d - 2*(ca*cb + sa*sb - d*(sa - sb))
        if tmp >= 4*sγ*sγ    # TODO: M-style connections
            th = atan2(cb - ca, d + sa - sb)
            t = mod2piF(-a + th)
            p = sqrt(tmp) - 2*sγ
            q = mod2piF(b - th)
            p = r*p
            if (allow_short_turn_1 || t >= θ_lim) && (allow_short_turn_3 || q >= θ_lim)
                c = turn_length(t) + p + turn_length(q)
                if c < cmin
                    cmin = c
                    ctrl = [
                        turn_control(t);
                        SVector(zero(VelocityCurvRateStep{T}),
                                StepControl(p, VelocityCurvRateControl(T(1), T(0))),
                                zero(VelocityCurvRateStep{T}));
                        turn_control(q)
                    ]
                end
            end
        end
    end

    ### 2. RSR: a+γ, b-γ
    if RXX && XXR
        ca, sa, cb, sb = cap, sap, cbm, sbm
        tmp = 2 + d*d - 2*(ca*cb + sa*sb - d*(sb - sa))
        if tmp >= 4*sγ*sγ    # TODO: M-style connections
            th = atan2(ca - cb, d - sa + sb)
            t = mod2piF(a - th)
            p = sqrt(tmp) - 2*sγ
            q = mod2piF(-b + th)
            p = r*p
            if (allow_short_turn_1 || t >= θ_lim) && (allow_short_turn_3 || q >= θ_lim)
                c = turn_length(t) + p + turn_length(q)
                if c < cmin
                    cmin = c
                    ctrl = [
                        turn_control(-t);
                        SVector(zero(VelocityCurvRateStep{T}),
                                StepControl(p, VelocityCurvRateControl(T(1), T(0))),
                                zero(VelocityCurvRateStep{T}));
                        turn_control(-q)
                    ]
                end
            end
        end
    end

    ### 3. RSL: a+γ, b+γ
    if RXX && XXL
        ca, sa, cb, sb = cap, sap, cbp, sbp
        tmp = d * d - 2 + 2 * (ca*cb + sa*sb - d * (sa + sb))
        if tmp >= 0
            p = sqrt(tmp + 4*sγ*sγ) - 2*sγ
            x = sqrt(tmp)
            th = atan2(ca + cb, d - sa - sb) - atan2(T(2), x)
            ε = acot(x/2) - acot((p/2 + sγ)/cγ)
            t = mod2piF(a - th - ε)
            q = mod2piF(b - th - ε)
            p = r*p
            if (allow_short_turn_1 || t >= θ_lim) && (allow_short_turn_3 || q >= θ_lim)
                c = turn_length(t) + p + turn_length(q)
                if c < cmin
                    cmin = c
                    ctrl = [
                        turn_control(-t);
                        SVector(zero(VelocityCurvRateStep{T}),
                                StepControl(p, VelocityCurvRateControl(T(1), T(0))),
                                zero(VelocityCurvRateStep{T}));
                        turn_control(q)
                    ]
                end
            end
        end
    end

    ### 4. LSR: a-γ, b-γ
    if LXX && XXR
        ca, sa, cb, sb = cam, sam, cbm, sbm
        tmp = -2 + d * d + 2 * (ca*cb + sa*sb + d * (sa + sb))
        if tmp >= 0
            p = sqrt(tmp + 4*sγ*sγ) - 2*sγ
            x = sqrt(tmp)
            th = atan2(-ca - cb, d + sa + sb) - atan2(-T(2), x)
            ε = acot(x/2) - acot((p/2 + sγ)/cγ)
            t = mod2piF(-a + th - ε)
            q = mod2piF(-b + th - ε)
            p = r*p
            if (allow_short_turn_1 || t >= θ_lim) && (allow_short_turn_3 || q >= θ_lim)
                c = turn_length(t) + p + turn_length(q)
                if c < cmin
                    cmin = c
                    ctrl = [
                        turn_control(t);
                        SVector(zero(VelocityCurvRateStep{T}),
                                StepControl(p, VelocityCurvRateControl(T(1), T(0))),
                                zero(VelocityCurvRateStep{T}));
                        turn_control(-q)
                    ]
                end
            end
        end
    end

    ### 5. RLR: a+γ, b-γ
    if RXX && XXR
        ca, sa, cb, sb = cap, sap, cbm, sbm
        tmp = (6 - d * d + 2 * (ca*cb + sa*sb + d * (sa - sb))) / 8
        if abs(tmp) < 1
            p_1 = mod2piF(acos(tmp) - 2*γ)
            p_2 = mod2piF(-p_1 - 4*γ)
            t_2 = mod2piF(a - atan2(ca - cb, d - sa + sb) + p_2/2)
            q_2 = mod2piF(a - b - t_2 + p_2)
            t_1 = mod2piF(t_2 + T(pi) - p_2 - 2*γ)
            q_1 = mod2piF(q_2 + T(pi) - p_2 - 2*γ)
            if (allow_short_turn_1 || t_1 >= θ_lim) && (allow_short_turn_3 || q_1 >= θ_lim)
                c_1 = turn_length(t_1) + turn_length(p_1) + turn_length(q_1)
                if c_1 < cmin
                    cmin = c_1
                    ctrl = [
                        turn_control(-t_1);
                        turn_control(p_1);
                        turn_control(-q_1)
                    ]
                end
            end
            if (allow_short_turn_1 || t_2 >= θ_lim) && (allow_short_turn_3 || q_2 >= θ_lim)
                c_2 = turn_length(t_2) + turn_length(p_2) + turn_length(q_2)
                if c_2 < cmin
                    cmin = c_2
                    ctrl = [
                        turn_control(-t_2);
                        turn_control(p_2);
                        turn_control(-q_2)
                    ]
                end
            end
        end
    end

    ### 6. LRL: a-γ, b+γ
    if LXX && XXL
        ca, sa, cb, sb = cam, sam, cbp, sbp
        tmp = (6 - d * d + 2 * (ca*cb + sa*sb - d * (sa - sb))) / 8
        if abs(tmp) < 1
            p_1 = mod2piF(acos(tmp) - 2*γ)
            p_2 = mod2piF(-p_1 - 4*γ)
            t_2 = mod2piF(-a + atan2(-ca + cb, d + sa - sb) + p_2/2)
            q_2 = mod2piF(b - a - t_2 + p_2)
            t_1 = mod2piF(t_2 + T(pi) - p_2 - 2*γ)
            q_1 = mod2piF(q_2 + T(pi) - p_2 - 2*γ)
            if (allow_short_turn_1 || t_1 >= θ_lim) && (allow_short_turn_3 || q_1 >= θ_lim)
                c_1 = turn_length(t_1) + turn_length(p_1) + turn_length(q_1)
                if c_1 < cmin
                    cmin = c_1
                    ctrl = [
                        turn_control(t_1);
                        turn_control(-p_1);
                        turn_control(q_1)
                    ]
                end
            end
            if (allow_short_turn_1 || t_2 >= θ_lim) && (allow_short_turn_3 || q_2 >= θ_lim)
                c_2 = turn_length(t_2) + turn_length(p_2) + turn_length(q_2)
                if c_2 < cmin
                    cmin = c_2
                    ctrl = [
                        turn_control(t_2);
                        turn_control(-p_2);
                        turn_control(q_2)
                    ]
                end
            end
        end
    end

    ctrl = scalespeed.(ctrl, v)
    cmin, ctrl
end

@inline function CC_turn_length(β::T, κ_max::T, σ_max::T, θ_lim::T, r::T, γ::T) where {T}    # θ_lim = κ_max^2/σ_max
    b = abs(β)
    if b < T(1e-5)
        2*r*sin(γ) + b*r*cos(γ)
    elseif b < θ_lim
        bp = sqrt(b/T(pi))
        σ = (cos(b/2)*fresnelC(bp) + sin(b/2)*fresnelS(bp)) / r / sin(b/2 + γ)
        σ = T(pi)*σ*σ
        2*sqrt(b/σ)
    else
        (b+θ_lim)/κ_max
    end
end

@inline function CC_turn_control(β::T, κ_max::T, σ_max::T, θ_lim::T, r::T, γ::T) where {T}    # θ_lim = κ_max^2/σ_max
    b = abs(β)
    if b < T(1e-5)
        sγ, cγ = sincos(γ)
        SVector(
            StepControl(r*sγ + b*r*cγ/2, VelocityCurvRateControl(T(1), b/(r*r*sγ*sγ))),
            zero(VelocityCurvRateStep{T}),
            StepControl(r*sγ + b*r*cγ/2, VelocityCurvRateControl(T(1), -b/(r*r*sγ*sγ)))
        )
    elseif b < θ_lim
        bp = sqrt(b/T(pi))
        σ = (cos(b/2)*fresnelC(bp) + sin(b/2)*fresnelS(bp)) / r / sin(b/2 + γ)
        σ = T(pi)*σ*σ
        SVector(
            StepControl(sqrt(b/σ), VelocityCurvRateControl(T(1), flipsign(σ, β))),
            zero(VelocityCurvRateStep{T}),
            StepControl(sqrt(b/σ), VelocityCurvRateControl(T(1), -flipsign(σ, β)))
        )
    else
        SVector(
            StepControl(κ_max/σ_max, VelocityCurvRateControl(T(1), flipsign(σ_max, β))),
            StepControl((b-θ_lim)/κ_max, VelocityCurvRateControl(T(1), T(0))),
            StepControl(κ_max/σ_max, VelocityCurvRateControl(T(1), -flipsign(σ_max, β)))
        )
    end
end