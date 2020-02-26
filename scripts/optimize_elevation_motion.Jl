@everywhere using Healpix, Stripeline, Plots, Dates, Statistics, Distributed
@everywhere using SharedArrays, FITSIO
@everywhere const Sl = Stripeline

pars_arg = ARGS

@everywhere get_fsky(map) = count(map.pixels .> 0) / length(map.pixels) * 100

@everywhere function get_χ²(hitmap, mean_hits)
    vals = filter(x -> x.>0, hitmap.pixels)
    χ² = 0.0
    for hit in vals
        χ² += (hit - mean_hits)^2 / mean_hits
    end
    χ²
end

@everywhere function hit_requirement(;
                         k=1/sqrt(2), T_cmb=1.822, T_atm=9, T_noise=26.5,
                         sampling_rate=50, sensitivity=1.12e-6, band=7.74e9, l=1,
                         d=0.35)
    """
    It is the radiometer equation expressed as a function of:
        - T_sys = T_cmb + T_atm + T_noise
        - sampling_rate
        - sensitivity requirement
        - squared pixels of side l [deg] (area [rad**2])
        - duty cycle
    """
    N = deg2rad(l)^2 / Healpix.nside2pixarea(256)
    return k^2 * (T_cmb + T_atm + T_noise)^2 * sampling_rate / sensitivity^2 /
        band / N / d
end

@everywhere function mean_hit(map)
    it, c = 0, 0
    for p in map.pixels
        if p > 0.0
            it += p
            c += 1
        end
    end
    return it / c
end

@everywhere function bezier1D(t, elevation_angle_deg, tilt_angle_deg, T, a, b) 
    t /= T
    k = t - floor(t) 
    y = 3a * k * (1 - k)^2 + 3b * k^2 * (1- k) + k^3
    return elevation_angle_deg + tilt_angle_deg/2 * sin(2π*y)
end

@everywhere function modulate_sin(t, elevation_angle_deg, tilt_angle_deg, T)
    elevation_angle_deg + tilt_angle_deg/2 * (sin(2π/T*t) + sin(π/T*t)) 
end

@everywhere function modulate_evenly_elevation(t, elevation_angle_deg,
                                               tilt_angle_deg, T)
    elevation_angle_deg + tilt_angle_deg/2 * sin(2π/T*t)
end

@everywhere function modulate_unevenly_elevation(t, elevation_angle_deg,
                                                 tilt_angle_deg, T;
                                     a=0.45, b=2*0.46, ϕ=-30.0)
    elevation_angle_deg += 5
    tilt_angle_deg += 16
    o = elevation_angle_deg + tilt_angle_deg/2 * sin(2π/T)
    w = o^-(2/3) / ((elevation_angle_deg + tilt_angle_deg/2)^-(2/3))
    offset = o*w - elevation_angle_deg
    el = elevation_angle_deg + tilt_angle_deg/2 * sin(2π*t/T)
    weight = el^-(2/3) / ((elevation_angle_deg + tilt_angle_deg/2)^-(2/3))
    el * weight - offset*(a*cos.(2π/T*t + ϕ) + b)
end

@everywhere function get_hits_map(elevation_angle::Float64,
                      spin_velocity::Float64,
                      sampling_rate::Float64,
                      time_duration::Int,
                      t_start::DateTime,
                      NSIDE::Int,
                      latitude::Float64,
                      longitude::Float64,
                      height::Int,
                      pointings::Array{Float64,1},
                      effects::Dict,
                      modulation::Function;
                      tilt_angle::Float64=0.0,
                      T::Float64=36000.0,
                      a::Float64=0.1,
                      b::Float64=0.1,
                      chunks_len::Int=86400)

    τ_s = 1 / sampling_rate
    times = 0:τ_s:time_duration
    
    skydirs = Array{Float64}(undef, length(times), 2)
    skyψ = Array{Float64}(undef, length(times), 1)
    hits_map = Healpix.Map{Float64, Healpix.RingOrder}(NSIDE)
    
    (skydirs, skyψ) = Sl.genpointings(pointings,
                                      times,
                                      t_start;
                                      latitude_deg=latitude,
                                      longitude_deg=longitude,
                                      height_m=height,
                                      precession=effects["precess"],
                                      nutation=effects["nutate"],
                                      aberration=effects["aberrate"],
                                      refract=effects["refract"]) do t_s
                                      return (
                                          0.0,
                                          deg2rad(modulation(t_s,
                                                             elevation_angle,
                                                             tilt_angle,
                                                             T,
                                                             a,
                                                             b)),
                                          Sl.timetorotang(t_s,
                                                          spin_velocity))
                                      end
        
    for i in eachindex(skyψ)
        hits_map.pixels[Healpix.ang2pix(hits_map, π/2 - skydirs[i, 1],
                                        skydirs[i, 2])] += 1.0
    end

    hits_map
end

@everywhere function get_optimizer(A, B;
                       NSIDE=256,
                       elevation_angle_deg=20.0,
                       modulation_period=600.0,
                       tilt_angle_deg=30.0,
                       sampling_rate_Hz=1.0,
                       t_start=DateTime(2021, 4, 1, 0, 0, 0),
                       t_stop=DateTime(2021, 4, 2, 0, 0, 0),
                       only_I0=false)

    effects = Dict("precess" => true,
                   "nutate" => true,
                   "aberrate" => true,
                   "refract" => true)
    
    time_duration_s = ((t_stop - t_start) / 1000).value 
    
    tenerife_latitude_deg = Sl.TENERIFE_LATITUDE_DEG 
    tenerife_longitude_deg = Sl.TENERIFE_LONGITUDE_DEG
    tenerife_height_m = Sl.TENERIFE_HEIGHT_M 
    
    modulation = bezier1D
    
    spin_velocity_rpm = 1.0
    
    db = Sl.InstrumentDB()
        
    if only_I0
        keys_dbfp = ["I0"]
    else
        keys_dbfp = collect(keys(db.focalplane))
        for (j, f) in enumerate(keys_dbfp)
            if startswith(f, "W")
                delete!(db.focalplane, f)
                deleteat!(keys_dbfp, j)
            end
        end
    end
        
    hitmap = Healpix.Map{Float64, Healpix.RingOrder}(NSIDE)
    hitmap.pixels = zeros(12NSIDE^2)
    for k in keys_dbfp
        oned_map = get_hits_map(elevation_angle_deg,
                                spin_velocity_rpm,
                                sampling_rate_Hz,
                                time_duration_s,
                                t_start,
                                NSIDE,
                                tenerife_latitude_deg,
                                tenerife_longitude_deg,
                                tenerife_height_m,
                                db.focalplane[k].orientation,    
                                effects,
                                modulation,
                                T=modulation_period,
                                tilt_angle=tilt_angle_deg,
                                a=A,
                                b=B)
        hitmap.pixels += oned_map.pixels
    end
    mean_hits = mean_hit(hitmap)
    return get_χ²(hitmap, mean_hits)
end


nside = 256
elevation = 20.0
T = 3600.0
tilt = 30.0
sampling_rate = 10.0
onlyI0 = false
tstart = DateTime(2021, 4, 1, 0, 0, 0)
tstop = DateTime(2021, 4, 2, 0, 0, 0)

duration = ((tstop - tstart) / 1000).value 

grid_resolution = 0.1

X = 0:grid_resolution:1
Y = 0:grid_resolution:1

χ² = Array{Float64}(undef, length(X), length(Y))
χ²_ = SharedArray{Float64}(length(X)*length(Y))
S = Array{Tuple{Float64, Float64}}(undef, length(X)*length(Y))

for (i, x) in enumerate(X)
    for (j, y) in enumerate(Y)
        S[length(X)*(i-1) + j] = (x, y)
    end
end

@sync @distributed for (i, s) in collect(enumerate(S))
    x = s[1]
    y = s[2]
    χ²_[i] = get_optimizer(x, y;
                           NSIDE=nside,
                           elevation_angle_deg=elevation,
                           modulation_period=T,
                           tilt_angle_deg=tilt,
                           sampling_rate_Hz=sampling_rate,
                           t_start=tstart,
                           t_stop=tstop,
                           only_I0=onlyI0)
end

for i = 1:length(X)
    for j = 1:length(Y)
        χ²[i, j] = χ²_[length(X)*(i-1) + j]
    end
end

f = FITS(join([pars_arg[1], ".fits"]), "w")

write(f, χ²)

write_key(f[1], "NSIDE", nside, "Healpix nside")
write_key(f[1], "ANG_EL", elevation, "Elevation angle [Deg]")
write_key(f[1], "TILT_ANG", tilt, "Tilt angle [Deg]")
write_key(f[1], "S_RATE", tilt, "Sampling rate [Hz]")
write_key(f[1], "PERIOD", T, "Period of the modulation [sec]")
write_key(f[1], "GRID_RES", grid_resolution,
          "X any Y range from 0 to 1 with this resolution")

write_key(f[1], "DURATION", duration, "Simulation duration[sec]")

if onlyI0 == true
    write_key(f[1], "DETECTOR", "I0", "Only detector I0 used")
else
    write_key(f[1], "DETECTOR", "all detectors", "Sum of all detectors")
end

close(f)
