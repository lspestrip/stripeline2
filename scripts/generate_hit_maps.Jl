using Stripeline, Healpix, ColorSchemes, Plots, FITSIO, Dates, MPI, JSON
const Sl = Stripeline

MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm) 

function bezier1D(t, elevation_angle_deg, tilt_angle_deg, T; a=0.45, b=0.0) 
    t /= T
    k = t - floor(t) 
    y = 3a * k * (1 - k)^2 + 3b * k^2 * (1- k) + k^3
    return elevation_angle_deg + tilt_angle_deg/2 * sin(2π*y)
end

function modulate_evenly_elevation(t, elevation_angle_deg, tilt_angle_deg, T)
    elevation_angle_deg + tilt_angle_deg/2 * sin(2π/T*t)
end

function modulate_unevenly_elevation(t, elevation_angle_deg, tilt_angle_deg, T;
                                     a=0.45, b=2*0.46, ϕ=-0.30)
    elevation_angle_deg += 5
    tilt_angle_deg += 16
    o = elevation_angle_deg + tilt_angle_deg/2 * sin(2π/T)
    w = o^-(2/3) / ((elevation_angle_deg + tilt_angle_deg/2)^-(2/3))
    offset = o*w - elevation_angle_deg
    el = elevation_angle_deg + tilt_angle_deg/2 * sin(2π*t/T)
    weight = el^-(2/3) / ((elevation_angle_deg + tilt_angle_deg/2)^-(2/3))
    el * weight - offset*(a*cos.(2π/T*t + ϕ) + b)
end

function get_hits_map(elevation_angle::Float64,
                      spin_velocity::Float64,
                      sampling_rate::Float64,
                      time_duration::Int,
                      t_start::DateTime,
                      NSIDE::Int,
                      latitude::Float64,
                      longitude::Float64,
                      height::Int,
                      pointings::Array{Float64,1},
<<<<<<< HEAD
                      effects::Dict,
                      modulation::Function;
=======
                      effects::Dict;
>>>>>>> e4dd169f9952e26b16292dccd44ce64cf69db67e
                      tilt_angle::Float64=0.0,
                      T::Float64=36000.0,
                      chunks_len::Int=86400)

    τ_s = 1 / sampling_rate
    times = 0:τ_s:time_duration
    
    skydirs = Array{Float64}(undef, length(times), 2)
    skyψ = Array{Float64}(undef, length(times), 1)
    hits_map = Healpix.Map{Float64, Healpix.RingOrder}(NSIDE)

    nchunks = round(Int, times.len / chunks_len)

    actual_times = times[1:chunks_len]
    actual_t_start = t_start
    actual_t_stop = actual_t_start + Dates.Second(chunks_len-1)
    
    for n = 1:nchunks
        
        if n == nchunks
            actual_times = times[(n-1)*chunks_len+1:end]
            actual_t_stop = t_stop
        end
        
        (skydirs, skyψ) = Sl.genpointings(pointings,
                                          actual_times,
                                          actual_t_start;
                                          latitude_deg=latitude,
                                          longitude_deg=longitude,
                                          height_m=height,
                                          precession=effects["precess"],
                                          nutation=effects["nutate"],
                                          aberration=effects["aberrate"],
<<<<<<< HEAD
                                          refract=effects["refract"]) do t_s
                                          return (
                                              0.0,
                                              deg2rad(modulation(t_s,
                                                                 elevation_angle,
                                                                 tilt_angle,
                                                                 T)),
                                              Sl.timetorotang(t_s,
                                                              spin_velocity))
                                          end
=======
                                          refract=effects["refract"]) do time_s
                                              return (
                                                  0.0,
                                                  deg2rad(el_angle_mod(time_s)),
                                                  Sl.timetorotang(time_s,
                                                                  spin_velocity))
                                      end
>>>>>>> e4dd169f9952e26b16292dccd44ce64cf69db67e
        
        for i in eachindex(skyψ)
            
            hits_map.pixels[Healpix.ang2pix(hits_map, π/2 - skydirs[i, 1],
                                            skydirs[i, 2])] += 1.0
        end
        
        if n != nchunks
            actual_t_start = actual_t_stop + Dates.Second(1)
            actual_t_stop = actual_t_start + Dates.Second(chunks_len-1)
            actual_times = times[n*chunks_len+1:(n+1)*chunks_len]
        end
    end
    hits_map
end

pars_arg = ARGS

NSIDE = 256

effects = Dict("precess" => true,
               "nutate" => true,
               "aberrate" => true,
               "refract" => true)

<<<<<<< HEAD
modulation = bezier1D  # choose the function to modulate the elevation
elevation_angle_deg = 20.0     
modulation_period = parse(Float64, pars_arg[1]) 
tilt_angle_deg = 30.0

    
=======
>>>>>>> e4dd169f9952e26b16292dccd44ce64cf69db67e
day_start = parse(Int, pars_arg[3])
month_start = parse(Int, pars_arg[4])
year_start = parse(Int, pars_arg[5])
day_stop = parse(Int, pars_arg[6])
month_stop = parse(Int, pars_arg[7])
year_stop = parse(Int, pars_arg[8])

t_start = DateTime(year_start, month_start, day_start, 0, 0, 0)
t_stop = DateTime(year_stop, month_stop, day_stop, 0, 0, 0)
time_duration_s = ((t_stop - t_start) / 1000).value 

tenerife_latitude_deg = Sl.TENERIFE_LATITUDE_DEG 
tenerife_longitude_deg = Sl.TENERIFE_LONGITUDE_DEG
tenerife_height_m = Sl.TENERIFE_HEIGHT_M 

spin_velocity_rpm = 1.0
sampling_rate_hz = 50.0

db = Sl.InstrumentDB()
keys_dbfp = collect(keys(db.focalplane))

for (j, f) in enumerate(keys_dbfp)
    if startswith(f, "W")
        delete!(db.focalplane, f)
        deleteat!(keys_dbfp, j)
    end
end

if size > length(keys_dbfp)
    N = Sl.split_into_n(length(keys_dbfp), length(keys_dbfp))
else
    N = Sl.split_into_n(length(keys_dbfp), size)
end

for r in rank
    z = r + 1
    if z > length(keys_dbfp)
    else
        subkeys = keys_dbfp[sum(N[1:z])-N[z]+1:sum(N[1:z])]
        for k in subkeys
            hits_map = get_hits_map(elevation_angle_deg,
                                    spin_velocity_rpm,
                                    sampling_rate_hz,
                                    time_duration_s,
                                    t_start,
                                    NSIDE,
                                    tenerife_latitude_deg,
                                    tenerife_longitude_deg,
                                    tenerife_height_m,
                                    db.focalplane[k].orientation,
<<<<<<< HEAD
                                    effects,
                                    modulation,
                                    T=modulation_period,
                                    tilt_angle=tilt_angle_deg)
=======
                                    effects)
>>>>>>> e4dd169f9952e26b16292dccd44ce64cf69db67e
                        
            save_to_path = join([pars_arg[2],
                                 db.focalplane[k].name,
                                 "_",
                                 pars_arg[1],
                                 "°_",
                                 Dates.format(t_start, "yyyy-mm-dd"),
                                 "_",
                                 Dates.format(t_stop, "yyyy-mm-dd"), 
                                 ".fits"])
            
            Healpix.saveToFITS(hits_map, save_to_path, typechar = "D")
            
            f = FITS(save_to_path, "r+")
                        
            write_key(f[2], "LATITUDE", tenerife_latitude_deg,
                      "Tenerife Latitude [Deg]")
            write_key(f[2], "LONGITUD", tenerife_longitude_deg,
                      "Tenerife Longitude [Deg]")
            write_key(f[2], "HEIGHT", tenerife_height_m, "Tenerife Height [m]")
            
            write_key(f[2], "DURATION", time_duration_s, "[sec]")
            write_key(
                f[2], "TSTART", Dates.format(t_start, "yyyy-mm-dd HH:MM:SS"),
                "Starting Date and time of the observation: yyyy-mm-dd HH:MM:SS")
            write_key(
                f[2], "TSTOP", Dates.format(t_stop, "yyyy-mm-dd HH:MM:SS"),
                "Ending Date and time of the observation: yyyy-mm-dd HH:MM:SS")

            if tilt_angle_deg == 0.0
                write_key(f[2], "ANG_EL", elevation_angle_deg,
                          "Elevation angle [Deg]")
            else
                write_key(f[2], "ANG_EL", string(modulation),
                          "Elevation angle [Deg]")
            end
             
            write_key(f[2], "SPIN_VEL", spin_velocity_rpm, "Spin velocity [rpm]")
            write_key(f[2], "S_RATE", sampling_rate_hz, "Sampling rate [Hz]")
            
            write_key(f[2], "DETECTOR", db.focalplane[k].name, "Detector name")
            
<<<<<<< HEAD
            write_key(f[2], "PRECESS", effects["precess"],
                      "If true precession is taken into account")
            write_key(f[2], "NUTATE", effects["nutate"], "If true nutation is 
                      taken into account")
            write_key(f[2], "ABERRATE", effects["aberrate"],
                      "If true aberration is taken into account")
            write_key(f[2], "REFRACT", effects["refract"],
                      "If true aberration istaken into account")
=======
            write_key(f[2], "PRECESS", effects["precess"], "If true precession is 
                      taken into account")
            write_key(f[2], "NUTATE", effects["nutate"], "If true nutation is 
                      taken into account")
            write_key(f[2], "ABERRATE", effects["aberrate"], "If true aberration is 
                      taken into account")
            write_key(f[2], "REFRACT", effects["refract"], "If true aberration is 
                      taken into account")
>>>>>>> e4dd169f9952e26b16292dccd44ce64cf69db67e

            close(f)
        end
    end
end
    
MPI.Finalize()
