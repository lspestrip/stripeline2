using Stripeline, Healpix, ColorSchemes, Plots, FITSIO
const Sl = Stripeline
const Hp = Healpix

searchdir(path, key) = filter(x->occursin(key, x), readdir(path))

pars_arg = ARGS

elevation_angle_deg = pars_arg[1]
path = pars_arg[2]
name = pars_arg[3]

files = searchdir(path, join([elevation_angle_deg, "°"]))
g = FITS(join([path, files[1]]), "r")
z = FITS(join([path, files[end]]), "r")
h = read_header(g[2])
y = read_header(z[2])

hits_map = Healpix.Map{Float64, Healpix.RingOrder}(h["NSIDE"])
hits_map.pixels = zeros(length(hits_map.pixels))

db = Sl.InstrumentDB()
keys_dbfp = collect(keys(db.focalplane))

for (j, f) in enumerate(keys_dbfp)
    if startswith(f, "W")
        delete!(db.focalplane, f)
        deleteat!(keys_dbfp, j)
    end
end

for file in files
    load_path = join([path, file])
    single_detector_map = Healpix.readMapFromFITS(load_path, 1, Float64)
    hits_map.pixels += single_detector_map.pixels
end

save_to_path = join([path,
                     name,
                     elevation_angle_deg,
                     "°_",
                     h["TSTART"][1:10],
                     "_",
                     h["TSTOP"][1:10],
                     ".fits"])

Healpix.saveToFITS(hits_map, save_to_path, typechar = "D")

f = FITS(save_to_path, "r+")

write_key(f[2], "LATITUDE", h["LATITUDE"], "Tenerife Latitude [Deg]")
write_key(f[2], "LONGITUD", h["LONGITUD"], "Tenerife Longitude [Deg]")
write_key(f[2], "HEIGHT", h["HEIGHT"], "Tenerife Height [m]")

write_key(f[2], "DURATION", h["DURATION"], "[sec]")
write_key(f[2], "TSTART", h["TSTART"],
          "Starting Date and time of the observations")    
write_key(f[2], "TSTOP", y["TSTOP"], "Ending Date and time of the observations")

write_key(f[2], "ANG_EL", h["ANG_EL"], "Elevation angle [Deg]")
write_key(f[2], "SPIN_VEL", h["SPIN_VEL"], "Spin velocity [rpm]")
write_key(f[2], "S_RATE", h["S_RATE"], "Sampling rate [Hz]")

write_key(f[2], "DETECTOR", "all detectors", "Sum of all detectors")

<<<<<<< HEAD
write_key(f[2], "PRECESS", h["PRECESS"],
          "If true precession is taken into account")
write_key(f[2], "NUTATE", h["NUTATE"], "If true nutation is taken into account")
write_key(f[2], "ABERRATE", h["ABERRATE"],
          "If true aberration is taken into account")
write_key(f[2], "REFRACT", h["REFRACT"],
          "If true aberration is taken into account")
=======
write_key(f[2], "PRECESS", h["PRECESS"], "If true precession is taken into account")
write_key(f[2], "NUTATE", h["NUTATE"], "If true nutation is taken into account")
write_key(f[2], "ABERRATE", h["ABERRATE"], "If true aberration is taken into account")
write_key(f[2], "REFRACT", h["REFRACT"], "If true aberration is taken into account")

close(f)
>>>>>>> e4dd169f9952e26b16292dccd44ce64cf69db67e

close(f)
