export fix_fknee, get_info_from_DB
import Random


"""
    function fix_fknee(f_samp_hz, tcmb_k, tatm_k, ttel_k, tnoise_k, β_hz, slope, fknee_bicocca)
    
    This function can be used in order to adapt the fknees measured during Bicocca's unit tests to the nominal conditions 
    of the future data taking in Tenerife.
    The signal "seen" by the polarimeter in Bicocca was about 21.3 K (22.7 on port A, 29.9 on port B of the Magic Tee as written in
    Claudio Pincella's Thesis), while in Tenerife will be of about 28.5 (tatm_k + ttel_k + tcmb_k).
    Moreover, the sampling frequency used in Bicocca was 25 Hz, at Tenerife will be 50 Hz.
    As a consequence, the white noise level will be different in Tenerife, ad so will be the fknee (the 1/f slope, instead is the same).

    This function corrects the fknee according to the following formula:

    fknee_corrected = fknee*(sigma_bicocca/sigma_tenerife)^(1/slope)

    Given:
    - the sampling frequency, in Hz
    - the CMB temperature (in K)
    - the atmospheric temperature, in K
    - the telescope temperature, in K
    - the noise temperature of the polarimeter, in K
    - the bandwidth of the polarimeter, in Hz
    - the 1/f slope of the polarimeter
    - the fknee of the polarimeter, as it is from the DB

    It returns the corrected fknee.
"""
function fix_fknee(f_samp_hz, tcmb_k, tatm_k, ttel_k, tnoise_k, β_hz, slope, fknee_bicocca)

    τ_s_tenerife = 1/f_samp_hz
    τ_s_bicocca = 1/25      

    sigma_bicocca = (1/sqrt(2))*(21.3+tnoise_k) / (sqrt(β_hz * τ_s_bicocca))   
    sigma_tenerife = (1/sqrt(2))*(tnoise_k + tatm_k + ttel_k + tcmb_k) / (sqrt(β_hz * τ_s_tenerife))  
                                                                           
    fknee_tenerife = fknee_bicocca*(sigma_bicocca/sigma_tenerife)^(1/slope)
end


"""
    function get_info_from_DB(db, horns, stokes, f_samp_hz, tcmb_k, tatm_k, ttel_k)
    
This function can be used to take information about specific horns
from the STRIP instrument database.

Given:
- The database (db = Stripeline.InstrumentDB())
- A vector containg the names of the horns of interest (e.g. ["Y0",
  "I1", "R6"])
- A string indicating the Stokes parameter object of the simulation
  (`"Q"` or `"U"`)
- the sampling frequency, in Hz
- the CMB temperature (in K)
- the atmospheric temperature, in K
- the telescope temperature, in K

It returns 5 arrays, containing for each horn:
- the 3D vector containing the orientation of the horn in the sky
- the unique ID of the polarimeter associated to the horn
- the estimated bandwidth, in Hz
- the estimated noise temperature, in K
- the estimated knee frequency, in Hz
- the estimeted 1/f slope
        
Some corrections are applied to the the values taken directly from the DB:
- 1/f slopes greater than 2 are manually put equal to 2 since the Noise Generator used in the simulation works up to ths value.
- Knee frequencies are corrected by using the function 'fix_fknee'.
- For those polarimeters without valid measurements of the knee
frequency, the function randomly picks a value from a vector
containing all the available fknee measurements.

"""
function get_info_from_DB(db, horns, stokes, f_samp_hz, tcmb_k, tatm_k, ttel_k)

    @assert stokes in ["Q", "U"]
    
    orientations=[]
    polarimeters = Array{Int64}(undef, length(horns))
    β_hz_array = Array{Float64}(undef,length(horns))
    tnoise_k_array = Array{Float64}(undef,length(horns))
    fknee_hz_array = Array{Float64}(undef,length(horns))
    slope_array = Array{Float64}(undef,length(horns))

    measured_fknees = Float64[]
    measured_slopes = Float64[]

    for i in 1:length(horns)
        polarimeters[i] = db.focalplane[horns[i]].polid
        β_hz = db.detectors[polarimeters[i]].bandshape.bandwidth_hz
        tnoise_k = db.detectors[polarimeters[i]].tnoise.tnoise_k

        if(db.detectors[polarimeters[i]].spectrum.fknee_q_hz !=0)
            if stokes == "Q"
                fknee_hz = db.detectors[polarimeters[i]].spectrum.fknee_q_hz
                slope = db.detectors[polarimeters[i]].spectrum.slope_q
            else stokes == "U"
                fknee_hz = db.detectors[polarimeters[i]].spectrum.fknee_u_hz
                slope = db.detectors[polarimeters[i]].spectrum.slope_u
            end
            
            slope = clamp(slope, 0.0, 2.0)
            fknee_corrected = fix_fknee(f_samp_hz, tcmb_k, tatm_k, ttel_k, tnoise_k, β_hz, slope, fknee_hz)

            push!(measured_fknees, fknee_corrected)
            push!(measured_slopes, slope)

        end
    end
    rng = Random.MersenneTwister(1234)

    for i in 1:length(horns)
        append!(orientations, [db.focalplane[horns[i]].orientation])
        β_hz_array[i] = db.detectors[polarimeters[i]].bandshape.bandwidth_hz
        tnoise_k_array[i] = db.detectors[polarimeters[i]].tnoise.tnoise_k
        
        if(db.detectors[polarimeters[i]].spectrum.fknee_q_hz ==0)
            fknee_hz_array[i] = rand(rng, measured_fknees)
            slope_array[i] = rand(rng, measured_slopes)
        else
            if stokes == "Q"
                fknee = db.detectors[polarimeters[i]].spectrum.fknee_q_hz
                slope = db.detectors[polarimeters[i]].spectrum.slope_q
            else stokes == "U"
                fknee = db.detectors[polarimeters[i]].spectrum.fknee_u_hz
                slope = db.detectors[polarimeters[i]].spectrum.slope_u
            end    
            slope_array[i] = clamp(slope, 0.0, 2.0)
            fknee_hz_array[i]  = fix_fknee(f_samp_hz, tcmb_k, tatm_k, ttel_k, tnoise_k_array[i], β_hz_array[i], slope_array[i], fknee)

        end
    end
    return orientations, polarimeters, β_hz_array, tnoise_k_array, fknee_hz_array, slope_array
end
