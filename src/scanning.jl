# -*- encoding: utf-8 -*-
#
# MIT License
#
# Copyright (c) 2016 Maurizio Tomasi
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#
# # TERMINOLOGY
#
# The code in this file tries to follow the conventions laid in
# document LSPE-STRIP-SP-017 ("STRIP Telescope Coordinate Systems
# Specification") as much as possible. As the reader might not be
# acquainted with it, here is a short summary of the relevant
# terminology:
#
# +-----+-----------------------------+------------------------------------------+
# | SCS | Site Coordinate System      | Centered on the GPS station              |
# +-----+-----------------------------+------------------------------------------+
# | MCS | Mount Coordinate System     | Centered at the basis of the telescope   |
# +-------+---------------------------+------------------------------------------+
# | ACS | A zimuth Coordinate System  | Moveable, centered in the middle of      |
# |     |                             | ground motor. Positive angles go from    |
# |     |                             | North to *East* (beware!)                |
# +-----+-----------------------------+------------------------------------------+
# | ECS | Elevation Coordinate System | Moveable, the Y axis is aligned with     |
# |     |                             | the motor axis                           |
# +-----+-----------------------------+------------------------------------------+
# | BCS | Boresight Coordinate System | It should have been moveable, but the    |
# |     |                             | current design does not allow it (alas!) |
# +-----+-----------------------------+------------------------------------------+
# | RDP | Reference Detector Plane    | Centered in the focus of the telescope.  |
# |     |                             | All the beam vectors are specified in    |
# |     |                             | this frame.                              |
# +-----+-----------------------------+------------------------------------------+
#
# We always use a right-handed coordinate system:
#
#                 ^ Z axis (aligned with "up" direction, e.g., Zenith)
#                 |
#                 |
#                 |
#                 |
#                 |
#                 |
#                 |
#                 |
#                 |
#          Origin O----------------> Y axis
#                /
#               /
#              /
#             /
#            /
#           v  X axis (aligned with North, when applicable)
#

import Healpix
import StaticArrays
import LinearAlgebra: ×, dot
import AstroLib
import Dates

export TENERIFE_LATITUDE_DEG, TENERIFE_LONGITUDE_DEG, TENERIFE_HEIGHT_M
export timetorotang, telescopetoground, groundtoearth
export genpointings, polarizationangle

const TENERIFE_LATITUDE_DEG = 28.3
const TENERIFE_LONGITUDE_DEG = -16.509722
const TENERIFE_HEIGHT_M = 2390

include("quaternions.jl")

"""
    timetorotang(time, rpm)

Convert a time into a rotation angle, given the number of rotations per minute.
The time should be expressed in seconds. The return value is in radians.
`time` can either be a scalar or a vector.
"""
function timetorotang(time_s, rpm)
    if rpm == 0
        0.0
    else
        2 * π * time_s * (rpm / 60)
    end
end


@doc raw"""
    telescopetoground(wheelanglesfn, time_s)

Return a quaternion of type `Quaternion{Float64}` representing the
coordinate transform from the focal plane to the ground of the
telescope. The parameter `wheelanglesfn` must be a function which
takes as input a time, `time_s`, in seconds, and it must return a
3-tuple containing the angles of the following motors, measured in
**radians**:

1. The boresight motor (rotation around the ``z`` axis,
   counterclockwise)

2. The altitude motor (rotation around the ``y`` axis,
   counterclockwise)

3. The ground motor (rotation around the ``z`` axis, **clockwise**:
   N→E→S→W)

# Example

`````julia
telescopetoground(3600.0) do
    # Boresight motor keeps a constant angle equal to 0°
    # Altitude motor remains at 20° from the Zenith
    # Ground motor spins at 1 RPM
    (0.0, deg2rad(20.0), timetorotang(time_s, 1))
end
`````
"""
function telescopetoground(wheelanglesfn, time_s)
    (wheel1ang, wheel2ang, wheel3ang) = wheelanglesfn(time_s)

    qwheel1 = qrotation_z(wheel1ang)
    qwheel2 = qrotation_y(wheel2ang)

    # The minus sign here takes into account the fact that the azimuth
    # motor requires positive angles to turn North into East
    qwheel3 = qrotation_z(-wheel3ang)

    qwheel3 * (qwheel2 * qwheel1)
end


"""
    groundtoearth(groundq, time_s, latitude_deg; day_duration_s=86400.0)

Return a quaternion of type `Quaternion{Float64}` representing the
coordinate transformation from the ground of the telescope to the
Equatorial coordinate system. The parameter `groundq` must be a
quaternion describing the coordinate transformation from the focal
plane of the telescope to the ground. The parameter `time_s` must be a
time in seconds, and `latitude_deg` is the latitude (in degrees, N is
positive) of the location where the observation is made.

The keyword `day_duration_s` specifies the length of a sidereal day in
seconds.
"""
function groundtoearth(groundq, time_s, latitude_deg; day_duration_s = 86400.0)
    locq = qrotation_y(deg2rad(90 - latitude_deg))
    earthq = qrotation_z(2π * time_s / day_duration_s)

    earthq * (locq * groundq)
end


"""
    vector2equatorial(dir, jd, latitude_deg, longitude_deg, height_m; 
                      prec=true, nut=true, aber=true)
<<<<<<< HEAD
=======

>>>>>>> e4dd169f9952e26b16292dccd44ce64cf69db67e
Transform the Healpix coordinates of a vector into equatorial coordinates.
The parameter `vector` is 3D vector, `jd` is the julian date. The
paramters `latitude_deg`, `longitude_deg` and `height_m` should contain the
latitude (in degrees, N is positive), the longitude (in degrees, counterclockwise
is positive) and the height (in meters) of the location where the observation is
made. The parameters `prec`, `nut` and `aber` allow to consider secondary effects
respectively of precession, nutation and aberration. 
"""
function vector2equatorial(vector, jd, latitude_deg, longitude_deg, height_m,
                           prec, nut, aber, ref)
    (θ, ϕ) = Healpix.vec2ang(vector...)
    alt_rad = π / 2 - θ
    az_rad = 2π - ϕ

    ra_deg, dec_deg, _ = AstroLib.hor2eq(rad2deg(alt_rad),
                                         rad2deg(az_rad),
                                         jd,
                                         latitude_deg,
                                         longitude_deg,
                                         height_m,
                                         precession = prec,
                                         nutate = nut,
                                         aberration = aber,
                                         refract = ref)
    (deg2rad(dec_deg), deg2rad(ra_deg))
end


"""
    polarizationangle(northdir, poldir)

Calculate the polarization angle projected in the sky in IAU
conventions.  The parameter `northdir` must be a versor that points
the North and `poldir` must be a versor that identify the polarization
direction projected in the sky.  The return value is in radians.
"""
function polarizationangle(northdir, poldir)
    cosψ = clamp(dot(northdir, poldir), -1, 1)
    crosspr = northdir × poldir
    sinψ = clamp(sqrt(dot(crosspr, crosspr)), -1, 1)
    ψ = atan(sinψ, cosψ)
    ψ
end


"""
    get_northdir(θ, ϕ)

Compute the North of a vector. The North for a vector v is -dv/dθ, as
θ is the colatitude and moves along the meridian.
"""
get_northdir(θ, ϕ) = StaticArrays.@SArray [-cos(θ) * cos(ϕ), -cos(θ) * sin(ϕ), sin(θ)]


"""
    quat_to_angles(boreaxis, polaxis, quat)

Transform the boresight direction `dir` and the polarization direction
`poldir` according to quaternion `quat`. Return the 3-tuple (θ, ϕ, ψ)
representing the colatitude, longitude, and polarization angle
(calculated northward).

This function is used internally by [`genpointings`](@ref).
"""
function quat_to_angles(boreaxis, polaxis, quat)
    rotmatr = rotationmatrix_normalized(quat)
    boresight = rotmatr * boreaxis
    poldir = rotmatr * polaxis

    (θ, ϕ) = Healpix.vec2ang(boresight...)

    northdir = get_northdir(θ, ϕ)
    (θ, ϕ, polarizationangle(northdir, poldir))
end


function genpointings(wheelanglesfn,
                      dir,
                      timerange_s;
                      polaxis = Float64[1.0, 0.0, 0.0],
                      latitude_deg = 0.0,
                      ground = false,
                      day_duration_s = 86400.0)

    if ground
        dirs = Array{Float64}(undef, length(timerange_s), 4)
        ψ = Array{Float64}(undef, length(timerange_s), 2)
    else
        dirs = Array{Float64}(undef, length(timerange_s), 2)
        ψ = Array{Float64}(undef, length(timerange_s))
    end

    for (idx, time_s) = enumerate(timerange_s)

        # This converts the RDP into the MCS (ground reference frame)
        groundq = telescopetoground(wheelanglesfn, time_s)
        # This converts the MCS into the celestial reference frame
        quat = groundtoearth(groundq, time_s, latitude_deg; day_duration_s = day_duration_s)

        θ, ϕ, curψ = quat_to_angles(dir, polaxis, quat)
        (dirs[idx, 1], dirs[idx, 2]) = (θ, ϕ)

        if ground
            # Re-run the transformation algorithm using the ground quaternion
            θ_ground, ϕ_ground, ψ_ground = quat_to_angles(dir, polaxis, groundq)

            (dirs[idx, 3], dirs[idx, 4]) = (θ_ground, ϕ_ground)
            (ψ[idx, 1], ψ[idx, 2]) = (curψ, ψ_ground)
        else
            ψ[idx] = curψ
        end
    end

    (dirs, ψ)
end


function genpointings(wheelanglesfn,
                      dir,
                      timerange_s,
                      t_start;
                      polaxis = Float64[1.0, 0.0, 0.0],
                      latitude_deg = 0.0,
                      longitude_deg = 0.0,
                      height_m = 0,
                      precession = true,
                      nutation = true,
                      aberration = true,
                      refract = true)

    skydirs = Array{Float64}(undef, length(timerange_s), 2)
    skyψ = Array{Float64}(undef, length(timerange_s))

    for (idx, time_s) = enumerate(timerange_s)
        groundq = telescopetoground(wheelanglesfn, time_s)
        rotmatr = rotationmatrix_normalized(groundq)
        vector = rotmatr * dir

        jd = AstroLib.jdcnv(t_start + Dates.Nanosecond(round(Int64, time_s*1e9)))
        Dec_rad, Ra_rad = vector2equatorial(vector,
                                            jd,
                                            latitude_deg,
                                            longitude_deg,
                                            height_m,
                                            precession,
                                            nutation,
                                            aberration,
                                            refract)

        poldir = rotmatr * polaxis
        northdir = get_northdir(π / 2 - Dec_rad, Ra_rad)
        
        skydirs[idx, 1] = Dec_rad
        skydirs[idx, 2] = Ra_rad
        skyψ[idx] = polarizationangle(northdir, poldir)
    end

    (skydirs, skyψ)
end


@doc raw"""
    genpointings(wheelanglesfn, dir, timerange_s; 
                 polaxis=Float64[1.0, 0.0, 0.0],
                 latitude_deg=0.0, ground=false)
    genpointings(wheelanglesfn, dir, timerange_s, t_start;
                 polaxis=Float64[1.0, 0.0, 0.0],
                 latitude_deg=0.0, longitude_deg=0.0, height_m=0.0)

Generate a set of pointings for a STRIP detector. The parameter
`wheelanglesfn` must be a function which takes as input a time in
seconds and returns a 3-tuple containing the angles (in radians) of
the three motors:

1. The boresight motor (rotation around the ``z`` axis,
   counterclockwise)

2. The altitude motor (rotation around the ``y`` axis,
   counterclockwise)

3. The ground motor (rotation around the ``z`` axis, **clockwise**:
   N→E→S→W)

The meaning of the parameters/keywords is the following:

- `dir` specifies the pointing direction of the mean (the boresight is
  [0, 0, 1]). It must be normalized.

- `timerange_s` is an enumerable type that specifies at which times
  (in seconds) pointings must be computed.

- `t_start` is a DateTime which tells the exact UTC date and time of the 
  observation.

- `latitude_deg` is the latitude of the location where the observation
  is made (in degrees, North is positive).

- `ground` is a Boolean: if `true`, the function will return a 4-tuple
  containing the colatitude and longitude measured in Equatorial
  coordinates (columns 1 and 2) and in ground coordinates (columns 3
  and 4). If false, only the Equatorial coordinates are computed.

- `polaxis` is the polarization axis; it must be normalized.

- `precession` is a boolen parameter that allow to include the Earth's precession
  effect; default is true.
<<<<<<< HEAD
- `nutation` is a boolen parameter that allow to include the Earth's nutation
  effect; default is true.
- `aberration` is a boolen parameter that allow to include the stellar aberration
  effect; default is true.
- `refract` is a boolen parameter that allow to include the refraction correction
  ; default is true.

=======

- `nutation` is a boolen parameter that allow to include the Earth's nutation
  effect; default is true.

- `aberration` is a boolen parameter that allow to include the stellar aberration
  effect; default is true.

- `refract` is a boolen parameter that allow to include the refraction correction
  ; default is true.
>>>>>>> e4dd169f9952e26b16292dccd44ce64cf69db67e

# Return values

In the first form, the function returns a 2-tuple containing the sky
directions (a N×2 array containing declination and right ascension, in
Equatorial coordinates) and the polarization angle for each time step.

In the second form, the function returns a 2-tuple (4-tuple)
containing the directions (a N×2 or Nx4 array containing the
colatitude and the longitude) and the polarization angles at each time
step.


# Examples

Here is an example of the first kind:

`````julia
dir, psi = genpointings([0, 0, 1], 0:0.1:1) do time_s
    # Boresight motor keeps a constant angle equal to 0°
    # Altitude motor remains at 20° from the Zenith
    # Ground motor spins at 1 RPM
    return (0.0, deg2rad(20.0), timetorotang(time_s, 1))
end
`````

And here is an example of the second kind. Note that, unlike the
previous example, we use a lambda function instead of the `do...end`
notation.

`````julia
import Dates

dirs, psi = genpointings(time_s -> (0, deg2rad(20),
                                    timetorotang(time_s, 1)),
                         [0, 0, 1],
                         0:0.1:1,
                         Dates.DateTime(2019, 01, 01, 0, 0, 0),
                         latitude_deg=10.0,
                         longitude_deg=20.0,
                         height_m = 1000) do time_s
`````
"""
genpointings
