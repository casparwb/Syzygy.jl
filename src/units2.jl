# module SyzygyUnits

abstract type SyzygyUnitSystem end

# struct SIUnits{tL, tT, tM} <: SyzygyUnitSystem
#     unit_length::tL
#     unit_time::tT
#     unit_mass::tM

#     function SIUnits()
#         units = (u"m", u"s", u"kg")
#         ts = typeof.(units)
#         return new{ts...}(units...)
# end

# struct CGSUnits{tL, tT, tM} <: SyzygyUnitSystem
#     unit_length::tL
#     unit_time::tT
#     unit_mass::tM

#     function CGSUnits()
#         units = (u"cm", u"s", u"g")
#         ts = typeof.(units)
#         return new{ts...}(units...)
# end

struct SIUnits{T} <: SyzygyUnitSystem
    name::String
    units::T

    function SIUnits()
        units = (u"m", u"s", u"kg")
        return new{typeof(units)}("SI", units)
    end
end

struct CGSUnits{T} <: SyzygyUnitSystem
    name::String
    units::T

    function CGSUnits()
        units = (u"cm", u"s", u"g")
        return new{typeof(units)}("CGS", units)
    end
end

function set_units!(units::SyzygyUnitSystem)
    @set_preferences!("units" => (units.name, units.units))
    @info "New default units set; restart your Julia session for this change to take effect!"

    return nothing
end
    


# let G = 6.6743015e-11u"m^3/kg/s^2", c = 299_792_458u"m/s"

#     mass_factor = u"kg"(1.0u"Msun")
#     length_factor = u"m"(1.0u"Rsun")
#     time_factor = ustrip(u"m"^(3/2)/u"kg"^(1/2), sqrt(length_factor^3/mass_factor))*u"s"

#     @unit nbmeter "nbmeter" NBodyMeter length_factor false
#     @unit nbmass "nbmass" NBodyMass mass_factor false
#     @unit nbtime "nbtime" NBodyTime time_factor false
# end

# function set_units(units)
#     allowed_units = ["SI", "Solar", "CGS", "Geometric", "Nbody", "Solar_s"]

#     if !(units in allowed_units)
#         throw(ArgumentError("Given unit system not allowed. Currently $(allowed_units) are allowed.")) 
#     end

#     @set_preferences!("units" => units)
#     @info "New default units set; restart your Julia session for this change to take effect!"
# end

# function units_from_unit_system(unit_system)
#     allowed_units = ["SI", "Solar", "SI", "Geometric", "Nbody", "Solar_s"]
    
#     units = if unit_system == "SI"
#         u"kg, m, s"
#     elseif unit_system == "CGS"
#         u"g, cm, s"
#     elseif unit_system == "Solar"
#         u"Rsun, Msun, yr"
#         elseif unit_system == "Solar_s"
#         u"Rsun, Msun, s"
#     elseif unit_system == "Geometric"
#         geometer, geomass, geotime
#     elseif unit_system == "Nbody"
#         nbmeter, nbmass, nbtime
#     else
#         throw(ArgumentError("Given unit system not allowed. Currently $(allowed_units) are allowed.")) 
#     end

#     return units
# end

const unit_system, units = @load_preference("units", ("SI", ("m", "s", "kg")))
# const units = units_from_unit_system(unit_system)
const units = uparse.(units)
const unit_length = only([u for u in units if Unitful.dimension(u) == Unitful.𝐋])
const unit_mass = only([u for u in units if Unitful.dimension(u) == Unitful.𝐌])
const unit_time = only([u for u in units if Unitful.dimension(u) == Unitful.𝐓])

Unitful.preferunits(units...)
const localpromotion = copy(Unitful.promotion)

# end