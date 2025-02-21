
let G = 6.6743015e-11u"m^3/kg/s^2", c = 299_792_458u"m/s"

    mass_factor = ustrip(u"kg/m", c^2/G)*u"kg"
    time_factor = ustrip(u"s/m", 1/c)*u"s"

    @unit geometer "geometer" GeometricMeter 1u"m" false
    @unit geomass "geomass" GeometricMass mass_factor false
    @unit geotime "geotime" GeometricTime time_factor false
end

let G = 6.6743015e-11u"m^3/kg/s^2", c = 299_792_458u"m/s"

    mass_factor = u"kg"(1.0u"Msun")
    length_factor = u"m"(1.0u"Rsun")
    time_factor = ustrip(u"m"^(3/2)/u"kg"^(1/2), sqrt(length_factor^3/mass_factor))*u"s"

    @unit nbmeter "nbmeter" NBodyMeter length_factor false
    @unit nbmass "nbmass" NBodyMass mass_factor false
    @unit nbtime "nbtime" NBodyTime time_factor false
end

function set_units(units)
    allowed_units = ["SI", "Solar", "CGS", "Geometric", "Nbody", "Solar_s"]

    if !(units in allowed_units)
        throw(ArgumentError("Given unit system not allowed. Currently $(allowed_units) are allowed.")) 
    end

    @set_preferences!("units" => units)
    @info "New default units set; restart your Julia session for this change to take effect!"
end

function units_from_unit_system(unit_system)
    allowed_units = ["SI", "Solar", "SI", "Geometric", "Nbody", "Solar_s"]
    
    units = if unit_system == "SI"
        u"kg, m, s"
    elseif unit_system == "CGS"
        u"g, cm, s"
    elseif unit_system == "Solar"
        u"Rsun, Msun, yr"
        elseif unit_system == "Solar_s"
        u"Rsun, Msun, s"
    elseif unit_system == "Geometric"
        geometer, geomass, geotime
    elseif unit_system == "Nbody"
        nbmeter, nbmass, nbtime
    else
        throw(ArgumentError("Given unit system not allowed. Currently $(allowed_units) are allowed.")) 
    end

    return units
end

const unit_system = @load_preference("units", "Solar")
const units = units_from_unit_system(unit_system)
const unit_length = only([u for u in units if Unitful.dimension(u) == Unitful.𝐋])
const unit_mass = only([u for u in units if Unitful.dimension(u) == Unitful.𝐌])
const unit_time = only([u for u in units if Unitful.dimension(u) == Unitful.𝐓])

Unitful.preferunits(units...)
const localpromotion = copy(Unitful.promotion)