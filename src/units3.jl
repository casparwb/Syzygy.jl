module SyzygyUnits
export unit_length, unit_mass, unit_time
using DynamicQuantities#, @us_str, @u_str, @register_unit
# using Syzygy.DynamicQuantities
# using Reexport
# @reexport using DynamicQuantities
# unit values following https://iopscience.iop.org/article/10.3847/0004-6256/152/2/41

# try
# catch err
#     @show err
# end
# function __init__()
# end


# # function set_units(unit_length=, 
#                    unit_time=, 
#                    unit_mass=)

#     @assert isone(ulength(uexpand(unit_length)))
#     @assert isone(umass(uexpand(unit_mass)))
#     @assert isone(utime(uexpand(unit_time))) 

#     return unit_length, unit_mass, unit_time
# end

if !(:Rsun in DynamicQuantities.UNIT_SYMBOLS) 
    @register_unit Rsun 6.975e8u"m" 
end
if !(:Msun in DynamicQuantities.UNIT_SYMBOLS) 
    @register_unit Msun 1.9884754153381438e30u"kg" 
end
if !(:Lsun in DynamicQuantities.UNIT_SYMBOLS) 
    @register_unit Lsun 3.828e26u"W" 
end
if !(:AU in DynamicQuantities.UNIT_SYMBOLS)
    @register_unit AU 149597870700u"m" 
end
# if !(:au in DynamicQuantities.UNIT_SYMBOLS) 
#     @register_unit au 149597870700u"m" 
# end



unit_length, unit_mass, unit_time = u"m", u"kg", u"s"#us"Rsun", us"Msun", us"yr"

end

    # @register_unit Rsun 6.975e8u"m"
    # @register_unit Msun 1.9884754153381438e30u"kg"
    # @register_unit Lsun 3.828e26u"W"