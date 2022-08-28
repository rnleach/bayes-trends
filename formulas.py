""" Meteorlogical formulas for doing weather calculations.

This module should not rely on any modules that are not part of the 
standard python distribution.

References:

American Meteorological Society, 2021: Glossary of Meteorology. 
    Accessed May 26th, 2021. https://glossary.ametsoc.org/wiki/Welcome

Alduchov, O. A., and R. E. Eskridge, Improved Magnus' Form Approximation of 
    Saturation Vapor Pressure. National Climatic Data Center, 21 pp.

Djurić, D., 1994, Weather Analysis. Prentice Hall, 304 pp.

Rogers, R. R., and M. K. Yau, 1989, A Short Course in Cloud Physics, 3rd ed. 
    Butterworth Heinemann, 290 pp. ISBN 0-7506-3215-1.
"""
from math import exp, log

# Acceleration due to gravity at the Earth's surface. (m /s^2)
g = -9.81

# The gas constant for dry air. (J / K / kg)
R = 287.04

# The gas constant for water vapor. (J / K / kg)
Rv = 461.55

# Specific heat of dry air at constant pressure. (J / K / kg)
cp = 1004.0

# Specific heat of dry air at constant volume. (J / K / kg)
cv = 717.0

# Ratio of R / Rv. (unitless)
epsilon = R / Rv

# Ratio of cp and cv. (unitless)
gamma = cp / cv


def c_to_f(t_c):
    """Convert from Celsius to Fahrenheit."""
    if t_c is None:
        return None
    return 1.8 * t_c + 32.0


def f_to_c(t_f):
    """Convert from Fahrenheit to Celsius."""
    if t_f is None:
        return None
    return (t_f - 32.0) / 1.8


def c_to_k(t_c):
    """Convert Celsius to Kelvin."""
    if t_c is None:
        return None
    return t_c + 273.15


def k_to_c(t_k):
    """Convert Kelvin to Celsius."""
    if t_k is None:
        return None
    return t_k - 273.15


def inHg_to_hPa(p_inHg):
    """Convert inches of mercury to hectopascals."""
    if p_inHg is None:
        return None
    return p_inHg * 33.86389


def hPa_to_inHg(p_hPa):
    """Convert hectopascals to inches of mercury."""
    if p_hPa is None:
        return None
    return p_hPa / 33.86389


def knots_to_mps(knots):
    """Convert knots to meters per second."""
    if knots is None:
        return None
    return knots * 0.51444


def mps_to_knots(mps):
    """Convert meters per second to knots."""
    if mps is None:
        return None
    return mps / 0.51444


def m_to_km(meters):
    """Convert meters to kilometers."""
    if meters is None:
        return None
    
    return meters / 1000


def theta_kelvin(pressure_hpa, t_c):
    """Calculate the potential temperature in Kelvin.

    Potential temperature technically requires a reference 
    pressure to be defined. 1000 hPa is the most commonly 
    used and is assumed in this function.

    Returns:
    The potential temperature in K.
    """
    if pressure_hpa is None or t_c is None:
        return None
    return c_to_k(t_c) * (1000.0 / pressure_hpa)**(R / cp)


def temperature_c_from_theta(theta_kelvin, pressure_hpa):
    """Temperature from pressure and potential temperature.

    Returns:
    The temperature in °C.
    """
    if theta_kelvin is None or pressure_hpa is None:
        return None
    return k_to_c(theta_kelvin * (pressure_hpa / 1000.0)**(R / cp))


def vapor_pressure_liquid_water(t_c):
    """Calculate the vapor pressure over liquid water in hPa.

    From page 11 of Alduchov and Eskridge. The authors claim the formula
    is accurate from -80°C to 50°C.

    Returns:
    The vapor pressure over water in hPa.
    """
    
    if t_c is None:
        return None
    return 6.1037 * exp(17.641 * t_c / (t_c + 243.27))


def dew_point_from_vapor_pressure_over_liquid(vp_hpa):
    """Dew point from vapor pressure over liquid water.
    
    This function is the inverse of vapor_pressure_liquid_water.

    Returns:
    The dew point in °C.
    """
    if vp_hpa is None:
        return None
    a = log(vp_hpa / 6.1037) / 17.641
    return a * 243.27 / (1.0 - a)


def calc_rh(t_c, dp_c):
    """Calculate relative humidity.

    Returns:
    The relative humidity as a percent.
    """
    if t_c is None or dp_c is None:
        return None
    return 100.0 * vapor_pressure_liquid_water(dp_c) / vapor_pressure_liquid_water(t_c)


def calc_dp(t_c, rh):
    """Calculate the dew point in Celsius.

    Arguments:
    t_c - the temperature in °C.
    rh - the relative humidity as a percent, (0-100)

    Returns:
    The dew point in °C.
    """
    sat_vp = vapor_pressure_liquid_water(t_c)
    vp = sat_vp * rh / 100.0
    a = log(vp / 6.1037) / 17.641
    return a * 243.27 / (1.0 - a)


def mixing_ratio(dp_c, pressure_hpa):
    """Calculate the mixing ratio from the dew point and pressure.

    Returns:
    The mixing ratio in units of g/g or kg/kg.
    """
    if dp_c is None or pressure_hpa is None:
        return None
    vp = vapor_pressure_liquid_water(dp_c)
    return epsilon * vp / (pressure_hpa - vp)


def dew_point_from_p_and_mw(pressure_hpa, mw):
    """Dew point from pressure and mixing ratio.

    The mixing ratio is in units of g/g or kg/kg, NOT g/kg.

    Returns:
    The dew point in °C.
    """
    if pressure_hpa is None or mw is None:
        return None
    vp = mw * pressure_hpa / (mw + R / Rv)
    return dew_point_from_vapor_pressure_over_liquid(vp)


def virtual_temperature_c(t_c, dp_c, pressure_hpa):
    """Virtual temperature in Celsius.

    Source: AMS Glossary of Meteorology.

    Returns:
    The virtual temperature in °C.
    """
    if t_c is None or dp_c is None or pressure_hpa is None:
        return None
    rv = mixing_ratio(dp_c, pressure_hpa)
    t_k = theta_kelvin(pressure_hpa, t_c)
    vt_k = t_k * (1.0 + rv / epsilon) / (1.0 + rv)
    
    return temperature_c_from_theta(vt_k, pressure_hpa)


def temperature_kelvin_at_lcl(t_c, dp_c):
    """Approximate temperature at the Lifting Condensation Level (LCL).

       Eq 5.17 from "Weather Analysis" by Dušan Dujrić

       Returns:
       temperature K at the lifting condensation level
    """
    if t_c is None or dp_c is None:
        return None
    
    if dp_c >= t_c:
        return c_to_k(t_c)
    
    celsius_lcl = dp_c - (0.001296 * dp_c + 0.1963) * (t_c - dp_c)
    return c_to_k(celsius_lcl)


def press_and_temp_k_at_lcl(t_c, dp_c, pressure_hpa):
    """Temperature and pressure at the lifting condensation level (LCL).

       Eqs 5.17 and 5.18 from "Weather Analysis" by Dušan Dujrić

       Returns:
       tuple (pressure hPa, temperature K)
    """
    if t_c is None or dp_c is None or pressure_hpa is None:
        return (None, None)
    
    if dp_c >= t_c:
        # It's either saturated or super saturated, so we're already at the LCL.
        return (pressure_hpa, c_to_k(t_c))
    
    t_lcl = temperature_kelvin_at_lcl(t_c, dp_c)
    t_kelvin = c_to_k(t_c)
    p_lcl = pressure_hpa * (t_lcl / t_kelvin)**(cp / R)
    return (p_lcl, t_lcl)


def latent_heat_of_condensation(t_c):
    """Latent heat of condensation for water.

    Polynomial curve fit to Table 2.1 on p. 16 in Rogers and Yau.

    Returns:
    The latent of condensation in J/kg.
    """
    if t_c is None:
        return None
    # The table has values from -40.0 to 40.0. So from -100.0 to -40.0 is
    # actually an exrapolation. I graphed the values from the extrapolation,
    # and the curve looks good, and is approaching the latent heat of
    # sublimation, but does not exceed it. This seems very reasonable to me,
    # especially considering that a common approximation is to just use a
    # constant value.
    if t_c < -100.0 or t_c > 60.0:
        return None
    
    return (2500.8 - 2.36 * t_c + 0.0016 * t_c * t_c - 0.00006 * t_c * t_c * t_c) * 1000.0


def specific_humidity(dp_c, pressure_hpa):
    """Calculate the specific humidity.

    Eqs 5.11 and 5.12 from "Weather Analysis" by Dušan Dujrić

    Returns:
    The specific humidity as a decimal, or in g/g or kg/kg and NOT g/kg.
    """
    if dp_c is None or pressure_hpa is None:
        return None
    vp = vapor_pressure_liquid_water(dp_c)
    assert vp > 0
    assert pressure_hpa > 0
    
    return vp / pressure_hpa * epsilon


def dew_point_from_p_and_specific_humidity(pressure_hpa, specific_humidity):
    """Dew point from specific humidity and pressure.
    
    This is the inverse of the specific_humidity equation above.

    Returns:
    The dew point in °C.
    """
    assert pressure_hpa > 0
    assert specific_humidity > 0

    if pressure_hpa is None or specific_humidity is None:
        return None
    vp = specific_humidity * pressure_hpa / epsilon
    return dew_point_from_vapor_pressure_over_liquid(vp)


def theta_e_kelvin(t_c, dp_c, pressure_hpa):
    """Calculate equivalent potential temperature.

    Eq 5.23 from "Weather Analysis" by Dušan Dujrić

    Returns:
    The equivalent potential temperature in K.
    """
    if t_c is None or dp_c is None or pressure_hpa is None:
        return None
    
    lc = latent_heat_of_condensation(t_c)
    if lc is None:
        return None
    
    theta = theta_kelvin(pressure_hpa, t_c)
    t_lcl = temperature_kelvin_at_lcl(t_c, dp_c)
    qs = specific_humidity(dp_c, pressure_hpa)
    return theta * (1.0 + lc * qs / (cp * t_lcl))


def theta_e_saturated_kelvin(pressure_hpa, t_c):
    """Equivalent potential temperature assuming saturation.

    Eq 5.23 from "Weather Analysis" by Dušan Dujrić

    Returns:
    The equivalent potential temperature in K.
    """
    if pressure_hpa is None or t_c is None:
        return None
    
    lc = latent_heat_of_condensation(t_c)
    if lc is None:
        return None
    
    theta = theta_kelvin(pressure_hpa, t_c)
    qs = specific_humidity(t_c, pressure_hpa)
    
    return theta * (1.0 + lc * qs / (cp * c_to_k(t_c)))


def temperature_c_from_theta_e_saturated_and_pressure(pressure_hpa, theta_e_k):
    """Temperature from equivalent potential temperature and pressure.
    
    Assume saturation. This is the inverse of theta_e_saturated_kelvin.
    """
    if pressure_hpa is None or theta_e_k is None:
        return None
    
    def func_to_minimize(t_c):
        te = theta_e_saturated_kelvin(pressure_hpa, t_c)
        if te is None:
            return None
        return te - theta_e_k
    
    try:
        return find_root(func_to_minimize, -80.0, 50.0)
    except Exception:
        return None


# FIXME: Update this to use Brent's method, which usually converges faster.
def find_root(func, low_val, high_val):
    """Bisection algorithm for finding the root of an equation.
    
    low_val is the left bracket, func(low_val) < 0
    high_val is the right bracket, func(high_val) > 0
    func is a function that takes a single argument and returns a 
        single floating point value or None.

    Assumes that there is only one root between low_val and high_val.

    Used when finding wet bulb temperature.
    """
    if func is None or low_val is None or high_val is None:
        return None
    
    MAX_IT = 50
    EPS = 1.0e-10
    
    if low_val > high_val:
        temp = high_val
        high_val = low_val
        low_val = temp
    
    f_low = func(low_val)
    f_high = func(high_val)
    if f_low is None or f_high is None:
        return None
    
    # Check to make sure we have bracketed a root.
    if f_high * f_low > 0.0:
        raise Exception(
            "Failed to bracket minimum! f_high(%f) = %f and f_low(%f) = %f" %
            (high_val, f_high, low_val, f_low)
        )
    
    mid_val = (high_val - low_val) / 2.0 + low_val
    f_mid = func(mid_val)
    if f_mid is None:
        return None
    
    for _ in range(MAX_IT):
        if f_mid * f_low > 0.0:
            low_val = mid_val
            f_low = f_mid
        else:
            high_val = mid_val
        
        if abs(high_val - low_val) < EPS:
            break
        
        mid_val = (high_val - low_val) / 2.0 + low_val
        f_mid = func(mid_val)
        if f_mid is None:
            return None
    
    return mid_val


if __name__ == "__main__":
    from math import isclose
    import numpy as np
    
    #
    # Test temperature conversions
    #
    for t in range(-40, 40, 1):
        assert isclose(t, f_to_c(c_to_f(t)), abs_tol=1e-8)
        assert isclose(t, k_to_c(c_to_k(t)), abs_tol=1e-8)
    
    #
    # Test potential temperature
    #
    for t in np.arange(-40, 40, 1):
        for p in np.arange(1060.0, 100.0, -10.0):
            assert isclose(t, temperature_c_from_theta(theta_kelvin(p, t), p), abs_tol=1e-8)
    
    #
    # Test vapor pressure
    #
    for t in np.arange(-40, 40, 1):
        assert isclose(
            t,
            dew_point_from_vapor_pressure_over_liquid(vapor_pressure_liquid_water(t)),
            abs_tol=1e-8
        )
    
    #
    # Test RH
    #
    for t in range(-39, 40, 1):
        for dp in range(-40, t, 1):
            assert isclose(dp, calc_dp(t, calc_rh(t, dp)), abs_tol=1e-8)
    
    #
    # Test mixing ratio
    #
    for dp in np.arange(-40, 40, 1):
        for p in np.arange(1060.0, 100.0, -10.0):
            assert isclose(dp, dew_point_from_p_and_mw(p, mixing_ratio(dp, p)), abs_tol=1e-8)
    
    #
    # Test Virtual Temperature
    #
    for t in np.arange(-39.0, 40.0, 1):
        for dp in np.arange(-40, t, 1):
            for p in np.arange(1060.0, 100.0, -10.0):
                vt = virtual_temperature_c(t, dp, p)
                assert vt >= t
    
    #
    # Test theta_e_saturated
    #
    for t in np.arange(-40, 40, 1):
        for p in np.arange(1060.0, 100.0, -10.0):
            theta_e = theta_e_saturated_kelvin(p, t)
            t_back = temperature_c_from_theta_e_saturated_and_pressure(p, theta_e)
            assert isclose(t, t_back, abs_tol=1e-8)
