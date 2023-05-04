import numpy as np

def get_cylindrical_th(internal_pressure, inner_diameter, design_stress, joint_factor=1, alt=False, cr=1):
    Pi = internal_pressure * pow(10, 5)
    Di = inner_diameter / 1000
    if alt:
        return (((Pi*(Di-cr*pow(10, -3)))/(2*design_stress*joint_factor)) + cr*pow(10, -3)) * 1000
    return ((Pi * Di) / ((2 * design_stress * joint_factor) - Pi)) * 1000

def get_flat_end_th(internal_pressure, nominal_diameter, design_stress, full_face=True):
    Pi = internal_pressure * pow(10, 5)
    De = nominal_diameter / 1000
    if full_face:
        Cp = 0.4
    else:
        Cp = 0.55
    return (Cp * De * np.sqrt(Pi/design_stress)) * 1000

def get_hemispherical_th(cylindrical_th):
    return cylindrical_th * 0.6  # C&R page 819 for more info

def get_ellipsoid_th(internal_pressure, inner_diameter, design_stress, joint_factor):
    # For a major and minor axis ratio of 2:1 (most common)
    Pi = internal_pressure * pow(10, 5)
    Di = inner_diameter / 1000
    return ((Pi * Di) / ((2 * design_stress * joint_factor) - (0.2 * Pi))) * 1000

def get_ts_thickness(shell_inner_diameter, wall_th, design_pressure, design_stress, bundle_diameter, min_thickness):
    F = 1
    G = pow(((shell_inner_diameter/2) - 0.8 - wall_th)/1000, 2) * np.pi
    P = design_pressure * pow(10, 5)
    S = design_stress
    K = 0.785
    n = 1 - (K/pow(1.25, 2))
    bending = (1/3) * F * G * np.sqrt(P / (n * S))
    C = (bundle_diameter / 1000) * np.pi
    A = np.pi * pow((bundle_diameter / 2000), 2)
    Dl = 4 * A / C
    shear = (0.31 * Dl * P) / (S * (1 - (1/1.25)))
    return max([bending, shear, min_thickness])


ss304_tensile_strength = 510 * pow(10, 6)
ss304_design_stress = 145 * pow(10, 6)
a20_tensile_strength = 560 * pow(10, 6)
a20_design_stress = 150 * pow(10, 6)
# Ellipsoid is the most common for my range
