import numpy as np
import image_to_poly as itp


# Conversion and integration functions
def f_to_c(temps):
    if type(temps) == list:
        temps = np.asarray(temps)
    return ((temps - 32) * 5) / 9


def thermal_i_to_m(values):
    if type(values) == list:
        values = np.asarray(values)
    return values * (1055.055853 / (60 * 60 * 0.3048 * (5 / 9)))


def yield_i_to_m(values):
    if type(values) == list:
        values = np.asarray(values)
    return 6894757.2931783 * values


def integrate(coefficients, uv, lv):
    output = 0
    order = len(coefficients)
    for _i in range(order):
        index = _i + 1
        output += ((pow(uv, index) / index) - (pow(lv, index) / index)) * coefficients[order - index]
    return output


def evaluate_polynomial(coefficients, value):
    order = len(coefficients) - 1
    output = 0
    for index in range(order + 1):
        output += coefficients[order - index] * pow(value, index)

    return output


tube_gap_count = 5
pipe_schedules = ["5S", "10S", "40S", "80S", "10", "20", "30", "40", "STD", "60", "XS", "80", "a", "b", "c"]
with open("ASME_pipe_sizes_1.txt", "r") as _file:
    pipe_sizes = {float(line.split(" ")[2]): [line.split(" ")[0],
                                              {pipe_schedules[index]: (float(val) if val != "-" else None)
                                               for index, val in enumerate(line.split(" ")[4::3])}]
                  for line in _file.readlines()}

# Data for 2-propanol
inlet_temp_ipa = 90
outlet_temp_ipa = 40
bulk_temp_ipa = (inlet_temp_ipa + outlet_temp_ipa) / 2

t_ipa = [0, 20, 25, 40, 60, 80, 100, 120, 140, 160, 180, 200]
cp_ipa_cal = np.asarray([0.541, 0.596, 0.613, 0.67, 0.742, 0.798, 0.848, 0.894, 0.939, 0.983, 1.024, 1.062])
cp_ipa = 1000 * 4.1833 * cp_ipa_cal
c_ipa = np.polyfit(t_ipa, cp_ipa, 6)
mean_cp_ipa = integrate(c_ipa, inlet_temp_ipa, outlet_temp_ipa) / (inlet_temp_ipa - outlet_temp_ipa)
mean_cp_ipa_bulk = evaluate_polynomial(c_ipa, bulk_temp_ipa)

rho_t_ipa = list(range(0, 110, 10))
rho_ipa_vl = np.asarray([0.8003, 0.7931, 0.7885, 0.7774, 0.7689, 0.7599, 0.7504, 0.7402, 0.7295, 0.7182, 0.7062]) * 1000
rho_c_ipa = np.polyfit(rho_t_ipa, rho_ipa_vl, 6)
rho_ipa = integrate(rho_c_ipa, inlet_temp_ipa, outlet_temp_ipa) / (inlet_temp_ipa - outlet_temp_ipa)
rho_ipa_bulk = evaluate_polynomial(rho_c_ipa, bulk_temp_ipa)

mfr_ipa = 85000 / (60 * 60)
vfr_ipa = mfr_ipa / rho_ipa

mu_t_ipa = [0, 20, 50, 80, 100]
mu_values_ipa = np.asarray([4.710, 2.432, 1.042, 0.524, 0.354]) / 1000
mu_c_ipa = np.polyfit(mu_t_ipa, mu_values_ipa, 4)
mu_ipa = integrate(mu_c_ipa, inlet_temp_ipa, outlet_temp_ipa) / (inlet_temp_ipa - outlet_temp_ipa)
mu_ipa_bulk = evaluate_polynomial(mu_c_ipa, bulk_temp_ipa)

k_t_ipa = [-25, 0, 25, 50, 75, 100]
k_values_ipa = [0.146, 0.141, 0.135, 0.129, 0.124, 0.118]
k_c_ipa = np.polyfit(k_t_ipa, k_values_ipa, 4)
k_ipa = integrate(k_c_ipa, inlet_temp_ipa, outlet_temp_ipa) / (inlet_temp_ipa - outlet_temp_ipa)
k_ipa_bulk = evaluate_polynomial(k_c_ipa, bulk_temp_ipa)

# Data for water
inlet_temp_wa = 25

t_wa = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
cp_wa = [4217.6, 4192.1, 4181.8, 4178.4, 4178.5, 4180.6, 4184.3, 4189.5, 4196.3, 4205, 4215.9]
c_wa = np.polyfit(t_wa, cp_wa, 6)

rho_values_wa = [999.84, 999.7, 998.21, 995.65, 992.22, 988.03, 983.2, 977.78, 971.82, 965.35, 958.4]
rho_c_wa = np.polyfit(t_wa, rho_values_wa, 6)

mu_values_wa = [0.001793, 0.001307, 0.001002, 0.0007977, 0.0006532, 0.000547, 0.0004665, 0.000404, 0.0003544, 0.0003145,
                0.0002818]
mu_c_wa = np.polyfit(t_wa, mu_values_wa, 6)

k_values_wa = [0.561, 0.58, 0.5984, 0.6154, 0.6305, 0.6435, 0.6543, 0.6631, 0.67, 0.6753, 0.6791]
k_c_wa = np.polyfit(t_wa, k_values_wa, 6)

# Data for Alloy 20
tc_temp = f_to_c(np.asarray([70, 100, 150, 200, 250, 300]))  # It gives a warning for some reason when kept as a list
tc_val = thermal_i_to_m([6.9, 6.9, 7.2, 7.5, 7.8, 8.0])
c_tc = np.polyfit(tc_temp, tc_val, 5)
mean_tc_a20 = integrate(c_tc, 90, 25) / (90 - 25)  # W/m K

ys_temp = tc_temp
ys_val = yield_i_to_m([35.0, 35.0, 32.0, 30.9, 30.2, 29.6])
c_ys = np.polyfit(ys_temp, ys_val, 5)
mean_ys_a20 = integrate(c_ys, 90, 25) / (90 - 25)  # Pa

# # Data for hastelloy-C276
# mean_tc_a20 = 9.8  # W/m K

# Pricing Weights and Densities
tube_pw = 0.8 * 26
shell_pw = 0.8 * 4.38
rho_tube = 8100
rho_shell = 8000

# # Pricing Weights and Densities for Hastelloy-C276
# tube_pw = 0.8 * 28
# shell_pw = 0.8 * 28
# rho_tube = 8890
# rho_shell = 8890

# Pipe Data
corrosion_allowance_index = 2
tube_correction = -2
corrosion_allowance = 0.5 * corrosion_allowance_index
tube_th = [[1.2, 1.6, 2, 2],
           [1.6, 2, 2.6, 2.6],
           [1.6, 2, 2.6, 3.2],
           [1.6, 2, 2.6, 3.2],
           [2, 2.6, 3.2, 3.2],
           [2, 2.6, 3.2, 3.2]]
dimensions = {16: tube_th[0][corrosion_allowance_index + tube_correction],
              20: tube_th[1][corrosion_allowance_index + tube_correction],
              25: tube_th[2][corrosion_allowance_index + tube_correction],
              30: tube_th[3][corrosion_allowance_index + tube_correction],
              38: tube_th[4][corrosion_allowance_index + tube_correction],
              50: tube_th[5][corrosion_allowance_index + tube_correction]}  # Outer Diameter : Thickness
outer_ds = list(dimensions.keys())
inner_ds = [val - (dimensions[val] * 2) for val in outer_ds]

# Tube Sheet Data
min_ts_thickness = [19.1, 19.1, 19.1, 22.2, 25.4, 31.8]  # From B-7.1.1 in TEMA 10th edition
ts_thickness = [(tst + corrosion_allowance) / 1000 for tst in min_ts_thickness]

# Fouling Resistances
ri = 0.0002
ro_max = 0.0003  # Cooling Tower Water
ro_min = 0.00017  # Cooling Tower Water. Would be 0.0001 for river water
ro_selected = (ro_min + ro_max) / 2
hf_selected = 1 / ro_selected
# Since the fouling resistances are low, only the >1000 values are included
if hf_selected < 1000 or hf_selected > 6000:
    raise Exception("Fouling resistance is too high or low.")
f_r_t_6 = [1.12, 1.38, 1.55]
f_r_t_2 = [1.37, 2.31, 2.96]
f_r_t_1 = [1.64, 3.44, 4.77]

f_r_l_6 = [1.06, 1.20, 1.28]
f_r_l_2 = [1.19, 1.44, 1.55]
f_r_l_1 = [1.32, 1.99, 2.38]
if hf_selected > 2000:
    fouled_ratio_t = [np.interp(hf_selected, [2000, 6000], [f_r_t_2[_index], f_r_t_6[_index]]) for _index in range(3)]
    fouled_ratio_l = [np.interp(hf_selected, [2000, 6000], [f_r_l_2[_index], f_r_l_6[_index]]) for _index in range(3)]
else:
    fouled_ratio_t = [np.interp(hf_selected, [1000, 2000], [f_r_t_1[_index], f_r_t_2[_index]]) for _index in range(3)]
    fouled_ratio_l = [np.interp(hf_selected, [1000, 2000], [f_r_l_1[_index], f_r_l_2[_index]]) for _index in range(3)]

sd_bs_ratio = [1, 2, 5]

# shell-side coefficient jh from graph
jh_re_s = np.log10(np.asarray([1.75, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                               20, 30, 40, 50, 60, 70, 80, 90, 100,
                               200, 300, 400, 500, 600, 700, 800, 900, 1000]) * 1000)
jh_pixels_s = [30, 35, 51, 63, 71, 79, 84, 89, 95, 98, 124, 136, 145,
               152, 157, 162, 165, 169, 172, 191, 200, 206,
               211, 215, 218, 220, 222, 224]
jh_values_s = [((278 - pixel) / (278 - 75)) - 3 for pixel in jh_pixels_s]
jh_c_s = np.polyfit(jh_re_s, jh_values_s, 6)
jh_flat_s = pow(10, jh_values_s[0])

# Baffle thickness data
baffle_th = [[3.2, 4.8, 6.4, 9.5, 9.5],
             [4.8, 6.4, 9.5, 9.5, 12.7],
             [6.4, 7.9, 9.5, 12.7, 15.9],
             [6.4, 9.5, 12.7, 15.9, 15.9],
             [9.5, 12.7, 15.9, 19.1, 19.1]]

# Fn values from graph
fn_ncv = [5, 10, 15, 20, 25, 30, 35]
fn_pixels = [312, 239, 193, 166, 150, 137, 128]
fn_values = [(((429 - pixel) / (429 - 26)) / 5) + 0.9 for pixel in fn_pixels]
fn_c = np.polyfit(fn_ncv, fn_values, 6)

# Fw values from graph
# It was found that many more data points were required.

# rw_1 = 0.1 * ((84-65)/(163-65))
rw_2 = 0.8 + 0.1 * ((917 - 847) / (163 - 65))

# Using a modified version of the image in Coulson and Richardson
fw_c = itp.get_fw_c()

# Fb values from (a modified version of) the graph
fb_c_l = itp.get_fb_c_l()
fb_c_h = itp.get_fb_c_h()

bl_c, bl_ar_max = itp.get_bl_c()
p_bl_ar, p_bl_val, p_bl_ar_max = itp.get_p_bl_c()

r_c, r_x, r_y = itp.get_r_c()
theta_c = itp.get_theta_c()

# Maximum unsupported spans from TEMA
lm_od = [6.4, 9.5, 12.7, 15.9, 19.1, 22.2, 25.4, 31.8, 38.1, 50.8]
lm_values = [660, 889, 1118, 1321, 1524, 1753, 1880, 2235, 2540, 3175]

# Data for jf friction factor from the graph
jf_x, jf_y = itp.get_jf_data()

# Data for jf friction factor for pressure from the graph
p_jf_x, p_jf_y = itp.get_p_jf_data()

# Data for the Fb factor for pressure from the graph
p_fb_c_l, p_fb_l_max = itp.get_p_fb_c_l()
p_fb_c_h, p_fb_h_max = itp.get_p_fb_c_h()

# Tags for data output
out_tags = ["Tube HTC", "Tube Velocity", "Tubes Per Pass", "Shell Inner Diameter",
            "Overall HTC", "HT Area", "Tube-Side Reynolds Number", "Baffle-Shell Clearance", "Bundle Diameter",
            "Tube Wall Resistance", "Tube Sheet Thickness", "Head Radial Wall Thickness"]
hs_tags = ["Shell HTC", "Total Shell F", "Ideal Shell HTC",
           "Tube Row CF (Fn)", "Window CF (Fw)", "Bypass CF (Fb)", "Leakage CF (Fl)",
           "Cross-Flow Area (As)", "Reynolds Number", "Number of Sealing Strips", "Shell Velocity", "Constriction Rows",
           "Tubes in Window Zone", "Ideal Baffle Ratio"]

max_iter = 128
