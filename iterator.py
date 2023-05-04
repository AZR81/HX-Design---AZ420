from math import ceil
from time import time

import matplotlib as mpl
import matplotlib.pyplot as plt

import mechanical_design as md
from constants import *

mpl.rc('font', family='Times New Roman')
session_id = int(time())
print(session_id)

outlet_temp_wa = 43
def recalculate_water_data():
    _bulk_temp_wa = (outlet_temp_wa + inlet_temp_wa) / 2
    _mean_cp_wa = integrate(c_wa, outlet_temp_wa, inlet_temp_wa) / (outlet_temp_wa - inlet_temp_wa)
    _rho_wa = integrate(rho_c_wa, outlet_temp_wa, inlet_temp_wa) / (outlet_temp_wa - inlet_temp_wa)
    _mu_wa = integrate(mu_c_wa, outlet_temp_wa, inlet_temp_wa) / (outlet_temp_wa - inlet_temp_wa)
    _mu_wa_bulk = evaluate_polynomial(mu_c_wa, _bulk_temp_wa)
    _k_wa = integrate(k_c_wa, outlet_temp_wa, inlet_temp_wa) / (outlet_temp_wa - inlet_temp_wa)
    return _bulk_temp_wa, _mean_cp_wa, _rho_wa, _mu_wa, _mu_wa_bulk, _k_wa


bulk_temp_wa, mean_cp_wa, rho_wa, mu_wa, mu_wa_bulk, k_wa = recalculate_water_data()
print(bulk_temp_wa, mean_cp_wa, rho_wa, mu_wa, mu_wa_bulk, k_wa, outlet_temp_wa)

print("Loaded all data.")


def get_baffle_th(shell_id, u_length):
    if shell_id < 152:
        raise Exception(f"Shell diameter is too small({shell_id}).")
    elif shell_id < 356:
        bti = 0
    elif shell_id < 711:
        bti = 1
    elif shell_id < 965:
        bti = 2
    elif shell_id < 1524:
        bti = 3
    elif shell_id < 2540:
        bti = 4
    else:
        raise Exception(f"Shell diameter is too large({shell_id}).")
    bth_data = baffle_th[bti]
    if u_length < 610:
        return bth_data[0]
    elif u_length < 914:
        return bth_data[1]
    elif u_length < 1219:
        return bth_data[2]
    elif u_length < 1524:
        return bth_data[3]
    return bth_data[4]


def get_shell_th_limit(shell_id):
    if shell_id < 200:
        return 3.2
    elif shell_id < 300:
        return 3.2
    elif shell_id < 580:
        return 3.2
    elif shell_id < 740:
        return 4.8
    elif shell_id < 990:
        return 6.4
    elif shell_id < 1520:
        return 6.4
    elif shell_id < 2030:
        return 7.9
    elif shell_id < 2540:
        return 9.5
    raise Exception(f"Shell diameter is too large ({shell_id}).")


def get_shell_th(shell_id, use_sheet=False, nozzle=False, design_pressure=0.0, max_tension=0):
    # Using the dimensions for pipes following ANSI/ASME B36.19M standards
    if nozzle:
        ideal_th = md.get_cylindrical_th(design_pressure, shell_id, max_tension) + corrosion_allowance
    else:
        ideal_th = get_shell_th_limit(shell_id)
    if use_sheet:
        return None, None, shell_id, ideal_th
    idea_od = shell_id + ideal_th
    nearest_od = []
    for sod in pipe_sizes.keys():
        if sod >= idea_od:
            nearest_od.append(sod)
    if not nearest_od:
        return None, None, shell_id, ideal_th
    sid_difference = 10000000
    best_size = None
    best_th = None
    output_sid = shell_id
    output_th = ideal_th
    stop_flag = False
    for shell_od in nearest_od:
        if stop_flag:
            break
        pipe_size_keys = list(pipe_sizes[shell_od][1].keys())
        sth = [shell_od - (2*val) if val is not None else shell_od for val in pipe_sizes[shell_od][1].values()]
        for sth_index, sid_val in enumerate(sth):
            th_val = pipe_sizes[shell_od][1][pipe_size_keys[sth_index]]
            if nozzle and (sid_val >= shell_id) and th_val is not None and \
                    (md.get_cylindrical_th(design_pressure, shell_id, max_tension) < th_val):
                best_size = pipe_sizes[shell_od][0]
                best_th = pipe_size_keys[sth_index]
                output_sid = sid_val
                output_th = th_val
                stop_flag = True
                break
            elif not nozzle:
                if (sid_val - shell_id > 0) and (
                        (sid_val - shell_id) < sid_difference) and th_val is not None and th_val >= ideal_th:
                    sid_difference = sid_val - shell_id
                    best_size = pipe_sizes[shell_od][0]
                    best_th = pipe_size_keys[sth_index]
                    output_sid = shell_id + sid_difference
                    output_th = th_val
    return best_size, best_th, output_sid, output_th


def calculate_water_flow_rate(outlet_temp):
    return (-85000 * mean_cp_ipa * (inlet_temp_ipa - outlet_temp_ipa)) / \
           (mean_cp_wa * (inlet_temp_wa - outlet_temp))  # 178064.1438634734


def plot_water_flow_rate(minimum_temp=41):
    xvs = np.asarray(np.linspace(minimum_temp, 45, 10 * (45 - minimum_temp + 1)))
    yvs = np.asarray([calculate_water_flow_rate(tmp) for tmp in xvs])
    plt.plot(xvs, yvs)
    plt.xlim([minimum_temp, 45])
    plt.ylim([160000, yvs[0]])
    plt.xticks([35, 37, 39, 41, 43, 45], [35, 37, 39, 41, 43, 45], fontname="Times New Roman")
    plt.yticks([160000, 180000, 200000, 220000, 240000, 260000, 280000, 300000, 320000],
               [160, 180, 200, 220, 240, 260, 280, 300, 320], fontname="Times New Roman")
    plt.xlabel("Water Outlet Temperature ($^\circ$C)", fontname="Times New Roman")
    plt.ylabel("Water Mass Flow Rate ($10^3$ kg hr$^{-1}$)", fontname="Times New Roman")
    plt.show()


def calculate_duty():
    return (mean_cp_ipa * (inlet_temp_ipa - outlet_temp_ipa) * 85000) / (60 * 60)


def calculate_area(Q, U, F, LMTD):
    return Q / (U * F * LMTD)


def calculate_log_mean(dt1, dt2):
    return (dt1 - dt2) / np.log(dt1 / dt2)


def get_f(t1, t2, T1, T2):
    R = (T1 - T2) / (t2 - t1)
    S = (t2 - t1) / (T1 - t1)
    num = np.sqrt((R ** 2) + 1) * np.log((1 - S) / (1 - R * S))
    den = (R - 1) * np.log((2 - S * (R + 1 - np.sqrt(1 + (R ** 2)))) / (2 - S * (R + 1 + np.sqrt(1 + (R ** 2)))))

    return num / den


def new_get_f(t1, t2, T1, T2, N):
    P = (T2 - T1) / (t1 - T1)
    R = (t1 - t2) / (T2 - T1)
    if R == 1:
        W1 = (N - (N * P)) / (N - (N * P) + P)
        sq2 = np.sqrt(2)
        return sq2 * (((1 - W1) / W1) / np.log(((W1 / (1 - W1)) + (1 / sq2)) / ((W1 / (1 - W1)) - (1 / sq2))))

    S = np.sqrt((R ** 2) + 1) / (R - 1)
    W = pow(((1 - P * R) / (1 - P)), 1 / N)
    return S * np.log(W) / np.log((1 + W - S + S * W) / (1 + W + S - S * W))


def plot_f_factor(outlet_temp=43, max_shells=9):
    xvs = list(range(1, max_shells + 1))
    yvs = [new_get_f(90, 40, 25, outlet_temp, vl) for vl in xvs]
    plt.plot(xvs, yvs)
    plt.show()


def get_dp_i(velocity, velocity_n, _length, inside_diameter, pass_count, reynolds, hi, u):  # m
    if reynolds < 10 or reynolds > pow(10, 6):
        raise Exception(f"Reynolds not in range 10 to 10^6 ({reynolds}).")
    inside_diameter = inside_diameter / 1000
    velocity_head = (rho_ipa * (velocity ** 2)) / 2
    nozzle_velocity_head = (rho_ipa * (velocity_n ** 2)) / 2
    jf = pow(10, np.interp(np.log10(reynolds), jf_x, jf_y))
    if reynolds < 2100:
        m = 0.25
    else:
        m = 0.14

    t_wall = ((u / hi) * (bulk_temp_wa - bulk_temp_ipa)) + bulk_temp_ipa
    viscosity_factor = mu_ipa_bulk / evaluate_polynomial(mu_c_ipa, t_wall)
    tube_dp = pass_count * velocity_head * ((8 * jf * (_length / inside_diameter) * pow(viscosity_factor, -m)) + 2.5)
    nozzle_dp = nozzle_velocity_head * (0.5 + 1 + 1.5)

    return nozzle_dp + tube_dp


def get_hi(velocity, di):
    # HMT Fundamentals Notes, at bulk temperature
    reynolds = (rho_ipa_bulk * velocity * (di / 1000)) / mu_ipa_bulk
    prandtl = (mean_cp_ipa_bulk * mu_ipa_bulk) / k_ipa_bulk
    if 2200 <= reynolds <= (5 * (10 ** 5)):
        # Blasius
        cf = 0.079 * pow(reynolds, -0.25)
        numerator = (cf / 2) * (reynolds - 1000) * prandtl
        denominator = 1 + 12.7 * np.sqrt(cf / 2) * (pow(prandtl, 0.67) - 1)
        nusselt = numerator / denominator
        return ((nusselt * k_ipa_bulk) / (di / 1000)), reynolds
    raise Exception(f"Reynolds is not in range ({reynolds}).")


def get_jh(reynolds):
    if reynolds < 10:
        raise Exception(f"Reynolds is not in range ({reynolds}).")
    elif 10 <= reynolds <= 90:
        return pow(10, np.log10(reynolds) * ((((361 - 204) / (361 - 155)) - 2) - (((361 - 74) / (361 - 155)) - 2)))
    elif 90 < reynolds <= 1000:
        adjustment1 = np.log10(90) * ((((361 - 204) / (361 - 155)) - 2) - (((361 - 74) / (361 - 155)) - 2))
        adjustment2 = np.log10(90) * ((((361 - 313) / (361 - 155)) - 2) - (((361 - 204) / (361 - 155)) - 2))
        adjustment3 = adjustment1 - adjustment2
        return pow(10, adjustment3 + np.log10(reynolds) * (
                (((361 - 313) / (361 - 155)) - 2) - (((361 - 204) / (361 - 155)) - 2)))
    elif 1000 < reynolds <= 1750:
        return jh_flat_s
    else:
        return pow(10, evaluate_polynomial(jh_c_s, np.log10(reynolds)))


def get_gs(As):
    return (calculate_water_flow_rate(outlet_temp_wa) / 3600) / As


def get_hoc(d_out, reynolds, vcf):
    prandtl = (mu_wa * mean_cp_wa) / k_wa
    hoc1 = (k_wa / d_out) * get_jh(reynolds) * reynolds * pow(prandtl, 1 / 3) * pow(vcf, 0.14)
    return hoc1


def get_fn(ncv, reynolds):  # needs mm
    if reynolds > 2000:
        if 5 <= ncv <= 35:
            return evaluate_polynomial(fn_c, ncv)
        raise Exception(f"Ncv is not in range ({ncv}).")
    elif 100 < reynolds <= 2000:
        return 1
    raise Exception(f"Reynolds is not in range ({reynolds}).")


def get_fw(rw):  # needs mm
    if 0 <= rw <= rw_2:
        return evaluate_polynomial(fw_c, rw)
    raise Exception(f"Rw is not in range ({rw}).")


def get_fb(reynolds, ds_in, db, baffle_s, cf_area, strips, constrictions):
    if reynolds < 100:
        alpha = 1.5
    else:
        alpha = 1.35

    Ab = (ds_in - db) * baffle_s * pow(10, -6)
    a_ratio = Ab / cf_area
    if strips != 0:
        if strips > (constrictions / 2):
            raise Exception(f"Strips and constrictions condition not met ({strips}, {constrictions}).")
        n_ratio = pow((2 * strips) / constrictions, 1 / 3)
        return np.exp(-alpha * a_ratio * (1 - n_ratio))
    if 0 <= a_ratio <= 0.4:
        if reynolds < 100:
            return evaluate_polynomial(fb_c_l, a_ratio)
        return evaluate_polynomial(fb_c_h, a_ratio)
    raise Exception(f"Clearance and CF Area conditions not met ({Ab}, {cf_area}, {a_ratio}).")


def get_fl(tb_c, sb_c, cf_area, tc, tiw, br, sid, heat_transfer=True):
    tb_area = (tb_c / 2) * 20 * np.pi * (tc - tiw) * pow(10, -6)
    if 0.15 <= br <= 0.45:
        theta = evaluate_polynomial(theta_c, br)
    else:
        raise Exception(f"Invalid value of baffle ratio ({br}).")
    sb_area = (sb_c / 2) * sid * (2 * np.pi - theta) * pow(10, -6)
    leakage_area = tb_area + sb_area
    if heat_transfer:
        bl_max = bl_ar_max
    else:
        bl_max = p_bl_ar_max
    if 0 <= (leakage_area / cf_area) <= bl_max:
        if heat_transfer:
            beta = evaluate_polynomial(bl_c, leakage_area / cf_area)
        else:
            beta = np.interp(leakage_area / cf_area, p_bl_ar, p_bl_val)
        return 1 - beta * ((tb_area + (2 * sb_area)) / leakage_area)
    raise Exception(f"Invalid area ratio ({leakage_area / cf_area}, {leakage_area}, {cf_area}, {bl_max}).")


def get_hs(pitch, d_out, ds_in, baffle_s, tubes, br, strips_ratio, tb_clearance, sb_clearance, db, vcf=1.0):
    As = (((pitch - d_out) * ds_in * baffle_s) / pitch) / 1000000
    reynolds = (get_gs(As) * (d_out / 1000)) / mu_wa

    baffle_cut_height = br * ds_in
    baffle_tips_distance = ds_in - 2 * baffle_cut_height
    constrictions = round(baffle_tips_distance / pitch)  # Square pitch
    if constrictions % 2 == 1:
        constrictions -= 1
    ideal_baffle_cut_height = (ds_in - constrictions * pitch) / 2
    ideal_baffle_ratio = round(ideal_baffle_cut_height / ds_in, 2)
    Hb = db / 2 - ds_in * (0.5 - br)
    bc = Hb / db
    if 0.15 <= bc <= 0.45:
        Rap = evaluate_polynomial(r_c, bc)
    else:
        raise Exception(f"Invalid value of baffle ratio ({bc}).")
    tiw = ceil(tubes * Rap)
    rw = (2 * tiw) / tubes
    strips = int(strips_ratio * constrictions)

    # Required Data:

    # Pitch distance (1.25 d_out?)
    # outside pipe diameter
    # inside shell diameter
    # Baffle spacing

    # Outside pipe diameter
    # Reynolds

    # number of tube rows between the baffle tips
    # Reynolds
    # Nc` (will most likely raise an exception)

    # ratio of the number of tubes in the window zones to the total number in the bundle

    # reynolds
    # clearance area between the bundle and the shell
    # maximum area for cross-flow, equation 12.21,
    # number of sealing strips encountered by the bypass stream in the cross-flow zone
    # the number of constrictions, tube rows, encountered in the cross-flow section.

    # tube to baffle clearance area, per baffle
    # shell to baffle clearance area, per baffle

    hoc = get_hoc(d_out / 1000, reynolds, vcf)  # m (converted here)
    fn = get_fn(constrictions, reynolds)  # mm
    fw = get_fw(rw)  # mm
    fb = get_fb(reynolds, ds_in, db, baffle_s, As, strips, constrictions)  # mm
    fl = get_fl(tb_clearance, sb_clearance, As, tubes, tiw, br, ds_in)  # mm

    return [hoc * fn * fw * fb * fl, fn * fw * fb * fl, hoc, fn, fw, fb, fl, As, reynolds, strips,
            (get_gs(As) / rho_wa),
            constrictions, tiw, ideal_baffle_ratio]


def get_fb_prime(reynolds, ds_in, db, baffle_s, cf_area, strips, constrictions):
    if reynolds < 100:
        alpha = 5
    else:
        alpha = 4

    Ab = (ds_in - db) * baffle_s * pow(10, -6)
    a_ratio = Ab / cf_area
    if strips != 0:
        if strips > (constrictions / 2):
            raise Exception(f"Strips and constrictions condition not met ({strips}, {constrictions}).")
        n_ratio = pow((2 * strips) / constrictions, 1 / 3)
        return np.exp(-alpha * a_ratio * (1 - n_ratio))

    if reynolds < 100:
        if 0 <= a_ratio <= p_fb_l_max:
            return evaluate_polynomial(p_fb_c_l, a_ratio)
        raise Exception(f"A ratio was outside the range ({p_fb_l_max}, {a_ratio}, {Ab}, {cf_area}).")

    if reynolds >= 100:
        if 0 <= a_ratio <= p_fb_h_max:
            return evaluate_polynomial(p_fb_c_h, a_ratio)
        raise Exception(f"A ratio was outside the range ({p_fb_h_max}, {a_ratio}, {Ab}, {cf_area}).")
    raise Exception(f"Clearance and CF Area conditions not met ({Ab}, {cf_area}, {a_ratio}).")


def get_dp_cf(As, u, hs, reynolds, ds_in, db, baffle_s, cf_area, strips, constrictions, tb_clearance, sb_clearance,
              tubes, tiw, br):
    if reynolds < 10 or reynolds > pow(10, 6):
        raise Exception(f"Reynolds not in range 10 to 10^6 ({reynolds}).")
    if reynolds < 1000:
        print("Reynolds is below 1000. Pressure value error might be large.")
    us = get_gs(As) / rho_wa
    velocity_head = (rho_wa * (us ** 2)) / 2

    t_wall = ((u / hs) * (bulk_temp_wa - bulk_temp_ipa)) + bulk_temp_wa
    viscosity_factor = mu_wa_bulk / evaluate_polynomial(mu_c_wa, t_wall)
    jf = pow(10, np.interp(np.log10(reynolds), p_jf_y, p_jf_y))
    dp_ideal = 8 * jf * constrictions * velocity_head * pow(viscosity_factor, -0.14)  # replaced ncv with constrictions
    fb_prime = get_fb_prime(reynolds, ds_in, db, baffle_s, cf_area, strips, constrictions)
    fl_prime = get_fl(tb_clearance, sb_clearance, As, tubes, tiw, br, ds_in, heat_transfer=False)

    return dp_ideal * fb_prime * fl_prime, dp_ideal, fb_prime, fl_prime


def get_dp_window(pitch, bundle_d, shell_d, baffle_cut, dt_out, velocity_s, fl_prime, nw):
    hb = (bundle_d / 2) - (shell_d * (0.5 - baffle_cut))  # height from the baffle chord to the top of the bundle
    nwv = ceil(hb / pitch)  # Number of restrictions for cross flow in window zone
    if 0.15 <= baffle_cut <= 0.45:
        ra = evaluate_polynomial(r_c, baffle_cut)
    else:
        raise Exception(f"Baffle cut value not in range 0.15 to 0.45 ({baffle_cut}).")

    # Window area - tubes area
    aw = ((np.pi * (shell_d ** 2) * 0.25 * ra) - (nw * np.pi * 0.25 * (dt_out ** 2))) * pow(10, -6)
    ws = calculate_water_flow_rate(outlet_temp_wa) / (60 * 60)  # Water mass flow rate
    uw = ws / (aw * rho_wa)  # Velocity in window zone
    uz = np.sqrt(uw * velocity_s)
    velocity_head = (rho_wa * (uz ** 2)) / 2
    return fl_prime * (2 + (0.6 * nwv)) * velocity_head, nwv


def get_dp_s(As, u, hs, reynolds, ds_in, db, baffle_s, strips, constrictions, tb_clearance,
             sb_clearance, tubes, tiw, br, pitch, dt_out, velocity_s, baffle_count, nozzle_radius, bsg):
    # Cross flow zone pressure drop
    dp_cf, dp_ideal, fb_prime, fl_prime = get_dp_cf(As, u, hs, reynolds, ds_in, db, baffle_s, As, strips,
                                                    constrictions, tb_clearance, sb_clearance, tubes, tiw, br)
    # Window zone pressure drop
    dp_window, nwv = get_dp_window(pitch, db, ds_in, br, dt_out, velocity_s, fl_prime, tiw)
    # End zone pressure drop
    dp_end = dp_ideal * fb_prime * ((nwv + constrictions) / constrictions)

    # Nozzle pressure drop
    nozzle_area = np.pi * nozzle_radius * nozzle_radius
    tube_gap_area = (baffle_s / 1000) * (((pitch - dt_out) / 1000) * tube_gap_count + (bsg / 500))
    chosen_area = min([nozzle_area, tube_gap_area])
    volumetric_flow_rate = (calculate_water_flow_rate(outlet_temp_wa) / 3600) / rho_wa
    nozzle_velocity = volumetric_flow_rate / chosen_area
    velocity_head = (rho_wa * nozzle_velocity * nozzle_velocity) / 2
    nozzle_loss = (1.5 + 0.5) * velocity_head
    return 2 * dp_end + dp_cf * (baffle_count - 1) + baffle_count * dp_window + nozzle_loss


def get_minimum_radius(mass_flow_rate, density, max_head, max_velocity=0.0):
    if max_velocity > 0:
        max_head = density * max_velocity * max_velocity
    volumetric_flow_rate = mass_flow_rate / density
    return pow((density * pow(volumetric_flow_rate, 2)) / (max_head * np.pi * np.pi), 1 / 4)


def htc_iterator(trial_u, shell_passes, min_tube_sheet_thickness, tube_outer_diameter, bundle_shell_gap,
                 tube_inner_diameter,
                 tube_pitch, baffle_count, baffle_opening_ratio, strips_to_constrictions_ratio,
                 tube_to_baffle_clearance,
                 shell_fouling_resistance, length_total, use_sheet=True):
    U = trial_u
    last_U = 0.1
    out = []
    converged = True
    iter_count = 0
    K1 = 0.156
    n1 = 2.291
    previous_ts_th = min_tube_sheet_thickness
    viscosity_correction_factor = 1.0
    while abs(last_U - U) > 0.1 and iter_count < max_iter:
        iter_count += 1
        try:
            last_U = U
            area = calculate_area(calculate_duty(),
                                  U,
                                  new_get_f(inlet_temp_ipa, outlet_temp_ipa, inlet_temp_wa, outlet_temp_wa,
                                            shell_passes),
                                  calculate_log_mean(inlet_temp_ipa - outlet_temp_wa, outlet_temp_ipa - inlet_temp_wa))

            a_t = (length_total - 2 * previous_ts_th) * np.pi * (tube_outer_diameter / 1000) * 2 * shell_passes
            tube_count = ceil((area / a_t))
            No = tube_count * 2
            bundle_diameter = ceil(tube_outer_diameter * pow((No / K1), (1 / n1)))

            _, _, shell_inner_diameter, _ = get_shell_th(bundle_diameter + bundle_shell_gap, use_sheet=use_sheet)
            if shell_inner_diameter > 635:
                baffle_to_shell_clearance = 4.8
            else:
                baffle_to_shell_clearance = 3.2

            vel = vfr_ipa / (tube_count * ((tube_inner_diameter / 2000) ** 2) * np.pi)
            hi, re_in = get_hi(vel, tube_inner_diameter)

            # Left a 10 mm gap for the floating tube sheet
            wall_th = np.round(
                (md.get_cylindrical_th(18 * 1.1, shell_inner_diameter - 10, md.a20_design_stress, alt=True,
                                       cr=corrosion_allowance)), 1)
            tube_sheet_thickness = np.round(
                (md.get_ts_thickness(shell_inner_diameter, wall_th, 18 * 1.1, md.a20_design_stress,
                                     bundle_diameter, min_tube_sheet_thickness)), 4)
            previous_ts_th = tube_sheet_thickness
            _epl = (length_total - 2 * tube_sheet_thickness)
            baffle_pitch = int(_epl * 1000) / (baffle_count + 1)
            hs_val = get_hs(tube_pitch,
                            tube_outer_diameter,
                            shell_inner_diameter,
                            baffle_pitch,
                            tube_count * 2,
                            baffle_opening_ratio,
                            strips_to_constrictions_ratio,
                            tube_to_baffle_clearance,
                            baffle_to_shell_clearance,
                            bundle_diameter,
                            viscosity_correction_factor)
            hs = hs_val[0]

            wall_temperature = ((U / hs) * (bulk_temp_wa - bulk_temp_ipa)) + bulk_temp_wa
            viscosity_correction_factor = mu_wa_bulk / evaluate_polynomial(mu_c_wa, wall_temperature)

            wall_res = ((tube_outer_diameter / 1000) * np.log(tube_outer_diameter / tube_inner_diameter)) / (
                    2 * mean_tc_a20)
            urev = (1 / hs +
                    shell_fouling_resistance +
                    wall_res +
                    (tube_outer_diameter / tube_inner_diameter) * ri +
                    (tube_outer_diameter / tube_inner_diameter) * (1 / hi))
            U = 1 / urev
            out = [hi, vel, hs_val, tube_count, shell_inner_diameter, U, area, re_in,
                   baffle_to_shell_clearance, bundle_diameter, wall_res, tube_sheet_thickness, wall_th]

        except Exception as exc:
            out.append(exc)
            last_U = U
            converged = False

    if iter_count == (max_iter - 1):
        converged = False
        return out, converged

    # Oversize the unit
    previous_vcf = -1
    iter_c_111 = 0
    try:
        tube_count = ceil(out[3] * 1.2)
        No = tube_count * 2
        bundle_diameter = ceil(tube_outer_diameter * pow((No / K1), (1 / n1)))

        _, _, shell_inner_diameter, _ = get_shell_th(bundle_diameter + bundle_shell_gap, use_sheet=use_sheet)
        if shell_inner_diameter > 635:
            baffle_to_shell_clearance = 4.8
        else:
            baffle_to_shell_clearance = 3.2

        vel = vfr_ipa / (tube_count * ((tube_inner_diameter / 2000) ** 2) * np.pi)
        wall_th = np.round((md.get_cylindrical_th(18 * 1.1, shell_inner_diameter - 10, md.a20_design_stress, alt=True,
                                                  cr=corrosion_allowance)), 1)
        tube_sheet_thickness = np.round(
            (md.get_ts_thickness(shell_inner_diameter, wall_th, 18 * 1.1, md.a20_design_stress,
                                 bundle_diameter, min_tube_sheet_thickness)), 4)
        _epl = (length_total - 2 * tube_sheet_thickness)
        baffle_pitch = int(_epl * 1000) / (baffle_count + 1)
        hi, re_in = get_hi(vel, tube_inner_diameter)
        hs_val = 0
        wall_res = 0
        while abs(viscosity_correction_factor - previous_vcf) > 0.00001:  # U is very sensitive to VCF
            iter_c_111 += 1
            previous_vcf = viscosity_correction_factor
            hs_val = get_hs(tube_pitch,
                            tube_outer_diameter,
                            shell_inner_diameter,
                            baffle_pitch,
                            tube_count * 2,
                            baffle_opening_ratio,
                            strips_to_constrictions_ratio,
                            tube_to_baffle_clearance,
                            baffle_to_shell_clearance,
                            bundle_diameter,
                            viscosity_correction_factor)
            hs = hs_val[0]

            wall_temperature = ((U / hs) * (bulk_temp_wa - bulk_temp_ipa)) + bulk_temp_wa
            viscosity_correction_factor = mu_wa_bulk / evaluate_polynomial(mu_c_wa, wall_temperature)

            wall_res = ((tube_outer_diameter / 1000) * np.log(tube_outer_diameter / tube_inner_diameter)) / (
                    2 * mean_tc_a20)
            urev = (1 / hs +
                    shell_fouling_resistance +
                    wall_res +
                    (tube_outer_diameter / tube_inner_diameter) * ri +
                    (tube_outer_diameter / tube_inner_diameter) * (1 / hi))
            U = 1 / urev
        out = [hi, vel, hs_val, tube_count, shell_inner_diameter, U, out[6], re_in,
               baffle_to_shell_clearance, bundle_diameter, wall_res, tube_sheet_thickness, wall_th]

    except Exception as exc:
        out.append(exc)
        converged = False
    return out, converged


def get_price(shell_id, shell_t, shell_l, shell_p, end_id, end_t, end_plate_t,
              tube_od, tube_t, tube_l, tube_p_tot, tube_n,
              baffle_t, baffle_c, baffle_n,
              tube_sheet_thickness, head_t, io_length, io_thickness, strip_count, epl, t_pitch):
    shell_id = shell_id / 1000
    shell_t = shell_t / 1000
    end_id = end_id / 1000
    end_plate_t = end_plate_t / 1000
    tube_od = tube_od / 1000
    tube_t = tube_t / 1000
    baffle_t = baffle_t / 1000
    head_t = head_t / 1000
    io_length = io_length / 1000
    io_thickness = io_thickness / 1000
    tube_or = tube_od / 2
    tube_ir = tube_or - tube_t
    shell_ir = shell_id / 2
    shell_or = shell_ir + shell_t
    head_ir = shell_ir - head_t

    strip_volume = shell_p * strip_count * epl * t_pitch * baffle_t
    strip_adj_weight = rho_shell * strip_volume * shell_pw

    tube_volume = np.pi * (pow(tube_or, 2) - pow(tube_ir, 2)) * tube_l * tube_p_tot * tube_n
    tube_adj_weight = rho_tube * tube_volume * tube_pw
    # Since the inlet and outlet areas and tube sheets are made of the same materials as the tubes,
    # the same density and price weight is used for them

    # The gap between the tube sheet and the shell has been taken as being negligible
    t_sheet_volume = 2 * tube_sheet_thickness * (pow(shell_ir, 2) - (tube_n * 2 * pow(tube_or, 2))) * np.pi
    t_sheet_adj_weight = rho_tube * t_sheet_volume * tube_pw * shell_p

    # Inlet and outlet area has been approximated as a capped cylinder
    io_volume = io_length * np.pi * (pow(shell_or, 2) - pow(shell_ir, 2)) + np.pi * pow(shell_or, 2) * io_thickness
    io_adj_weight = rho_tube * io_volume * tube_pw * shell_p

    # Floating head approximated as a cylinder. Very much an approximation but good enough
    head_volume = np.pi * (pow(shell_ir, 2) - pow(head_ir, 2)) * shell_id * 0.6
    head_adj_weight = rho_tube * head_volume * tube_pw * shell_p

    # Baffles and end caps are made of the same material as the shell body
    shell_volume = np.pi * (pow(shell_or, 2) - pow(shell_ir, 2)) * shell_l
    adj_shell_weight = rho_shell * shell_volume * shell_pw * shell_p

    # End zone
    end_plate_volume = np.pi * pow(shell_or, 2) * end_plate_t
    end_adj_weight = (end_plate_volume * shell_p * rho_tube * tube_pw)

    theta = evaluate_polynomial(theta_c, baffle_c)
    baffle_area = baffle_n * baffle_t * (
            (np.pi * pow(shell_ir, 2)) - (0.5 * pow(shell_ir, 2) * (theta - np.sin(theta))))
    baffle_adj_weight = baffle_area * shell_p * rho_shell * shell_pw
    return (tube_adj_weight + t_sheet_adj_weight + io_adj_weight + head_adj_weight + adj_shell_weight +
            end_adj_weight + baffle_adj_weight + strip_adj_weight)


def get_pumping_cost(pressure, mass_flow_rate_hour, density):
    flow_rate = mass_flow_rate_hour / (60 * 60)
    head = (pressure * pow(10, 5)) / (9.81 * density)
    power = (flow_rate * 9.81 * head) / 0.9  # 90% efficient based on Debem's website
    return (power / 1000) * 20 * 365.25 * 24 * 0.185


bps = []
ohtcs = []
nts = []
dpts = []
dpss = []
sdis = []


def trial_run(length_total):
    minimum_clearance = 6.4
    bs_gap = 58

    # Baffle
    baffle_pitch = (length_total * 1000) / 20
    baffle_opening_ratio = 0.25
    strips_to_constrictions_ratio = 0
    tube_to_baffle_clearance = 0.8

    shell_passes = 2

    with open(f"./Output/trial_{session_id}.txt", "a") as file:
        file.write(f"Minimum Clearance:{minimum_clearance}\n")
        file.write(f"Bundle-Shell Gap:{bs_gap}\n")
        file.write(f"Baffle Pitch:{baffle_pitch}\n")
        file.write(f"Baffle Opening Ratio:{baffle_opening_ratio}\n")
        file.write(f"Strips to Constrictions Ratio:{strips_to_constrictions_ratio}\n")
        file.write(f"Tube to Baffle Clearance:{tube_to_baffle_clearance}\n")
        file.write(f"Shell Passes:{shell_passes}\n")
        file.write(f"Shell Fouling Resistance:{(ro_min + ro_max) / 2}\n")
        file.write("\n")

    for tube_outer_diameter in outer_ds:
        pitch_multiplier = 1.5
        tube_pitch = tube_outer_diameter * pitch_multiplier
        if (tube_pitch - tube_outer_diameter) < minimum_clearance:
            tube_pitch = tube_outer_diameter + minimum_clearance

        index = outer_ds.index(tube_outer_diameter)
        iter_results = htc_iterator(500, shell_passes, ts_thickness[index], tube_outer_diameter, bs_gap,
                                    inner_ds[index], tube_pitch, baffle_pitch, baffle_opening_ratio,
                                    strips_to_constrictions_ratio, tube_to_baffle_clearance, (ro_min + ro_max) / 2,
                                    length_total)

        hi, i_vel, hs_values, tube_count, shell_inner_diameter, ohtc, _, re_in, _bsc, _bd, _epl = iter_results[0]
        nozzle_velocity = 0.1
        dpi = get_dp_i(i_vel, nozzle_velocity, (length_total - 2 * ts_thickness[index]), inner_ds[index],
                       shell_passes * 2, re_in, hi, ohtc)

        min_nozzle_d = 2000 * get_minimum_radius(calculate_water_flow_rate(outlet_temp_wa) / (3600 * rho_wa), rho_wa, 0,
                                                 1.4)
        nozzle_sze, nozzle_sch, nozzle_d, nozzle_th = get_shell_th(min_nozzle_d, nozzle=True,
                                                                   design_pressure=(1.1 * 6 * pow(10, 5)),
                                                                   max_tension=md.ss304_design_stress)
        dps = 2 * get_dp_s(hs_values[7], ohtc, hi, hs_values[8], shell_inner_diameter, _bd,
                           baffle_pitch, hs_values[9], hs_values[11], tube_to_baffle_clearance,
                           _bsc, tube_count * 2, hs_values[12], baffle_opening_ratio,
                           tube_pitch, tube_outer_diameter, hs_values[10], (int((_epl * 1000) / baffle_pitch)) - 1,
                           nozzle_d / 2000, bs_gap)
        if iter_results[1]:
            ohtcs.append(ohtc)
            bps.append(tube_outer_diameter)
            nts.append(tube_count)
            dpts.append(dpi / 100)
            dpss.append(dps / 100)
            sdis.append(shell_inner_diameter)

        with open(f"./Output/trial_{session_id}.txt", "a") as file:
            out, converged = iter_results
            hs_val = out[2]
            out.pop(2)
            file.write(f"{tube_outer_diameter}\n")
            file.write(f"Converged:{converged}\n")
            for tag, value in zip(out_tags, out):
                file.write(f"{tag}:{value}\n")
            file.write(f"Tube Pitch:{tube_pitch}\n")
            for tag, value in zip(hs_tags, hs_val):
                file.write(f"{tag}:{value}\n")
            file.write("\n")

        print(tube_outer_diameter)
        print(iter_results)
        print()


def strips_run(length_total):
    minimum_clearance = 6.4
    bs_gap = 37

    # Baffle
    baffle_pitch = (length_total * 1000) / 20
    baffle_opening_ratio = 0.25
    tube_to_baffle_clearance = 0.8

    shell_passes = 2

    with open(f"./Output/strips_{session_id}.txt", "a") as file:
        file.write(f"Minimum Clearance:{minimum_clearance}\n")
        file.write(f"Bundle-Shell Gap:{bs_gap}\n")
        file.write(f"Baffle Pitch:{baffle_pitch}\n")
        file.write(f"Baffle Opening Ratio:{baffle_opening_ratio}\n")
        file.write(f"Tube to Baffle Clearance:{tube_to_baffle_clearance}\n")
        file.write(f"Shell Passes:{shell_passes}\n")
        file.write(f"Shell Fouling Resistance:{(ro_min + ro_max) / 2}\n")
        file.write("\n")

    for tube_outer_diameter in [16, 20, 25]:
        for strips_to_constrictions_ratio in np.linspace(0, 1, 11):
            pitch_multiplier = 1.5
            tube_pitch = tube_outer_diameter * pitch_multiplier
            if (tube_pitch - tube_outer_diameter) < minimum_clearance:
                tube_pitch = tube_outer_diameter + minimum_clearance

            index = outer_ds.index(tube_outer_diameter)
            iter_results = htc_iterator(500, shell_passes, ts_thickness[index], tube_outer_diameter, bs_gap,
                                        inner_ds[index], tube_pitch, baffle_pitch, baffle_opening_ratio,
                                        strips_to_constrictions_ratio, tube_to_baffle_clearance, (ro_min + ro_max) / 2,
                                        length_total)

            with open(f"./Output/strips_{session_id}.txt", "a") as file:
                out, converged = iter_results

                file.write(f"{tube_outer_diameter}, {strips_to_constrictions_ratio}\n")
                if len(out) == (len(out_tags) + 1):
                    hs_val = out[2]
                    out.pop(2)
                    file.write(f"Converged:{converged}\n")
                    for tag, value in zip(out_tags, out):
                        file.write(f"{tag}:{value}\n")
                    file.write(f"Tube Pitch:{tube_pitch}\n")
                    if len(hs_val) == len(hs_tags):
                        for tag, value in zip(hs_tags, hs_val):
                            file.write(f"{tag}:{value}\n")
                    else:
                        file.write(f"{hs_val}\n")
                else:
                    file.write(f"{out}\n")
                file.write("\n")

            print(tube_outer_diameter)
            print(iter_results)
            print()


def baffle_pitch_run(length_total):
    minimum_clearance = 6.4
    bs_gap = 58
    tube_outer_diameter = 20

    # Baffle
    baffle_opening_ratio = 0.25
    strips_to_constrictions_ratio = 0.3
    tube_to_baffle_clearance = 0.8

    shell_passes = 2

    with open(f"./Output/baffle_pitch_{session_id}.txt", "a") as file:
        file.write(f"Minimum Clearance:{minimum_clearance}\n")
        file.write(f"Bundle-Shell Gap:{bs_gap}\n")
        file.write(f"Tube Outer Diameter:{tube_outer_diameter}\n")
        file.write(f"Baffle Opening Ratio:{baffle_opening_ratio}\n")
        file.write(f"Strips to Constrictions Ratio:{strips_to_constrictions_ratio}\n")
        file.write(f"Tube to Baffle Clearance:{tube_to_baffle_clearance}\n")
        file.write(f"Shell Passes:{shell_passes}\n")
        file.write(f"Shell Fouling Resistance:{(ro_min + ro_max) / 2}\n")
        file.write("\n")

    for baffle_pitch in np.linspace(350, np.interp(20, lm_od, lm_values), 10):
        baffle_pitch = int(baffle_pitch)
        pitch_multiplier = 1.25
        tube_pitch = tube_outer_diameter * pitch_multiplier
        if (tube_pitch - tube_outer_diameter) < minimum_clearance:
            tube_pitch = tube_outer_diameter + minimum_clearance

        index = outer_ds.index(tube_outer_diameter)
        iter_results = htc_iterator(500, shell_passes, ts_thickness[index], tube_outer_diameter, bs_gap,
                                    inner_ds[index], tube_pitch, baffle_pitch, baffle_opening_ratio,
                                    strips_to_constrictions_ratio, tube_to_baffle_clearance, (ro_min + ro_max) / 2,
                                    length_total)

        hi, i_vel, hs_values, tube_count, shell_inner_diameter, ohtc, _, re_in, _bsc, _bd, _epl, _ = iter_results[0]
        nozzle_velocity = 0.1
        dpi = get_dp_i(i_vel, nozzle_velocity, (length_total - 2 * ts_thickness[index]), inner_ds[index],
                       shell_passes * 2, re_in, hi, ohtc)
        dps = 2 * get_dp_s(hs_values[7], ohtc, hi, hs_values[8], shell_inner_diameter, _bd,
                           baffle_pitch, hs_values[9], hs_values[11], tube_to_baffle_clearance,
                           _bsc, tube_count * 2, hs_values[12], baffle_opening_ratio,
                           tube_pitch, tube_outer_diameter, hs_values[10], (int((_epl * 1000) / baffle_pitch)) - 1,
                           123123, bs_gap)
        if iter_results[1]:
            ohtcs.append(ohtc)
            bps.append(baffle_pitch)
            nts.append(tube_count)
            dpts.append(dpi / 100)
            dpss.append(dps / 100)
            sdis.append(shell_inner_diameter)

        with open(f"./Output/baffle_pitch_{session_id}.txt", "a") as file:
            out, converged = iter_results

            file.write(f"{baffle_pitch}\n")
            if len(out) == (len(out_tags) + 1):
                hs_val = out[2]
                out.pop(2)
                file.write(f"Converged:{converged}\n")
                for tag, value in zip(out_tags, out):
                    file.write(f"{tag}:{value}\n")
                file.write(f"Tube Pitch:{tube_pitch}\n")
                file.write(f"Inside Pressure Drop:{dpi}\n")
                file.write(f"Shell Side Pressure Drop:{dps}\n")
                if len(hs_val) == len(hs_tags):
                    for tag, value in zip(hs_tags, hs_val):
                        file.write(f"{tag}:{value}\n")
                else:
                    file.write(f"{hs_val}\n")
            else:
                file.write(f"{out}\n")
            file.write("\n")


def price_comparison_run(tube_outer_diameter, dpi_max, dps_max, length_total, suffix="", sheet=True):
    bs_gap = 38
    tube_to_baffle_clearance = 0.8

    shell_passes = 2
    placeholder_best = 10 ** 16
    previous_best = placeholder_best
    best_value = None
    index = outer_ds.index(tube_outer_diameter)
    _epl = (length_total - 2 * ts_thickness[index])
    previous_shell_diameter = 50
    min_baffles = (int((_epl * 1000) / np.interp(20, lm_od, lm_values))) - 1
    max_baffles = (int((_epl * 1000) / previous_shell_diameter)) - 1
    not_converged = 0
    too_many_baffles = 0
    baffle_ratio_overlap = 0
    bad_ratio = 0
    inner_dp = 0
    outer_dp = 0
    outer_dp_not_converged = 0
    total_cases = 0
    with open(f"./Output/pc{suffix}_{session_id}.txt", "a") as file:
        file.write(str([tube_outer_diameter, dpi_max, dps_max, length_total]) + "\n")
        for baffle_count in range(min_baffles, max_baffles + 1, 1):
            for strips_to_constrictions_ratio in np.linspace(0, 0.5, 11):
                for baffle_opening_ratio in np.linspace(0.2, 0.3, 11):
                    total_cases += 1
                    baffle_count_flag = False
                    pitch_multiplier = 1.25
                    tube_pitch = tube_outer_diameter * pitch_multiplier
                    iter_results = htc_iterator(500, shell_passes, ts_thickness[index], tube_outer_diameter, bs_gap,
                                                inner_ds[index], tube_pitch, baffle_count, baffle_opening_ratio,
                                                strips_to_constrictions_ratio, tube_to_baffle_clearance, ro_selected,
                                                length_total, use_sheet=sheet)
                    if not iter_results[1]:
                        not_converged += 1
                        continue

                    hi, i_vel, hs_values, tube_count, shell_inner_diameter, ohtc, _, re_in, _bsc, _bd, _, tsth, hcwth =\
                        iter_results[0]
                    _epl = (length_total - 2 * tsth)
                    baffle_pitch = int(_epl * 1000) / (baffle_count + 1)
                    if baffle_pitch < (shell_inner_diameter / 5) or baffle_pitch < 51:
                        too_many_baffles += 1
                        continue
                    min_t_nozzle_d = 2000 * get_minimum_radius(85000 / 3600, rho_ipa, 0, 1.5)
                    t_nozzle_sze, t_nozzle_sch, t_nozzle_d, t_nozzle_th = get_shell_th(min_t_nozzle_d, nozzle=True,
                                                                                       design_pressure=(1.1 * 18),
                                                                                       max_tension=md.a20_design_stress)
                    nozzle_velocity = (85000 / (3600 * rho_ipa)) / (np.pi * pow(t_nozzle_d / 2000, 2))
                    dpi = get_dp_i(i_vel, nozzle_velocity, (length_total - 2 * tsth), inner_ds[index],
                                   shell_passes * 2, re_in, hi, ohtc)
                    if dpi < 0 or dpi > dpi_max:
                        inner_dp += 1
                        continue
                    try:
                        min_nozzle_d = 2000 * get_minimum_radius(calculate_water_flow_rate(outlet_temp_wa) / 3600,
                                                                 rho_wa, 0,
                                                                 1.4)
                        nozzle_sze, nozzle_sch, nozzle_d, nozzle_th = get_shell_th(min_nozzle_d, nozzle=True,
                                                                                   design_pressure=(1.1 * 6),
                                                                                   max_tension=md.ss304_design_stress)
                        dps = 2 * get_dp_s(hs_values[7], ohtc, hi, hs_values[8], shell_inner_diameter, _bd,
                                           baffle_pitch, hs_values[9], hs_values[11], tube_to_baffle_clearance,
                                           _bsc, tube_count * 2, hs_values[12], baffle_opening_ratio,
                                           tube_pitch, tube_outer_diameter, hs_values[10], baffle_count,
                                           nozzle_d / 2000, bs_gap)
                    except Exception:
                        outer_dp_not_converged += 1
                        continue

                    # Adjusting for fouling
                    sb_bs_r = shell_inner_diameter / baffle_pitch
                    if sb_bs_r < 1 or sb_bs_r > 5:
                        bad_ratio += 1
                        if baffle_count_flag:
                            baffle_ratio_overlap += 1
                        continue
                    if hs_values[8] > 1000:
                        dps = dps * np.interp(sb_bs_r, sd_bs_ratio, fouled_ratio_t)
                    else:
                        dps = dps * np.interp(sb_bs_r, sd_bs_ratio, fouled_ratio_l)

                    if dps < 0 or dps > dps_max:
                        outer_dp += 1
                        continue

                    shell_size, shell_gauge, shell_inner_diameter, shell_th = get_shell_th(shell_inner_diameter,
                                                                                           use_sheet=sheet)
                    tube_id = dimensions[tube_outer_diameter]
                    baffle_thickness = get_baffle_th(shell_inner_diameter, baffle_pitch)
                    io_depth = (pow(tube_outer_diameter - (2 * tube_id), 2) * tube_count) / shell_inner_diameter
                    # Change 30 to something else if changed
                    # Head thickness is 6.5 including corrosion allowance. change if changed
                    price = get_price(shell_inner_diameter, shell_th,
                                      _epl, shell_passes, shell_inner_diameter, shell_th + corrosion_allowance,
                                      30 + corrosion_allowance, tube_outer_diameter,
                                      tube_id, length_total, 2 * shell_passes, tube_count,
                                      baffle_thickness, baffle_opening_ratio, baffle_count, tsth,
                                      6.5, io_depth,
                                      shell_th + corrosion_allowance, hs_values[-4], _epl, tube_pitch / 1000)
                    file.write(str(price) + ", ")
                    file.write(str([baffle_count, strips_to_constrictions_ratio, baffle_opening_ratio,
                                    tube_outer_diameter, length_total, dpi, dps, shell_size, shell_gauge]) + "\n")
                    if price < previous_best:
                        previous_best = price
                        best_value = [baffle_count, strips_to_constrictions_ratio, baffle_opening_ratio,
                                      tube_outer_diameter, length_total, dpi, dps, shell_size, shell_gauge,
                                      baffle_pitch, baffle_thickness, shell_th,
                                      iter_results[0]]
        file.write("convergence: " + str([not_converged, too_many_baffles, bad_ratio, baffle_ratio_overlap,
                                          inner_dp, outer_dp, outer_dp_not_converged, total_cases]))
        file.write("\n")
    if placeholder_best != previous_best:
        return previous_best, best_value
    return None, None


def finalise_design(tube_outer_diameter=25, dpi_max=42500, dps_max=ceil(42500), length_total=7.32):
    bs_gap = 38
    tube_to_baffle_clearance = 0.8

    shell_passes = 2

    index = outer_ds.index(tube_outer_diameter)
    _epl = (length_total - 2 * ts_thickness[index])

    baffle_count = 12
    strips_to_constrictions_ratio = 0.5
    baffle_opening_ratio = 0.22

    pitch_multiplier = 1.25
    tube_pitch = tube_outer_diameter * pitch_multiplier
    iter_results = htc_iterator(500, shell_passes, ts_thickness[index], tube_outer_diameter, bs_gap,
                                inner_ds[index], tube_pitch, baffle_count, baffle_opening_ratio,
                                strips_to_constrictions_ratio, tube_to_baffle_clearance, ro_selected,
                                length_total, use_sheet=True)
    if not iter_results[1]:
        raise Exception("Not converged")

    hi, i_vel, hs_values, tube_count, shell_inner_diameter, ohtc, _, re_in, _bsc, _bd, _, tsth, hcwth = \
        iter_results[
            0]
    _epl = (length_total - 2 * tsth)
    baffle_pitch = 555
    if baffle_pitch < (shell_inner_diameter / 5) or baffle_pitch < 51:
        raise Exception("Too many baffles")
    min_t_nozzle_d = 2000 * get_minimum_radius(85000 / 3600, rho_ipa, 0, 1.5)
    t_nozzle_sze, t_nozzle_sch, t_nozzle_d, t_nozzle_th = get_shell_th(min_t_nozzle_d, nozzle=True,
                                                                       design_pressure=(1.1 * 18),
                                                                       max_tension=md.a20_design_stress)
    nozzle_velocity = (85000 / (3600*rho_ipa)) / (np.pi * pow(t_nozzle_d/2000, 2))
    print("T_NOZZLE_VELOCITY", nozzle_velocity, min_t_nozzle_d)
    dpi = get_dp_i(i_vel, nozzle_velocity, (length_total - 2 * tsth), inner_ds[index],
                   shell_passes * 2, re_in, hi, ohtc)
    if dpi < 0 or dpi > dpi_max:
        raise Exception("DPI issue")

    min_nozzle_d = 2000 * get_minimum_radius(calculate_water_flow_rate(outlet_temp_wa) / 3600, rho_wa, 0,
                                             1.4)
    nozzle_sze, nozzle_sch, nozzle_d, nozzle_th = get_shell_th(min_nozzle_d, nozzle=True,
                                                               design_pressure=(1.1 * 6),
                                                               max_tension=md.ss304_design_stress)
    print("NOZZLE", min_nozzle_d, nozzle_sze, nozzle_sch, nozzle_d, nozzle_th)
    dps = 2 * get_dp_s(hs_values[7], ohtc, hi, hs_values[8], shell_inner_diameter, _bd,
                       baffle_pitch, hs_values[9], hs_values[11], tube_to_baffle_clearance,
                       _bsc, tube_count * 2, hs_values[12], baffle_opening_ratio,
                       tube_pitch, tube_outer_diameter, hs_values[10], baffle_count, nozzle_d / 2000, bs_gap)
    print("DPS", dps)
    # Adjusting for fouling
    sb_bs_r = shell_inner_diameter / baffle_pitch
    if sb_bs_r < 1 or sb_bs_r > 5:
        raise Exception("Bad ratio")
    if hs_values[8] > 1000:
        dps = dps * np.interp(sb_bs_r, sd_bs_ratio, fouled_ratio_t)
    else:
        dps = dps * np.interp(sb_bs_r, sd_bs_ratio, fouled_ratio_l)

    if dps < 0 or dps > dps_max:
        print("DPS ISSUE", dps)
        raise Exception("DPS issue")

    shell_size, shell_gauge, shell_inner_diameter, shell_th = get_shell_th(shell_inner_diameter, use_sheet=True)
    tube_id = dimensions[tube_outer_diameter]
    baffle_thickness = get_baffle_th(shell_inner_diameter, baffle_pitch)
    io_depth = (pow(tube_outer_diameter - (2 * tube_id), 2) * tube_count) / shell_inner_diameter
    # Change 30 to something else if changed
    # Head thickness is 6.5 including corrosion allowance. change if changed
    price = get_price(shell_inner_diameter, shell_th,
                      _epl, shell_passes, shell_inner_diameter, shell_th + corrosion_allowance,
                      30 + corrosion_allowance, tube_outer_diameter,
                      tube_id, length_total, 2 * shell_passes, tube_count,
                      baffle_thickness, baffle_opening_ratio, baffle_count, tsth,
                      6.5, io_depth,
                      shell_th + corrosion_allowance, hs_values[-4], _epl, tube_pitch / 1000)

    return price, [baffle_count, strips_to_constrictions_ratio, baffle_opening_ratio,
                   tube_outer_diameter, length_total, dpi, dps, shell_size, shell_gauge,
                   baffle_pitch, baffle_thickness, shell_th,
                   iter_results[0], nozzle_sze, nozzle_sch, nozzle_d, nozzle_th,
                   t_nozzle_sze, t_nozzle_sch, t_nozzle_d, t_nozzle_th]


print("DESIGN_OUTPUT", finalise_design())
# trial_run()  # 16, 20, 25 seem good
# strips_run()  # 0.3 seems good (for 20)

# vls = []
# bvs = []
# x_axis = []
# # lengths = [1.83, 2.44, 3.66, 4.88, 6.10, 7.32]
# lengths = [4.88, 6.10, 7.32]  # It was found that the first 3 never converge and were removed for performance reasons.
# lowest = 10**16
# lowest_bv = None
# for length in lengths:
#     xavals = []
#     vlvals = []
#     bvss = []
#     for tod in outer_ds:
#         xavals.append(tod)
#         vl, bv = price_comparison_run(tod, 40000, 40000/1.1, length, sheet=True)
#         if vl is not None and vl < lowest:
#             lowest = vl
#             lowest_bv = bv
#         vlvals.append(vl)
#         bvss.append(bv)
#     x_axis.append(xavals)
#     vls.append(vlvals)
#     bvs.append(bvss)
# print(lowest, lowest_bv)
# print()
# print("\n\n".join(["" if bvsval[2] is None else f"{vls[_ind][2]} {bvsval[2]}" for _ind, bvsval in enumerate(bvs)]))
# print("\n\n")
# # print([None if bvsval[1] is None else bvsval[1][9][6] for bvsval in bvs])
# plt.scatter(x_axis[0], vls[0], label=(str(lengths[0]) + " (m)"), color="#56B4E9", marker="o")
# plt.scatter(x_axis[1], vls[1], label=(str(lengths[1]) + " (m)"), color="#D55E00", marker="X")
# plt.scatter(x_axis[2], vls[2], label=(str(lengths[2]) + " (m)"), color="#F0E442", marker="s")
# plt.legend()
# plt.xlabel("Tube Outer Diameter (mm)", fontname="Times New Roman")
# plt.ylabel("Approximate Cost of Raw Materials (£1000)", fontname="Times New Roman")
# maxvl = max([mvl if mvl is not None else 190000 for mvl in vls[2]])
# maxvl = ceil(maxvl / 10000)
# if maxvl % 2 == 0:
#     maxvl += 1
# maxvl = maxvl * 10000
# plt.ylim([50000, maxvl])
# plt.yticks(range(50000, maxvl + 1, 20000), range(50, (maxvl // 1000)+1, 20), fontname="Times New Roman")
# plt.show()

# vls = [[], []]
# x_axis = list(range(35, 46))
# lengths = [4.88, 6.10, 7.32]
# tods = [20, 25]
# ipa_pumping_cost = get_pumping_cost(18, 85000, rho_ipa)
# for optmp in x_axis:
#     outlet_temp_wa = optmp
#     bulk_temp_wa, mean_cp_wa, rho_wa, mu_wa, mu_wa_bulk, k_wa = recalculate_water_data()
#     # print(bulk_temp_wa, mean_cp_wa, rho_wa, mu_wa, mu_wa_bulk, k_wa, outlet_temp_wa)
#     # print(calculate_water_flow_rate(outlet_temp_wa))
#     vl1s = []
#     vl2s = []
#     print("P", get_pumping_cost(6, calculate_water_flow_rate(outlet_temp_wa), rho_wa) + ipa_pumping_cost)
#     pumping_cost = 0
#     for length in lengths:
#         vl1, bv1 = price_comparison_run(20, 42500, 42500, length, suffix="_temperature")
#         vl2, bv2 = price_comparison_run(25, 42500, 42500, length, suffix="_temperature")
#         if vl1 is not None:
#             vl1s.append(vl1 + pumping_cost)
#         if vl2 is not None:
#             vl2s.append(vl2 + pumping_cost)
#     # enablePrint()
#     if vl1s:
#         print("1", min(vl1s))
#         vls[0].append(min(vl1s))
#     if vl2s:
#         print("2", min(vl2s))
#         vls[1].append(min(vl2s))
#     print(optmp)
#     # blockPrint()
# # plt.scatter(x_axis[3], vls[3], label=(str(lengths[3]) + " (m)"), color="#332288", marker="o")
# # plt.scatter(x_axis[4], vls[4], label=(str(lengths[4]) + " (m)"), color="#44AA99", marker="X")
# # plt.scatter(x_axis[5], vls[5], label=(str(lengths[5]) + " (m)"), color="#DDCC77", marker="s")
# # Third colour = #F0E442
# # enablePrint()
# plt.plot(x_axis, vls[0], marker="o", label="20 (mm)", clip_on=False, color="#56B4E9")
# plt.plot(x_axis, vls[1], marker="X", label="25 (mm)", clip_on=False, color="#D55E00")
# plt.legend()
# minv0 = min(vls[0])
# minv1 = min(vls[1])
# minv = int(min([minv0, minv1]) / 10000)
# maxv0 = max(vls[0])
# maxv1 = max(vls[1])
# maxv = ceil(max([maxv0, maxv1]) / 10000)
# if (minv % 2 == 0) and (maxv % 2 == 0):
#     minv = minv * 10000
#     maxv = maxv * 10000
# elif (minv % 2 == 1) and (maxv % 2 == 1):
#     minv = minv * 10000
#     maxv = maxv * 10000
# else:
#     minv = minv * 10000
#     maxv = (maxv + 1) * 10000
#
# plt.xlabel("Water Outlet Temperature ($^\circ$C)")
# plt.ylabel("Approximate Cost of Raw Materials (£1000)")
# # plt.xlim([35, 45])
# plt.xticks(range(35, 46, 2), range(35, 46, 2))
# plt.ylim([minv, maxv])
# plt.yticks(range(minv, maxv+1, 10000), range(minv//1000, 1+(maxv//1000), 10))
# plt.show()


# vls = []
# bvs = []
# x_axis = []
#
# vls1 = []
# bvs1 = []
# x_axis1 = []
#
# lengths = [4.88, 6.10, 7.32]
# for tod in outer_ds:
#     tod_values = []
#     vl_values = []
#     bv_values = []
#     tod1_values = []
#     vl1_values = []
#     bv1_values = []
#     for length in lengths:
#         vl, bv = price_comparison_run(tod, 40000, 40000 / 1.1, length, sheet=False, suffix="_shell_c")
#         vl1, bv1 = price_comparison_run(tod, 40000, 40000 / 1.1, length, sheet=True, suffix="_shell_c")
#         if vl is not None:
#             vl_values.append(vl)
#             bv_values.append(bv)
#             tod_values.append(tod)
#         if vl1 is not None:
#             vl1_values.append(vl1)
#             bv1_values.append(bv1)
#             tod1_values.append(tod)
#     if vl_values:
#         best_index = vl_values.index(min(vl_values))
#         x_axis.append(tod_values[best_index])
#         vls.append(vl_values[best_index])
#         bvs.append(bv_values[best_index])
#         print("PIPE")
#         print(tod_values[best_index])
#         print(vl_values[best_index])
#         print(bv_values[best_index])
#         print()
#     if vl1_values:
#         best_index1 = vl1_values.index(min(vl1_values))
#         x_axis1.append(tod1_values[best_index1])
#         vls1.append(vl1_values[best_index1])
#         bvs1.append(bv1_values[best_index1])
#         print("SHEET")
#         print(tod1_values[best_index1])
#         print(vl1_values[best_index1])
#         print(bv1_values[best_index1])
#         print()
#
# plt.scatter(x_axis, vls, label="Pipe", color="#332288", marker="o")
# plt.scatter(x_axis1, vls1, label="Rolled Sheet", color="#44AA99", marker="X")
# plt.legend()
# plt.xlabel("Tube Outer Diameter (mm)")
# plt.ylabel("Approximate Cost of Raw Materials (£1000)")
#
# maxvl = max([mvl if mvl is not None else 190000 for mvl in vls])
# maxvl = ceil(maxvl / 10000)
# if maxvl % 2 == 0:
#     maxvl += 1
# maxvl = maxvl * 10000
#
# plt.ylim([50000, maxvl])
# plt.yticks(range(50000, maxvl+1, 20000), range(50, (maxvl//1000)+1, 20))
# plt.show()
