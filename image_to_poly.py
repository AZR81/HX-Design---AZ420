import cv2
import numpy as np
import matplotlib.pyplot as plt

def evaluate_polynomial(coefficients, value):
    order = len(coefficients) - 1
    output = 0
    for index in range(order + 1):
        output += coefficients[order - index] * pow(value, index)

    return output


def extract_data(name, xr, xo, yr, yo):
    img = cv2.imread(name)
    height, width, _ = img.shape
    fw_rw_pixels = []
    fw_pixels = []
    for x in range(width):
        for y in range(height):
            if img[y, x, 1] < 100:
                fw_rw_pixels.append(x)
                fw_pixels.append(y)
    fw_rw = [xo + (xr * (pixel/width)) for pixel in fw_rw_pixels]
    fw_values = [yo + (yr * ((height - pixel)/height)) for pixel in fw_pixels]
    return fw_rw, fw_values


def get_fw_c():
    rv = 0.8714285714285714
    fw_rw, fw_values = extract_data("fw_values.png", rv, 0, 0.6, 0.6)
    return np.polyfit(fw_rw, fw_values, 20)


def get_fb_c_l():
    fb_ar, fb_values = extract_data("fb_low_values.png", 0.4, 0, 0.5, 0.5)
    return np.polyfit(fb_ar, fb_values, 4)


def get_fb_c_h():
    fb_ar, fb_values = extract_data("fb_high_values.png", 0.4, 0, 0.5, 0.5)
    return np.polyfit(fb_ar, fb_values, 4)


def get_bl_c():
    fb_ar, fb_values = extract_data("bl_values.png", 0.8, 0, 0.5, 0)
    return np.polyfit(fb_ar, fb_values, 12), fb_ar[-1]


def get_r_c():
    fb_ar, fb_values = extract_data("cut_r_values.png", 0.3, 0.15, 0.45, 0)
    return np.polyfit(fb_ar, fb_values, 6), fb_ar, fb_values


def get_theta_c():
    fb_ar, fb_values = extract_data("cut_theta_values.png", 0.3, 0.15, 1.8, 1.4)
    return np.polyfit(fb_ar, fb_values, 6)


def get_jf_data():
    return extract_data("jf_values.png", 5, 1, 3, -3)


def get_p_jf_data():
    return extract_data("pressure_jf_values.png", 5, 1, 3, -2)


def get_p_fb_c_l():
    fb_ar, fb_values = extract_data("pressure_fb_low_values.png", 0.4, 0, 0.8, 0.2)
    return np.polyfit(fb_ar, fb_values, 9), fb_ar[-1]


def get_p_fb_c_h():
    fb_ar, fb_values = extract_data("pressure_fb_high_values.png", 0.4, 0, 0.8, 0.2)
    return np.polyfit(fb_ar, fb_values, 9), fb_ar[-1]


def get_p_bl_c():
    fb_ar, fb_values = extract_data("p_bl_values.png", 0.8, 0, 0.7, 0)
    return fb_ar, fb_values, fb_ar[-1]


def compare_values(xo, yo, c, x0, y0, xn, yn):
    fw_rwval = np.linspace(x0, xn, 100)
    fw_vals = [evaluate_polynomial(c, val) for val in fw_rwval]
    plt.plot(xo, yo)
    plt.plot(fw_rwval, fw_vals)
    plt.xlim([x0, xn])
    plt.ylim([y0, yn])
    plt.show()


def compare_values_interp(xo, yo, x0, y0, xn, yn):
    fw_rwval = np.linspace(x0, xn, 100)
    fw_vals = np.interp(fw_rwval, xo, yo)
    plt.plot(xo, yo)
    plt.plot(fw_rwval, fw_vals)
    plt.xlim([x0, xn])
    plt.ylim([y0, yn])
    plt.show()
