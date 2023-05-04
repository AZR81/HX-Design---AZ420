import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rc('font', family='Times New Roman')

def analyse_log():
    values = np.zeros(66)

    with open("./Output/pc_temperature_1681748171.txt", "r") as file:
        lines = file.readlines()
        init = lines[0]
        index = -1
        vls = []
        for line in lines[:-1]:
            if line[0] == "[":
                index += 1
            elif line[0] != "c":
                ln = line.split(", ")
                vls.append(float(ln[0]))
            else:
                if vls:
                    values[index] = min(vls)
                    vls = []
                else:
                    values[index] = None

    # [4.88, 6.10, 7.32]
    temps = list(range(35, 46))
    for i in range(0, 66, 6):
        print(temps[i//6])
        print(values[i], values[i+2], values[i+4])
        print(values[i+1], values[i+3], values[i+5])
        print()


def make_presentable():
    out_tags = ["Tube HTC", "Tube Velocity", "Tubes Per Pass", "Shell Inner Diameter",
                "Overall HTC", "HT Area", "Tube-Side Reynolds Number", "Baffle-Shell Clearance", "Bundle Diameter",
                "Tube Wall Resistance", "Tube Sheet Thickness", "Head Radial Wall Thickness"]
    hs_tags = ["Shell HTC", "Total Shell F", "Ideal Shell HTC",
               "Tube Row CF (Fn)", "Window CF (Fw)", "Bypass CF (Fb)", "Leakage CF (Fl)",
               "Cross-Flow Area (As)", "Reynolds Number", "Number of Sealing Strips", "Shell Velocity",
               "Constriction Rows", "Tubes in Window Zone", "Ideal Baffle Ratio"]
    comparison_tags = ["Baffle Count", "Strips Ratio", "Baffle Ratio", "Tube OD", "Total Length", "DPI", "DPS",
                       "Shell Size", "Shell Schedule", "Baffle Pitch", "Baffle Thickness", "Shell Thickness",
                       "Shell Nozzle Size", "Shell Nozzle Schedule", "Shell Nozzle ID", "Shell Nozzle Thickness",
                       "Tube Nozzle Size", "Tube Nozzle Schedule", "Tube Nozzle ID", "Tube Nozzle Thickness"]

#([]])

    # price = 114889.88566219574
    # comparison_data = [12, 0.31, 0.2, 25, 7.32, 6331.281126869506, 36005.09858664505, None, None]
    # out = [3663.749730975416, 0.7391094169102934, 4956.978827696482, 1.0134221276426791, 1.0993186843574136,
    #        0.9427201528963384, 0.7037395390730485, 0.08282307692307692, 20146.074424007445, 4, 0.6006493761645959, 14,
    #        40]
    # text = [720.2674545978099, 0.5243372897360588, 162, 740, 405.11576526732216, 308.19029208130604, 10887.824453501458,
    #         4.8, 702, 0.0001389569682896662, 0.0225, 5.8]
    # [, []]

    # price = 111484.40969996927
    # comparison_data = [12, 0.5, 0.21000000000000002, 25, 7.32, 6490.607049906753, 35534.814909360626, None, None, 559.6153846153846, 6.4, 4.8]
    # out = [3811.8649398676093, 0.7664200934290387, 4973.597342435207, 1.0093167070967735, 1.087030433757028, 0.9908664762956737, 0.7049893909154283, 0.08237538461538463, 20255.56395892052, 6, 0.603913774948099, 13, 43]
    # text = [728.8595944154226, 0.5308915058577596, 160, 736, 409.99479784075504, 303.9132493386017, 11023.922259170226, 4.8, 698, 0.0001389569682896662, 0.0223, 5.8]

    # price = 111484.40969996927
    # comparison_data = [12, 0.5, 0.21, 25, 7.32, 7354.575106543272, 40921.617876917764, None, None, 559.6153846153846, 6.4, 4.8, '8', '5S', 213.56, 2.77, '8', '5S', 213.56, 2.77]
    # out = [3811.8649398676093, 0.7664200934290387, 4973.597342435207, 1.0093167070967735, 1.087030433757028, 0.9908664762956737, 0.7049893909154283, 0.08237538461538463, 20255.56395892052, 6, 0.603913774948099, 13, 43]
    # text = [728.8595944154226, 0.5308915058577596, 160, 736, 409.99479784075504, 303.9132493386017, 11023.922259170226, 4.8, 698, 0.0001389569682896662, 0.0223, 5.8]

    # price = 106746.65818200306
    # comparison_data = [12, 0.5, 0.21, 25, 7.32, 7826.565337704809, 41995.66681882218, None, None, 559.7692307692307, 6.4, 4.8, '8', '5S', 213.56, 2.77, '8', '5S', 213.56, 2.77]
    # out = [3816.530774253936, 0.7667714696633318, 4977.403209759055, 1.0093167070967735, 1.0874715889660294, 0.9907030154901612, 0.7051428044820884, 0.08094263076923078, 20614.10477841648, 6, 0.6146035656749295, 13, 41]
    # text = [788.769700286765, 0.5551806597205329, 153, 723, 431.1778353766237, 290.02191755128376, 12226.008984477121, 4.8, 685, 0.0001389569682896662, 0.0215, 5.7]

    # # Adjusted
    # price = 107613.71615378193
    # comparison_data = [12, 0.5, 0.24, 25, 7.32, 7661.277884872951, 36933.56857274917, None, None, 559.6923076923077, 6.4, 4.8, '8', '5S', 213.56, 2.77, '8', '5S', 213.56, 2.77]
    # out = [3696.471036196064, 0.745285720931811, 4959.80391463084, 1.0047787100744412, 1.0467569990063255, 1.0, 0.7086087369904408, 0.08137926153846155, 20503.502245827145, 6, 0.6113059831879492, 12, 52, 0.24]
    # text = [779.1696817926456, 0.5480170383047841, 155, 727, 426.3206785861148, 294.2029992804314, 12068.254029838707, 4.8, 689, 0.0001389569682896662, 0.0217, 5.7]

    # # Original 0.22
    # price = 106289.8081825576
    # comparison_data = [12, 0.5, 0.22, 25, 7.32, 7826.413536048627, 40023.73466373172, None, None, 559.7692307692307, 6.4, 4.8, '8', '5S', 213.56, 2.77, '8', '5S', 213.56, 2.77]
    # out = [3801.461374641631, 0.7637715671164118, 4977.2229529228125, 1.0047787100744412, 1.0762011466392951, 1.0, 0.7063169233716858, 0.08094263076923078, 20614.10477841648, 6, 0.6146035656749295, 12, 44, 0.24]
    # text = [788.769700286765, 0.5551806597205329, 153, 723, 430.9848186263497, 290.25695603120016, 12226.008984477121, 4.8, 685, 0.0001389569682896662, 0.0215, 5.7]

    # Correct Orientation
    # price = 106289.8081825576
    # comparison_data = [12, 0.5, 0.22, 25, 7.32, 7826.225708787481, 40399.41194832101, None, None, 559.7692307692307, 6.4, 4.8, '8', '5S', 213.56, 2.77, '8', '5S', 213.56, 2.77]
    # out = [3782.959072825934, 0.7600882335515281, 4976.9999137468285, 1.0093167070967735, 1.0762011466392951, 0.9907030154901612, 0.7063169233716858, 0.08094263076923078, 20614.10477841648, 6, 0.6146035656749295, 13, 44, 0.22]
    # text = [788.769700286765, 0.5551806597205329, 153, 723, 430.7459678103365, 290.25695603120016, 12226.008984477121, 4.8, 685, 0.0001389569682896662, 0.0215, 5.7]

    # Final values
    price = 106289.8081825576
    comparison_data = [12, 0.5, 0.22, 25, 7.32, 7826.413536048627, 40201.44241750104, None, None, 555, 6.4, 4.8, '8', '5S', 213.56, 2.77, '8', '5S', 213.56, 2.77]
    out = [3801.461374641631, 0.7637715671164118, 4977.2229529228125, 1.0047787100744412, 1.0762011466392951, 1.0, 0.7063169233716858, 0.08094263076923078, 20614.10477841648, 6, 0.6146035656749295, 12, 44, 0.24]
    text = [788.769700286765, 0.5551806597205329, 153, 723, 430.9848186263497, 290.25695603120016, 12226.008984477121, 4.8, 685, 0.0001389569682896662, 0.0215, 5.7]

    print("Price:", price)
    print()
    for ct, cd in zip(comparison_tags, comparison_data):
        print(ct, ":", cd)
    print()
    for ot, od in zip(out_tags, text):
        print(ot, ":", od)
    print()
    for ht, hd in zip(hs_tags, out):
        print(ht, ":", hd)


def regenerate_pump_graph():
    x_axis = []
    vls = [[], []]
    with open("final_w_t_out_data_1.txt", "r") as file:
        lines = [line.split(" ") for line in file.readlines()]
        for i in range(0, len(lines), 4):
            pump_cost = float(lines[i][1][:-1])
            x_axis.append(int(lines[i+3][0][:-1]))
            vls[0].append(pump_cost + float(lines[i+1][1][:-1]))
            vls[1].append(pump_cost + float(lines[i+2][1][:-1]))
    plt.plot(x_axis, vls[0], marker="o", label="20 (mm)", clip_on=False, color="#56B4E9")
    plt.plot(x_axis, vls[1], marker="X", label="25 (mm)", clip_on=False, color="#D55E00")
    plt.legend()
    # minv = 80000
    # maxv = 130000
    # plt.ylim([minv, maxv])
    # plt.yticks(range(minv, maxv + 1, 10000), range(minv // 1000, 1 + (maxv // 1000), 10))
    plt.xlabel("Water Outlet Temperature ($^\circ$C)")
    plt.ylabel("Approximate Cost of Raw Materials (£10$^6$)")
    plt.xticks(range(35, 46, 2), range(35, 46, 2))
    minv = 3000000
    maxv = 4200000
    plt.ylim([minv, maxv])
    plt.yticks(range(minv, maxv + 1, 200000), np.round(np.linspace(minv / 1000000, maxv / 1000000, 7), 1))
    plt.show()


def regenerate_pipe_graph():
    x_axis = [16, 20, 25, 30]
    vls = [71692.96760528909, 94738.39370786907, 120535.37564382514, 191498.64417032793]
    vls1 = [67717.55824141034, 89122.87257512954, 111484.40969996927, 183823.30678059708]
    plt.scatter(x_axis, vls, label="Pipe", color="#56B4E9", marker="o")
    plt.scatter(x_axis, vls1, label="Rolled Sheet", color="#D55E00", marker="X")
    plt.legend()
    plt.xlabel("Tube Outer Diameter (mm)")
    plt.ylabel("Approximate Cost of Raw Materials (£1000)")

    maxvl = 200000

    plt.ylim([60000, maxvl])
    plt.yticks(range(60000, maxvl+1, 20000), range(60, (maxvl//1000)+1, 20))
    plt.show()


def compare_radius():
    radius = np.linspace(0.1, 0.2, 200)
    mdot = 178100/3600
    vdot = mdot / 994.02637005113716
    def get_velocity(rad):
        area = np.pi * pow(rad, 2)
        return vdot / area
    velocities = np.apply_along_axis(get_velocity, 0, radius)
    plt.plot(radius, velocities)
    plt.show()


def get_minimum_radius(mass_flow_rate, density, max_head, max_velocity=0.0):
    if max_velocity > 0:
        max_head = density * max_velocity * max_velocity
    volumetric_flow_rate = mass_flow_rate/density
    return pow((density * pow(volumetric_flow_rate, 2))/(max_head*np.pi*np.pi), 1/4)


make_presentable()
