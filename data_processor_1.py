f = open("model_output.txt", "r")
fls = f.readlines()
f.close()

shell = {}
tubes = {}
ends = {}
plates = {}

for fl in fls:
    if "Shell: " in fl:
        material = fl.split("Shell: ")[1][:-1]
        if material not in shell:
            shell.update({material: 1})
        else:
            shell[material] += 1

    if "Tubes: " in fl:
        material = fl.split("Tubes: ")[1][:-1]
        if material not in tubes:
            tubes.update({material: 1})
        else:
            tubes[material] += 1

    if "End Bonnets: " in fl:
        material = fl.split("End Bonnets: ")[1][:-1]
        if material not in ends:
            ends.update({material: 1})
        else:
            ends[material] += 1

    if "Nameplate: " in fl:
        material = fl.split("Nameplate: ")[1][:-1]
        if material not in plates:
            plates.update({material: 1})
        else:
            plates[material] += 1

print(shell)
print(tubes)
print(ends)
print(plates)
