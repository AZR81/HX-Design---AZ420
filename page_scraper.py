import requests

URL = "https://www.shell-tube.com/americanindustrial/EAB/EAB-701-B4-FP-B-Z.html"
page = requests.get(URL)
with open("model_output.txt", "w") as f:
    with open("model_links.txt", "r") as sr_f:
        sr_fls = sr_f.readlines()

    for line in sr_fls:
        if "href=" in line:
            print(line)
            le = line.split('"')[1]
            print(le)
            URL = f"https://www.shell-tube.com/americanindustrial/EAB/{le}"
            page = requests.get(URL)
            pls = page.text.split("\n")
            for pl in pls:
                if '<td class="products_body_list">' in pl:
                    f.write(pl.split('<td class="products_body_list">')[1].split("<")[0])
                    f.write("\n")
            f.write("\n")
