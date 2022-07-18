import csv

def get_length(vars):
    ret = 16569
    for var in vars:
        if "." in var:
            ret += int(var.split(".")[1][0])
        if "-" in var and var.endswith("D"):
            del_size = int(var.split("-")[1].split("D")[0]) - int(var.split("-")[0])
            ret -= del_size
    return str(ret)

def write_hsd(x, vars):
    writer = csv.writer(open("../data/hsds/"+str(x)+".hsd", "w"), delimiter="\t")
    row = [x]
    row.append("1-"+get_length(vars))
    for var in vars:
        if var[-1] == "D" or "." in var:
            row.append(var)
        else:
            row.append(var[1:])
    writer.writerow(row)

var_dict = {}
with open("../data/variants.txt", "rt") as f:
    lines = f.readlines()
    for line in lines:
        key = line.split("\t")[0]
        if "[" in key or "]" in key:
            continue
        if "/" in key:
            key = key.split("/")[0]
        if key == '"""M4""""67"""':
            key = "M4''67"
        var_list = [x.replace("\n", "").replace("t", "T").replace("c", "C").replace("a", "A").replace("g", "G").replace("(", "").replace(")", "").replace(",", "").replace("!", "").replace("[", "").replace("]", "") for x in line.split("\t")[1:]]
        var_dict[key] = var_list

    for x in var_dict.keys():
        write_hsd(x, var_dict[x])
