import params
from utils import utils
def getFDASEs():
    fin = open("%s/polyDrugADR.txt" % params.DATA_DIR)
    se1 = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("|")
        ses = parts[-1].split(",")
        for se in ses:
            utils.add_dict_counter(se1, se)
    fin.close()
    fin = open("%s/CADER.txt" % params.CAD_OUT)
    se2 = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("$")
        ses = parts[-1].split(",")
        for se in ses:
            utils.add_dict_counter(se2, se)

    kvs1 = utils.sort_dict(se1)
    kvs2 = utils.sort_dict(se2)
    print(len(kvs1), len(kvs2))
    k1 = set()
    k2 = set()
    MIN_T = 130
    for kv in kvs1:
        k, v = kv
        if v >= MIN_T:
            k1.add(k)

    for kv in kvs2:
        k, v = kv
        if v >= MIN_T:
            k2.add(k)

    n1 = 0
    n2 = 0
    for k in k1:
        if k not in k2:
            n1 += 1

    for k in k2:
        if k not in k1:
            n2 += 1
    print(len(k1), len(k2), n1, n2, n1 / len(k1), n2 / len(k2),)

def getFDADrug():
    fin = open("%s/polyDrugADR.txt" % params.DATA_DIR)
    se1 = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("|")
        ses = parts[0].split(",")
        for se in ses:
            utils.add_dict_counter(se1, se)
    fin.close()
    fin = open("%s/CADER.txt" % params.CAD_OUT)
    se2 = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("$")
        ses = parts[1].split(",")
        for se in ses:
            utils.add_dict_counter(se2, se)

    kvs1 = utils.sort_dict(se1)
    kvs2 = utils.sort_dict(se2)
    print(len(kvs1), len(kvs2))
    k1 = set()
    k2 = set()
    MIN_T = 5
    for kv in kvs1:
        k, v = kv
        if v >= MIN_T:
            k1.add(k)

    for kv in kvs2:
        k, v = kv
        if v >= 60:
            k2.add(k)

    n1 = 0
    n2 = 0
    for k in k1:
        if k not in k2:
            n1 += 1

    for k in k2:
        if k not in k1:
            n2 += 1
    print(len(k1), len(k2), n1, n2, n1 / len(k1), n2 / len(k2),)
if __name__ == "__main__":
    # getFDASEs()
    getFDADrug()

