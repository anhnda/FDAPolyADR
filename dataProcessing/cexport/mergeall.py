import params
from utils import utils

def loadFile(path):
    d = dict()
    fin = open(path)
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("$")
        d[parts[0]] = parts[1]
    fin.close()
    return d
def mergeAll():
    dId2Drugs = loadFile("%s/ReportDrug2.txt" % params.CAD_OUT)
    dId2Indications = loadFile("%s/Indications1.txt" % params.CAD_OUT)
    dId2Ses = loadFile("%s/Reactions1.txt" % params.CAD_OUT)

    fout = open("%s/CADER.txt" % params.CAD_OUT, "w")
    for k, v in dId2Drugs.items():
        ses = utils.get_dict(dId2Ses, k, -1)
        indc = utils.get_dict(dId2Indications, k, "")
        if ses != -1 and len(ses) > 1:
            fout.write("%s$%s$%s$%s\n"  % (k, v, indc, ses))
    fout.close()
if __name__ == "__main__":
    mergeAll()