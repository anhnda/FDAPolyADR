from utils import utils
import params

def statsx():
    fin = open("%s/SUB/F2" % params.FADER_OUT)
    nr = 0
    dDrug2NReports = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.split("$")
        dDrug = parts[0].split(",")
        ses = parts[1].split(",")

