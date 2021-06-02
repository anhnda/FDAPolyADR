import params
from utils import utils
import numpy as np
import codecs
import csv
import io

CAD_FOLDER_INP = "/home/ubuntu/Downloads/extract_extrait"

invalidSes = ["drug", "death", "product", "off label use"]
def exportReactionsFile():
    fin = codecs.open("%s/reactions.txt" % CAD_FOLDER_INP)
    fout = open("%s/Reactions1.txt" % params.CAD_OUT, "w")

    dId2Ses = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        ios = io.StringIO(line.strip().lower())
        vv = list(csv.reader(ios, delimiter='$'))[0]
        # print( vv)
        sId = vv[1]
        seName = vv[5]
        isInValid = False
        for invalidSe in invalidSes:
            if seName.__contains__(invalidSe):
                isInValid = True
                break
        if isInValid:
            continue

        seList = utils.get_insert_key_dict(dId2Ses, sId, set())
        seList.add(seName)

        # print(cId,  currentDrugs)
    fin.close()

    for k,v in dId2Ses.items():
        fout.write("%s$%s\n" % (k, ",".join(list(v))))
    fout.close()


def exportIndicationFile():
    fin = codecs.open("%s/report_drug_indication.txt" % CAD_FOLDER_INP)
    fout = open("%s/Indications1.txt" % params.CAD_OUT, "w")

    dId2Ses = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        ios = io.StringIO(line.strip().lower())
        vv = list(csv.reader(ios, delimiter='$'))[0]
        # print( vv)
        sId = vv[1]
        indcName = vv[4]

        indcList = utils.get_insert_key_dict(dId2Ses, sId, set())
        indcList.add(indcName)

        # print(cId,  currentDrugs)
    fin.close()

    for k,v in dId2Ses.items():
        fout.write("%s$%s\n" % (k, ",".join(list(v))))
    fout.close()
if __name__ == "__main__":
    exportReactionsFile()
    # exportIndicationFile()