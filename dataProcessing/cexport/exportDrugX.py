import params
from utils import utils
import numpy as np
import codecs
import csv
import io

CAD_FOLDER_INP = "/home/ubuntu/Downloads/extract_extrait"


def loadDrugBank():
    fin = open("%s/Data/DrugBank/DrugBankNames.txt" % params.C_DIR)
    dName2Inchi = dict()
    dSyn2Name = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()

        parts = line.split("||")
        if parts[1] != "small molecule":
            continue
        drugName = parts[0]

        inchi = parts[3]
        dName2Inchi[drugName] = inchi
        syns = parts[-1].split("|")
        # if line.__contains__("phisohex"):
        #     print("Here")
        #     print(line)
        #     print(syns)
        #     exit(-1)
        for syn in syns:
            dSyn2Name[syn] = drugName
    fin.close()
    print("Drugbank: ", len(dName2Inchi))
    return dSyn2Name, dName2Inchi

def exportReportDrugFile():
    dSyn2Name, _ = loadDrugBank()
    fin = codecs.open("%s/report_drug.txt" % CAD_FOLDER_INP)
    fout = open("%s/ReportDrug1.txt" % params.CAD_OUT, "w")
    pId = -1
    currentDrugs = set()
    cc = 0
    nMiss = 0
    missingDrugs = set()
    print(dSyn2Name['phisohex'])
    while True:
        line = fin.readline()
        if line == "":
            break
        ios = io.StringIO(line.strip().lower())
        vv = list(csv.reader(ios, delimiter='$'))[0]
        # print( vv)
        cId = vv[1]
        if len(cId) != 9:
            continue
        # print(cId)
        drugName = vv[3]
        cc += 1
        if cc % 100 == 0:
            print("\r%s" % cc, end= "")
        if cId != pId:
            if pId != -1:
                isSelected = True
                drugBankNames = []
                for drug in currentDrugs:
                    drugBankName = utils.get_dict(dSyn2Name, drug, -1)
                    if drugBankName != -1:
                        drugBankNames.append(drugBankName)
                    else:
                        missingDrugs.add(drug)
                        # print("Miss: ", drugName)
                        isSelected = False
                        nMiss += 1

                        #  print("Skip: ", cId, drug)
                        if nMiss % 10000 == 0:
                            # print("Miss: ", nMiss, cc)
                            pass
                        break
                if isSelected:
                    # print("Write file")
                    fout.write("%s$%s\n" % (pId, ",".join(sorted(drugBankNames))))

            pId = cId
            currentDrugs = set()
        currentDrugs.add(drugName)
        # print(cId,  currentDrugs)
    fin.close()
    fout.close()
    print(list(missingDrugs))

def checkDupR():
    fin = open("%s/ReportDrug1.txt" % params.CAD_OUT)
    dCout = dict()
    nError = 0
    cc = 0
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("$")
        idx = parts[0]
        utils.add_dict_counter(dCout, idx)
        cc += 1

    print("Total: ", nError, cc)
    kvs = utils.sort_dict(dCout)
    fout = open("%s/S1.txt" % params.CAD_OUT, "w")
    for kv in kvs:
        k, v = kv
        fout.write("%s\t%s\n" % (k, v))
    fout.close()

def merger():
    fin = open("%s/ReportDrug1.txt" % params.CAD_OUT)
    fout = open("%s/ReportDrug2.txt" % params.CAD_OUT, "w")

    dCout = dict()
    nError = 0
    cc = 0
    dId2Drugs = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        line = line.strip()
        parts = line.split("$")
        id = parts[0]
        drugs = parts[1].split(",")
        drugSet = utils.get_insert_key_dict(dId2Drugs, id, set())
        for drug in drugs:
            drugSet.add(drug)
    fin.close()
    for k,v in dId2Drugs.items():
        fout.write("%s$%s\n" %(k, ",".join(sorted(list(v)))))
    fout.close()


if __name__ == "__main__":
    # exportReportDrugFile()
    # checkDupR()
    merger()
