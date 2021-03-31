import params
from utils import utils

N_SE = 800
N_SEG = 40
N_FILE = int(N_SE / N_SEG)

def exportPolySEs():
    drugDesMap = utils.load_obj("%s/DrugBank/DrugMorganDes" % params.DATA_DIR)
    seDict = dict()
    dComb2Se = dict()
    fin = open("%s/FTest/FileMap.txt" % params.FADER_OUT)
    hashFiles = fin.readlines()
    ln = min(N_FILE, len(hashFiles))
    hashFiles = hashFiles[:ln]
    for hashId in hashFiles:
        parts = hashId.strip().split("\t")
        hashId = parts[0]
        ses = parts[1].split("__")
        for se in ses:
            utils.get_update_dict_index(seDict, se)
        path = "%s/FTest/%s" % (params.FADER_OUT, hashId)
        print("Reading... ", path)
        polySes = open(path).readlines()
        for polySe in polySes:
            polySe = polySe.strip().split("_")
            drugComb = polySe[0]
            se = polySe[1].split("\t")[0]
            drugs = drugComb.split(",")
            isValidComb = True
            # print(drugs)
            for drug in drugs:
                if drug not in drugDesMap:
                    isValidComb = False
                    break

            if isValidComb:
                # print(drugComb)
                sel = utils.get_insert_key_dict(dComb2Se, drugComb, [])
                sel.append(se)

    fout = open("%s/PolySes.txt" % params.FADER_OUT, "w")
    for drugComb, ses in dComb2Se.items():
        fout.write("%s\t%s\n" % (drugComb, ",".join(ses)) )
    fout.close()

if __name__ == "__main__":
    exportPolySEs()



