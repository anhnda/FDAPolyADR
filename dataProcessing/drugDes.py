import params
import numpy as np
from utils import utils
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

MORGIN_SIZE = 2048

def genMorganBitVecFromSmiles(s, radius=3):
    m = Chem.MolFromSmiles(s)
    info = {}
    fp = AllChem.GetMorganFingerprintAsBitVect(m, radius, bitInfo=info)
    v = np.zeros(MORGIN_SIZE, dtype=int)
    for i in info.keys():
        v[i] = 1
    return v

def exportMorginFingerprint():
    fin = open("%s/DrugBank/DrugBankNames.txt" % params.DATA_DIR)
    fsmileMissng = open("%s/DrugBank/MissingSMILEs.txt" % params.DATA_DIR, "w")
    dName2Morgan = dict()
    while True:
        line = fin.readline()
        if line == "":
            break
        parts = line.strip().split("||")
        drugName = parts[0].lower()
        # print(drugName)
        t = parts[1]
        smile = parts[-2]
        if t == "small molecule" and len(smile) > 2:
            try:

                v = genMorganBitVecFromSmiles(smile)
                # print(smile)
                dName2Morgan[drugName] = v
            except:
                fsmileMissng.write("%s\n" % drugName)
                continue
    fin.close()
    utils.save_obj(dName2Morgan, "%s/DrugBank/DrugMorganDes" % params.DATA_DIR)

def fillMissingSMILEs():
    fin = open("%s/DrugBank/MissingSMILEsF.txt" % params.DATA_DIR)
    lines = fin.readlines()
    d = utils.load_obj("%s/DrugBank/DrugMorganDes" % params.DATA_DIR)
    for line in lines:
        line = line.strip()
        parts = line.split("||")
        try:
            v = genMorganBitVecFromSmiles(parts[1])
        except:
            print(parts[1])
        d[parts[0].lower()] = v

    utils.save_obj(d, "%s/DrugBank/DrugMorganDes" % params.DATA_DIR)

if __name__ == "__main__":
    exportMorginFingerprint()
    # fillMissingSMILEs()
    pass