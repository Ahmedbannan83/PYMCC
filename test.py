import pymcc

tpr = {
    "X06MV": 0.681,
    "X06_FFF": 0.676,
    "X06_FFF_CK": 0.662,
    "X10": 0.738,
    "X10_FFF": 0.725,
    "X18": 0.777,
}
# read mcc file and return list of measurement objects (PDD and/or Profiles)
mymcc = pymcc.readmcc.read_file("PDD.mcc")

mcc_dict = {}
for i in mymcc:
    mcc_dict[i.curve_type] = i.calc_results()

# provide object dict for composite tests
#pdd = mcc_dict
print(mcc_dict)