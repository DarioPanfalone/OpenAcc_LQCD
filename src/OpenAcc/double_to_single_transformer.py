# This script does some work for converting 
# 

#!/usr/bin/python

fileNames = [\
'alloc_vars.c',\
'alloc_vars.h',\
'cayley_hamilton.h',\
'fermion_force.c',\
'fermion_force.h',\
'fermion_force_utilities.c',\
'fermion_force_utilities.h',\
'fermionic_utilities.c',\
'fermionic_utilities.h',\
'fermion_matrix.c',\
'fermion_matrix.h',\
'inverter_full.c',\
'inverter_full.h',\
'inverter_multishift_full.c',\
'inverter_multishift_full.h',\
'ipdot_gauge.c',\
'ipdot_gauge.h',\
'matvecmul.h',\
'md_integrator.c',\
'md_integrator.h',\
'plaquettes.c',\
'plaquettes.h',\
'random_assignement.c',\
'random_assignement.h',\
'rettangoli.c',\
'rettangoli.h',\
'single_types.h',\
'stouting.c',\
'stouting.h',\
'struct_c_def.h',\
'su3_measurements.c',\
'su3_measurements.h',\
'su3_utilities.c',\
'su3_utilities.h']


dpFunctionNames = [] # new function names will just be dp function names + '_f' at the end
dpVariableNames = [] # new function names will just be dp function names + '_f' at the end
# as in struct_c_def.h
dpTypes = ['double_soa','dcomplex_soa','vec3_soa','vec3','su3_soa','thmat_soa','tamat_soa']
# corresponding types in in sp_struct_c_def.h
spTypes = ['float_soa','fcomplex_soa','vec3_soa_f','f_vec3','f_su3_soa','f_thmat_soa','f_tamat_soa']





dpToSpDict = dict(zip(dpTypes,spTypes))

# this function is needed to populate the list 'dpFunctionNames'


def findFunctionNames(lineRaw):
    returnTypes=['int','double','void']
    foundFunction = False
    line = lineRaw.strip()
    for returnType in returnTypes:
        initId = line.find(returnType) # looks for '(type)'
        endId = line.find('(')         # looks for '('
        checkEqualId = line.find('=')  # looks for '='
        checkCommentId = line.find('//')  # looks for '//'
        # if the line is in the format '(type) function_name( ****'
        # we found a function name
        foundFunction = initId == 0 and endId != -1 and initId < endId 
        foundFunction = foundFunction and ( checkEqualId > endId or checkEqualId == -1 )
        foundFunction = foundFunction and ( checkCommentId > endId or checkCommentId == -1 )
        if foundFunction:
            initId +=  len(returnType)
            dpFunctionName = line[initId:endId].strip() # removes spaces before and after
            if ' ' in dpFunctionName:  # so that this fails if, e.g.,
                                       # there is 'void static inline' instead of 
                                       # 'static inline void'
                print "ERROR in function name! " , dpFunctionName 
                foundFunction = False
            break
    if foundFunction:
        if dpFunctionName not in dpFunctionNames:
            print "Found function " + dpFunctionName
            dpFunctionNames.append(dpFunctionName)



def findFirstVarName(text):
    foundVariable = False 
    foundSoaType = []
    initId = len(text)
    for soaType in dpTypes:
        initId2 = text.find(soaType)
        con = initId2 > -1 
        con = con and (text[initId2+len(soaType):].strip(' *')[0] not in [')',';'])
        con = con and (text[initId2+len(soaType)+1] in [' ','*'])
        if con:
            foundVariable = True
            if initId2 < initId:
                initId = initId2
                foundSoaType = soaType

    if foundVariable:
        stringsToIgnore = [' ','*','const','__restrict']
        goOn = True 
        newText = str(text[initId+len(foundSoaType):])
        while goOn:
            newText2 = str(newText)
            for stringToIgnore in stringsToIgnore:
                if newText2[:len(stringToIgnore)] == stringToIgnore:
                    newText2 = newText2[len(stringToIgnore):]
            goOn = not (newText == newText2)
            newText = newText2
        foundVariableName = newText.split()[0].strip(',; /')
        print "Found \'" + foundVariableName + "\'" 
        dpVariableNames.append(foundVariableName)
        return text.find(foundVariableName) + len(foundVariableName) 

    else:
        return -1


# collecting all function names
for fileName in fileNames:
    print "----------\nChecking file" , fileName , '\n-----------'
    f  =open(fileName)
    fileLines = f.readlines()
    for line in fileLines:
        findFunctionNames(line)
    f.close()

functionNamesFile = open('functions_found.txt','w')
for foundFunction in dpFunctionNames:
    functionNamesFile.write(foundFunction + '\n')
functionNamesFile.close()
print "Total functions found: ", len(dpFunctionNames)

# collecting all global variable names
print "\n\n-----------------\nLooking for variable names\n-----------------\n"
f = open('alloc_vars.c')
text = f.read()
foundSomething = True
index = 0
while foundSomething:
     res = findFirstVarName(text[index:])
     if res == -1:
         foundSomething = False
     else:
         index += res
f.close()

globalVarNamesFile = open('glvar_found.txt','w')
for foundGlobalVar in dpVariableNames:
    globalVarNamesFile.write(foundGlobalVar + '\n')
globalVarNamesFile.close()

dpVariableNames.append('phases') # for U1 used in dirac matrix

print "Total global variables found: ", len(dpVariableNames)



# collected all double precision function names, now just replacing:
# - all function names
# - all relevant types

for fileName in fileNames:
    f = open(fileName,'r')
    text = f.read()
    newText = str(text)
    # replacing function names
    for dpFunctionName in dpFunctionNames:
        newText = newText.replace(dpFunctionName, dpFunctionName + '_f')
    # changing relevant (soa-like, arrays) types
    # note : this step could also change function names
    for dptype in dpToSpDict:
        newText = newText.replace(dptype,dpToSpDict[dptype])
    # changing filenames in '#includes'
    for fileName2 in fileNames:
        newText = newText.replace(fileName2,'sp_'+fileName2)
    #taking care of header guards
    headerGuard = fileName.upper().replace('.','_')
    if headerGuard not in newText:
        print "Warning, header guard \'" + headerGuard + "\' not found in file " + fileName
    newText = newText.replace(headerGuard, "SP_"+headerGuard)
    # taking care of global variables
    for dpVariableName in dpVariableNames:
        newText = newText.replace(dpVariableName,dpVariableName+'_f')

    newFile = open('sp_'+fileName,'w')
    newFile.write(newText)
    newFile.close()
    f.close()


 








     





