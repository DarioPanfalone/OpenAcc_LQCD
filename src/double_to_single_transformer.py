#!/usr/bin/python2 
# This script does some work for converting 

from sys import argv
import os 


fileNames = [\
'OpenAcc/alloc_vars.c',\
'OpenAcc/alloc_vars.h',\
'OpenAcc/cayley_hamilton.h',\
'OpenAcc/fermion_force.c',\
'OpenAcc/fermion_force.h',\
'OpenAcc/fermion_force_utilities.c',\
'OpenAcc/fermion_force_utilities.h',\
'OpenAcc/fermionic_utilities.c',\
'OpenAcc/fermionic_utilities.h',\
'OpenAcc/fermion_matrix.c',\
'OpenAcc/fermion_matrix.h',\
'OpenAcc/inverter_full.c',\
'OpenAcc/inverter_full.h',\
'OpenAcc/inverter_multishift_full.c',\
'OpenAcc/inverter_multishift_full.h',\
'OpenAcc/ipdot_gauge.c',\
'OpenAcc/ipdot_gauge.h',\
'OpenAcc/matvecmul.h',\
'OpenAcc/md_integrator.c',\
'OpenAcc/md_integrator.h',\
'OpenAcc/plaquettes.c',\
'OpenAcc/plaquettes.h',\
'OpenAcc/rettangoli.c',\
'OpenAcc/rettangoli.h',\
'OpenAcc/single_types.h',\
'OpenAcc/stouting.c',\
'OpenAcc/stouting.h',\
'OpenAcc/struct_c_def.h',\
'OpenAcc/su3_measurements.c',\
'OpenAcc/su3_measurements.h',\
'OpenAcc/su3_utilities.c',\
'OpenAcc/su3_utilities.h',\
'DbgTools/dbgtools.h',\
'DbgTools/dbgtools.c',\
'Mpi/communications.h',\
'Mpi/communications.c',\
]

# files where we can change all 'd_complex' with 'f_complex'
# without fears
filesALLDtoF=[\
'OpenAcc/alloc_vars.c',\
'OpenAcc/alloc_vars.h',\
'OpenAcc/cayley_hamilton.h',\
'OpenAcc/fermion_matrix.c',\
'OpenAcc/fermion_matrix.h',\
'OpenAcc/ipdot_gauge.c',\
'OpenAcc/ipdot_gauge.h',\
'OpenAcc/matvecmul.h',\
'OpenAcc/plaquettes.c',\
'OpenAcc/plaquettes.h',\
'OpenAcc/rettangoli.c',\
'OpenAcc/rettangoli.h',\
'OpenAcc/single_types.h',\
'OpenAcc/stouting.c',\
'OpenAcc/stouting.h',\
'OpenAcc/struct_c_def.h',\
'OpenAcc/su3_measurements.c',\
'OpenAcc/su3_measurements.h',\
'OpenAcc/su3_utilities.c',\
'OpenAcc/su3_utilities.h',
'Mpi/communications.h',\
'Mpi/communications.c',\
'DbgTools/dbgtools.h',\
'DbgTools/dbgtools.c',\
]




dpFunctionNames = [] # new function names will just be dp function names + '_f' at the end
dpVariableNames = [] # new function names will just be dp function names + '_f' at the end
# as in struct_c_def.h
dpTypes = ['double_soa','dcomplex_soa','vec3_soa','vec3','su3_soa','thmat_soa','tamat_soa',\
        'global_vec3_soa','global_su3_soa']
# corresponding types in in sp_struct_c_def.h
spTypes = ['float_soa','fcomplex_soa','vec3_f_soa','vec3_f','su3_soa_f','thmat_soa_f',\
        'tamat_soa_f','global_vec3_f_soa','global_su3_soa_f']





dpToSpDict = dict(zip(dpTypes,spTypes))

# this function is needed to populate the list 'dpFunctionNames'


def findFunctionNames(lineRaw):
    returnTypes=['int','double','void','d_complex']
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


dpFunctionNames.sort(key = len, reverse = True )

functionNamesFile = open('functions_found.txt','w')
for foundFunction in dpFunctionNames:
    functionNamesFile.write(foundFunction + '\n')
functionNamesFile.close()
print "Total functions found: ", len(dpFunctionNames)

# collecting all global variable names
print "\n\n-----------------\nLooking for variable names\n-----------------\n"
f = open('OpenAcc/alloc_vars.c')
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

dpVariableNames.sort(key = len , reverse = True)

globalVarNamesFile = open('glvar_found.txt','w')
for foundGlobalVar in dpVariableNames:
    globalVarNamesFile.write(foundGlobalVar + '\n')
globalVarNamesFile.close()

dpVariableNames.append('phases') # for U1 used in dirac matrix

print "Total global variables found: ", len(dpVariableNames)



# collected all double precision function names, now just replacing:
# - all function names
# - all relevant types

ans = '' # answer for 'overwrite file
for fileName in fileNames:
    if (len(argv) > 1 and fileName in argv) or len(argv)==1:
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
            newText = newText.replace(os.path.basename(fileName2),\
                    'sp_'+os.path.basename(fileName2))
        #taking care of header guards
        headerGuard = os.path.basename(fileName).upper().replace('.','_')
        if headerGuard not in newText:
            print "Warning, header guard \'" + headerGuard + "\' not found in file " + fileName
        newText = newText.replace(headerGuard, "SP_"+headerGuard)
        # taking care of global variables
        for dpVariableName in dpVariableNames:
            newText = newText.replace(dpVariableName,dpVariableName+'_f')
    
        newText = newText.replace('deltas_Omelyan','deltas_Omelyan_f')
    
    
        # it may happen that two transformations appear on the same symbol, 
        # and an '_f_f' is appended instead of just '_f'
        newText = newText.replace('_f_f','_f')
        newText = newText.replace('_f_f','_f')
    
    
        if fileName in filesALLDtoF:
            newText = newText.replace('d_complex','f_complex')
            newText = newText.replace('double','float')
            newText = newText.replace('%lf','%f')
        
    
    
    
        newFileName = os.path.dirname(fileName)+'/sp_'+os.path.basename(fileName)
        writeIt = True
        
        if os.path.exists(newFileName) and ans != 'a':
            ans = ''
            while ans not in ['y','n','a']:
                print "Overwrite file \'"+newFileName+"\'? (y=yes,n=no,a=yes to all)"
                ans = raw_input().lower()
                if ans == 'n':
                    writeIt = False
        if writeIt : 
            print "Writing file ", newFileName , " ..."
            newFile = open(newFileName,'w')
            newFile.write(newText)
            newFile.close()
        f.close()


 


print "REMEMBER TO ADD \n\ntypedef double complex d_complex\n\nto sp_struct_c_def.h !!"





     





