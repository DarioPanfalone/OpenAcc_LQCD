#!/usr/bin/python

#modify this to change compiler, linker etc options.
PGIcls = 'COMPILER=pgcc\n\
COMPILER_FLAGS=-O3 -acc -Minfo=accel -v -ta=tesla:cc35,cuda7.0 \n\
LINKER_FLAGS=-acc  -Minfo=accel -O3 -v -ta=tesla:cc35,cuda7.0 \n' 

PGIclsSLOW ='COMPILER=pgcc\n\
COMPILER_FLAGS=-O0 \n\
LINKER_FLAGS=-O0 \n' 

GNUcls = 'COMPILER=gcc\n\
COMPILER_FLAGS=-O3 -std=c99\n\
LINKER_FLAGS=-lm\n' 

compiler_linker_settings = ''

main_files = [] # file containing main()

import os.path as path
from sys import exit,argv,stderr,stdout

# standard libraries that are not in the repo
# and should not be cheked by make


def header_scanner(filename):
    f = open(filename)
    filesfound = []
    standards = []
    for line in f.readlines():
        if 'main(' in line and ')' in line:
            main_files.append(filename)
        if '#include' in line :
            init_filename_position = line.find('<')
            end_filename_position = line.rfind('>')
            if init_filename_position == -1:
                init_filename_position = line.find('"')
                end_filename_position = line.rfind('"')
            if init_filename_position != -1 and end_filename_position != -1:
                filefound = line[init_filename_position+1:end_filename_position]
#                stderr.write("File found : " + line + ' \n')
                if filefound[0] is '.':
                    filesfound.append(filefound)
#                    stderr.write(filefound+\
#                            ' is in the package, adding it.\n')
                else:
#                   stderr.write(filefound +' is not in this package, not adding it.\n')
                    standards.append(filefound) # libraries like stdio.h
            else:
                stderr.write( "#include line incorrect:" + filename + ' ' + line + '\n')

    return filesfound,standards

class file_node:
    def get_all_dependences_raw(self):
        res = self.direct_dependences
        for son in self.sons:
            son_dependences = son.get_all_dependences_raw()
            if son_dependences is not None:
                res = res + son_dependences
        return res
    def get_all_standard_dependences_raw(self):
        res = self.direct_standard_dependences
        for son in self.sons:
            son_dependences = son.get_all_standard_dependences_raw()
            if son_dependences is not None:
                res = res + son_dependences
        return res

    def __init__(self,filename,tancestors):
        self.ancestors = tancestors
        self.direct_standard_dependences = [] #libraries like stdio.h
        self.direct_dependences_relative = []
        self.direct_dependences = []
        self.all_dependences = []
        self.all_standard_dependences = []
        self.sons = []
        self.sons_dict = dict()
        self.name = path.abspath(filename)
#        stderr.write("Scanning " + filename + "...\n" )
        self.direct_dependences_relative, self.direct_standard_dependences\
                = header_scanner(filename)
        for dependence in self.direct_dependences_relative:
            dependencem = path.abspath(path.dirname(self.name)+ '/'+dependence)
            if dependencem not in self.direct_dependences and\
                    dependencem not in self.ancestors: 
                self.direct_dependences.append(dependencem)
                son = file_node(dependencem, self.ancestors + [ self.name ] )
                self.sons.append(son)
                self.sons_dict[dependencem] = son
        self.all_dependences_raw = self.get_all_dependences_raw();
        self.all_standard_dependences_raw = \
                self.get_all_standard_dependences_raw();
        if self.all_standard_dependences_raw is not None:
            for dependence in self.all_standard_dependences_raw:
                if dependence not in self.all_standard_dependences:
                    self.all_standard_dependences.append(dependence)
        if self.all_dependences_raw is not None:
            for dependence in self.all_dependences_raw:
                if dependence not in self.all_dependences:
                    self.all_dependences.append(dependence)
    def showtree(self,n):
        prestring = ' ' * n
        for son,dependence in zip(self.sons,self.direct_dependences):
            print(prestring+dependence)
            son.showtree(n+1)
    def generate_make_string(self):
        makestring = ''
        if '.h' in self.name:
            makestring += self.name + " :"
        elif '.c' in self.name:
            makestring += path.basename(self.name)[:-2] + '.o : ' + self.name
        else:
            stderr.write("Filename " + self.name + " not valid.\n")
            return ''

        for dependence in self.direct_dependences:
            makestring += ' ' + dependence
        makestring += '\n\t'
        if '.h' in self.name:
            makestring += 'touch ' + self.name + '\n'
        elif '.c' in self.name:
            makestring += '$(COMPILER) -c $(COMPILER_FLAGS) ' +\
                    self.name + '\n\n'
        else:
            stderr.write("Filename " + self.name + " not valid.\n")
            return ''
        return makestring


def generate_makefile(targv):
    res = ''
    res += compiler_linker_settings + '\n'

    mainlinking_string = "main : " 
    for filename in targv[1:]:
        if '.c' in filename:
            mainlinking_string += ' ' + path.basename(filename)[:-2] + '.o'
        node = file_node(filename,[])
#        stderr.write(filename + '\n')
        makestring = node.generate_make_string()
        res += makestring
    mainlinking_string += '\n\t$(COMPILER) -o main *.o $(LINKER_FLAGS) \n\tif ! [ -d run ] ; then mkdir run; fi ; cp main run/\n'
    res += mainlinking_string
    
    for main_file in main_files:
        onefilecomp_string  = 'main_onefilecomp: random.o\n'
        onefilecomp_string += '\t$(COMPILER) $(COMPILER_FLAGS) -c '+\
                '-DONE_FILE_COMPILATION ' + main_file + ' \n'
        maino_name = path.basename(main_file)[:-2] + '.o'
        onefilecomp_string += '\t$(COMPILER) -o '+\
                'main_onefilecomp ' + maino_name + ' random.o $(LINKER_FLAGS)\n'
        onefilecomp_string += '\tcp main_onefilecomp run/\n'
        res += onefilecomp_string

    makeclean_string='clean:\n\trm -f *.o main main_onefilecomp\n'
    res += makeclean_string
    
    return res 

if __name__ == '__main__':

    if 'PGISLOW' not in argv and 'GNU' not in argv and 'PGI' not in argv:
        stderr.write("Please specify one compiler: either PGISLOW, GNU or PGI\n")
        exit(1)
    
    clsset = False
    if 'PGISLOW' in argv : 
        compiler_linker_settings = PGIclsSLOW;
        argv.remove('PGISLOW')
        clsset = True
    if 'GNU' in argv :
        if clsset :
            stderr.write("Please specify one compiler: either PGISLOW, GNU or PGI\n")
            argv.remove('GNU')
            exit(1)
        else:
            compiler_linker_settings = GNUcls;
            argv.remove('GNU')
            clsset = True
    if 'PGI' in argv :
        if clsset :
            stderr.write("Please specify one compiler: either PGISLOW, GNU or PGI\n")
            argv.remove('PGI')
            exit(1)
        else:
            compiler_linker_settings = PGIcls;
            argv.remove('PGI')
            clsset = True

    makefile = generate_makefile(argv)
    stdout.write(makefile)


