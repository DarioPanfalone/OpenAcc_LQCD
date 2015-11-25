#!/usr/bin/python

#modify this to change compiler, linker etc options.
compiler_linker_settings = \
        'COMPILER=gcc\n\
        COMPILER_FLAGS=-O3 -std=c99\n\
        LINKER_FLAGS=\"-lm\"\n' 




import os.path as path
from sys import exit,argv,stderr,stdout

# standard libraries that are not in the repo
# and should not be cheked by make

silent = True

def header_scanner(filename):
    f = open(filename)
    filesfound = []
    for line in f.readlines():
        if '#include' in line:
            init_filename_position = line.find('<')
            end_filename_position = line.rfind('>')
            if init_filename_position == -1:
                init_filename_position = line.find('"')
                end_filename_position = line.rfind('"')
            if init_filename_position != -1 and end_filename_position != -1:
                filefound = line[init_filename_position+1:end_filename_position]
                stderr.write("File found : " + line + ' \n')
                if filefound[0] is '.':
                    filesfound.append(filefound)
                    stderr.write(filefound+\
                            ' is in the package, adding it.\n')
                else:
                    stderr.write(filefound +\
                            ' is not in this package, not adding it.\n')
            else:
                stderr.write( "#include line incorrect:" + filename + ' ' + line + '\n')

    return filesfound

class file_node:
    def get_all_dependences_raw(self):
        res = self.direct_dependences
        for son in self.sons:
            son_dependences = son.get_all_dependences_raw()
            if son_dependences is not None:
                res = res + son_dependences
        return res

    def __init__(self,filename):
        self.direct_dependences = []
        self.all_dependences = []
        self.sons = []
        self.name = filename
        stderr.write("Scanning " + filename + "...\n" )
        self.direct_dependences = header_scanner(filename)
        for dependence in self.direct_dependences:
            dependencem = path.dirname(self.name)+ '/'+dependence
            son = file_node(dependencem)
            self.sons.append(son)
        self.all_dependences_raw = self.get_all_dependences_raw();
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
            dependencem = path.dirname(self.name)+ '/'+dependence
            makestring += ' ' + dependencem 
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


if __name__ == '__main__':
    stdout.write(compiler_linker_settings)

    mainlinking_string = "main : " 
    for filename in argv[1:]:
        if '.c' in filename:
            mainlinking_string += ' ' + path.basename(filename)[:-2] + '.o'
        node = file_node(filename)
        stderr.write(filename + '\n')
        makestring = node.generate_make_string()
        stdout.write(makestring)
    mainlinking_string += '\n\t$(COMPILER) -o main $(LINKER_FLAGS) *.o\n'
    stdout.write(mainlinking_string)
    makeclean_string='clean:\n\trm -f *.o main\n'
    stdout.write(makeclean_string)




