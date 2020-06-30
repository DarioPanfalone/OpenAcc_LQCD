```
   ____                    _____ _        _____  _      ______ 
  / __ \                  / ____| |      |  __ \| |    |  ____|
 | |  | |_ __   ___ _ __ | (___ | |_ __ _| |__) | |    | |__   
 | |  | | '_ \ / _ \ '_ \ \___ \| __/ _` |  ___/| |    |  __|  
 | |__| | |_) |  __/ | | |____) | || (_| | |    | |____| |____ 
  \____/| .__/ \___|_| |_|_____/ \__\__,_|_|    |______|______|
        | |                                                    
        |_|                                                    

        OpenACC STAggered Parallel LatticeQCD Everywhere

By E. Calore, M. Mesiti, F. Negro, G. Silvi 
(and others, which can be added if they feel offended)
```
# An humorous intro to this code

## 0. Create the `configure` and `Makefile.in`s in the repo
This step is done with 
```
autoreconf --install
```
This will generate a number of files (among which the configure and the
Makefile.in files), starting from the Makefile.am in the various directories
and configure.ac.

## 1. Build - General ideas

At build time the geometry of the lattice must be known, and must be passed to 
the "configure" script (notice that there is also a "configure_wrapper" script, 
which is a python script, which makes easier to call "configure": more on that is written
below). At present, "configure" reads a file named "geom_defines.txt" 
which contains the lattice dimensions and the dimensions of the OpenAcc tiles.
A commented example of the file can be found in the 'doc' directory.

In case the system you are working on uses environment modules, you need to load them before 
invoking "configure".

Example:
```
mkdir build
cd build

<add the geom_defines.txt file>

../configure_wrapper pgi cc35 cuda8.0 # if you like the configure_wrapper script

make && make install  # you may want to use make -j32 to be faster.
```
### 1.1 Usage of the plain "configure" script

You must also pass to configure:
- The compiler : `CC=<the compiler you want to use>`
- prefix : The place where you want the code to be installed. Note that it defaults to 
  a system directory you're likely not to have access to, but even if you have access to it 
  you should not install this software in system directories.
- CFLAGS : compiler flags. Will be discussed later, depends on the compiler used.
- LDFLAGS : the linker flags

### 1.2 The `configure_wrapper` script

This simple wrapper allows the user to call call configure without remembering all the
options. If you don't trust this wrapper, you can look into the code and look at the 
options you need.
It is called, e.g., in this way:
```
../configure_wrapper pgi cc35 cuda8.0
```
where you must change the options according to your goals. If invoked in ths way, it will 
produce the command
```
./configure CC=pgcc CFLAGS="-acc=noautopar -v -O3 -Minfo=all -ta=tesla:cc35,cuda8.0 -DUSE_MPI_CUDA_AWARE -I${MPIINC}" LDFLAGS="-acc=noautopar -v -O3 -Minfo=all -ta=tesla:cc35,cuda8.0 -lmpi -L${MPILIB}"  --prefix=$(pwd)
```
and you will be asked if you want to execute it or not. Note that, for simplicity, 
the install directory (`--prefix`) is set to the working directory (the programs will be
installed in `$PWD/bin`).  If you want to change that, just specify 
`--prefix=<your explicit installing directory>`
After that, you need to run 
```
make -jN # with N > 1 to make it faster, 
```
And
```
make install 
```
to finally install in the `--prefix` directory (see above).

## 2. Benchmarks

In order to run all the necessary benchmarks, it is advised to create an aptly named 
directory. In such directory, create a "geom_defines.txt" file (see an example in 'docs').
Then build the software, using the commands described in the previous section.
You can also produce all the necessary slurm scripts and setting files using the script
'prepare_slurm_benchmarks.sh' in 'scripts'. For example,
```
mkdir my_benchmark
cd my_benchmark
<create the geom_defines.txt file>
../configure_wrapper pgi cc35 cuda8.0
make && make install 
../scripts/prepare_tbps.sh geom_defines.txt benchmark
```
Notice that 'scripts/prepare_tpbs.sh' assumes the executables are in ./bin.

## 3. Documentation
Do you want to know more? A 26-page latex document is waiting for you.
How to create it:
```
cd doc
pdflatex general_guide.tex
```
