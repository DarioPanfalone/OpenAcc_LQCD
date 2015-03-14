./compila_su_cudawn7.sh

while :
do
    sleep 1
    
    if [ -f successful_compilation ]
    then
	rm successful_compilation
	echo "Compilation.....   Yes"
	break
    fi
    
    if [ -f unsuccessful_compilation ]
    then
	rm unsuccessful_compilation
	echo "Compilation.....   No"
	break	
    fi
    
done
