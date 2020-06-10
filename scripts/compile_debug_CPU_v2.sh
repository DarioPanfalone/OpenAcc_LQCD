module load pgi/15.9
module load cuda/7.0


export PGI_ACC_BUFFERSIZE=3221225472

rm prog_deb_v2
rm *.o


pgcc  -v -c  Rand/random.c                   2> msg_err_0 #-Mcuda=maxregcount:128 
pgcc  -v -c  OpenAcc/fermion_force_debugger_main_v2.c      2> msg_err_1 #-Mcuda=maxregcount:128 
pgcc  -v  *.o -o prog_deb_v2                    2> msg_err_3 #-Mcuda=maxregcount:128  
