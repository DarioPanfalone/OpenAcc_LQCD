module load pgi/15.9
module load cuda/7.0


export PGI_ACC_BUFFERSIZE=3221225472

rm prog_solo
rm *.o


pgcc  -O0  -c  Rand/random.c                   2> msg_err_0 #-Mcuda=maxregcount:128 
pgcc  -O0  -c  OpenAcc/include_all_main.c      2> msg_err_1 #-Mcuda=maxregcount:128 
pgcc  -O0  *.o -o prog_solo                    2> msg_err_3 #-Mcuda=maxregcount:128  
