module load pgi/15.9
module load cuda/7.0


export PGI_ACC_BUFFERSIZE=3221225472

rm prog_deb
rm *.o


pgcc    -c  Rand/random.c                   2> msg_err_0 #-Mcuda=maxregcount:128 
pgcc    -c  DbgTools/test_staples_SigmaPToSigma.c      2> msg_err_1 #-Mcuda=maxregcount:128 
pgcc    *.o -o prog_deb_SPtoS                    2> msg_err_3 #-Mcuda=maxregcount:128  
