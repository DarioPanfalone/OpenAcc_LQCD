module load pgi/15.9
module load cuda/7.0


export PGI_ACC_BUFFERSIZE=3221225472

rm test_stout_ingredients
rm *.o


#pgcc    -c  Rand/random.c                   2> msg_err_0 #-Mcuda=maxregcount:128 
pgcc    -c  DbgTools/test_stout_ingredients.c      2> msg_err_1 #-Mcuda=maxregcount:128 
pgcc    *.o -o test_stout_ingredients               2> msg_err_3 #-Mcuda=maxregcount:128  
