#ifndef INVERTER_TRICKS_H_
#define INVERTER_TRICKS_H_

typedef struct inv_tricks_t{

	int singlePInvAccelMultiInv;
	int useMixedPrecision;
	double mixedPrecisionDelta;
	int restartingEvery;

} inv_tricks;

extern inv_tricks inverter_tricks;


#endif
