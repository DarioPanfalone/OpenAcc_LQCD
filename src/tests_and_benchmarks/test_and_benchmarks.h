#ifndef TESTS_AND_BENCHMARK_H_
#define TESTS_AND_BENCHMARK_H_


typedef struct test_info_t{

    int deoDoeIterations; 
    int multiShiftInverterRepetitions;
    double fakeShift;
    int benchmarkMode;
    int saveResults; 
    int parametersAreSet;




} test_info;

typedef struct gaugeMdTimes_t{

    double calcIpdotTimeBorder;
    double momSumMultTimeBorder;
    double momExpTimesConfTimeBorder;

    double calcIpdotTimeBulk;
    double momSumMultTimeBulk;
    double momExpTimesConfTimeBulk;

    double communicationsStartTime;
    double communicationsTime; // Measurement makes sense only withouy async communication.

    unsigned int count;

} gaugeMdTimeContainer;

typedef struct diracTimes_t{
 
    double totTransferTime; // Measurement makes sense only without async communications.
    unsigned int count;

} diracTimeContainer;

extern gaugeMdTimeContainer gauge_mdtimes;
extern gaugeMdTimeContainer gauge_mdtimes0;
extern diracTimeContainer dirac_times;

extern test_info  test_settings;

void gaugeMdCountersReset(gaugeMdTimeContainer*);
void diracCountersReset(diracTimeContainer*);












#endif
