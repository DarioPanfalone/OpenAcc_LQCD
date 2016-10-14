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
    double communicationsTime;

    int    count;

} gaugeMdTimeContainer;


extern gaugeMdTimeContainer gauge_mdtimes;
extern gaugeMdTimeContainer gauge_mdtimes0;

extern test_info  test_settings;













#endif
