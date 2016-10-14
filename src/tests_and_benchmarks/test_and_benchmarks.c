#ifndef TESTS_AND_BENCHMARK_C_
#define TESTS_AND_BENCHMARK_C_

#include "./test_and_benchmarks.h"



test_info test_settings = {.deoDoeIterations = 1000, 
    .multiShiftInverterRepetitions =10,
    .fakeShift = 1e-14,
    .benchmarkMode = 1,
    .saveResults = 0,
    .parametersAreSet = 0};


gaugeMdTimeContainer gauge_mdtimes0={
    .calcIpdotTimeBorder = 0,
    .momSumMultTimeBorder = 0,
    .momExpTimesConfTimeBorder = 0,
    .calcIpdotTimeBulk = 0,
    .momSumMultTimeBulk = 0,
    .momExpTimesConfTimeBulk = 0,
    .communicationsStartTime = 0,
    .communicationsTime = 0,
    .count = 0};


gaugeMdTimeContainer gauge_mdtimes;


#endif
