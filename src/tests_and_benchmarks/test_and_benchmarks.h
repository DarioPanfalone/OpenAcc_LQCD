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

extern test_info  test_settings;













#endif
