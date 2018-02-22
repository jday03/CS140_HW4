#include "innerproduct.h"
#include <cilk/reducer_opadd.h>
#ifdef CILKPAR
#include <cilk.h>
#else
#define cilk_for for
#define cilk_main main
#define cilk_spawn
#define cilk_sync
#endif
static const int COARSENESS = 2;
double rec_cilkified(double *a, double *b, int n)
{
    double sum = 0;
	if(n == COARSENESS){
        for(int count = 0; count < n; ++count){
            sum += a[count] * b[count];
        }

    }
    else{
        int split1 = n/2;
        int split2 = n-split1;
        double sum1 = cilk_spawn rec_cilkified(a,b,split1);
        double sum2 = cilk_spawn rec_cilkified(a+split1,b+split1,split2);
        cilk_sync;
        sum = sum1 + sum2;
    }
	return sum;
}

double loop_cilkified(double *a, double *b, int n)
{
    int outerCountMax = n/COARSENESS;
    int extraValues = n%COARSENESS;
    double * miniSums= new double [outerCountMax];
    double sum = 0;

    cilk_for (int outerCount = 0; outerCount < outerCountMax; ++outerCount){
        miniSums[outerCount] = 0;
        for (int innerCount = 0; innerCount < COARSENESS; ++innerCount){
            miniSums[outerCount] += a[outerCount * COARSENESS + innerCount] * b[outerCount * COARSENESS + innerCount];

            // This loop is only for extra values due to n/COARSENESS rounding down with division.
            if(extraValues > 0 && outerCount == outerCountMax - 1 && innerCount == COARSENESS - 1){
                for(int extraCount = 0; extraCount < extraValues;++extraCount){
                    miniSums[outerCount] += a[outerCount * COARSENESS + innerCount + extraCount] * b[outerCount * COARSENESS + innerCount + extraCount];
                }

            }
        }

    }

    // summing the values
    for (int sumCount = 0; sumCount < outerCountMax; ++sumCount){
        sum += miniSums[sumCount];
    }

	return sum;
}

double hyperobject_cilkified(double *a, double *b, int n)
{
    int outerCountMax = n/COARSENESS;
    int extraValues = n%COARSENESS;
    double * miniSums= new double [outerCountMax];
    double sum = 0;

    cilk::reducer< cilk::opadd <double> > sum;
    cilk_for (int outerCount = 0; outerCount < outerCountMax; ++outerCount){

        for (int innerCount = 0; innerCount < COARSENESS; ++innerCount){
            sum += a[outerCount * COARSENESS + innerCount] * b[outerCount * COARSENESS + innerCount];

            // This loop is only for extra values due to n/COARSENESS rounding down with division.
            if(extraValues > 0 && outerCount == outerCountMax - 1 && innerCount == COARSENESS - 1){
                for(int extraCount = 0; extraCount < extraValues;++extraCount){
                    sum  += a[outerCount * COARSENESS + innerCount + extraCount] * b[outerCount * COARSENESS + innerCount + extraCount];
                }

            }
        }

    }

    // summing the values
    for (int sumCount = 0; sumCount < outerCountMax; ++sumCount){
        sum += miniSums[sumCount];
    }

    return sum;
}

