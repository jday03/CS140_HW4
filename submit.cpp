#include "innerproduct.h"
#include <cilk/reducer_opadd.h>
static const int COARSENESS = 3000;

double rec_cilkified(double *a, double *b, int n)
{
    double sum = 0;
	if(n <= COARSENESS){
        for(int count = 0; count < n; ++count){
            sum += a[count] * b[count];
        }

    }
    else{
        int split1 = n/2;
        int split2 = n-split1;
        double sum1 = cilk_spawn rec_cilkified(a,b,split1);
        double sum2 = cilk_spawn rec_cilkified((a+split1),(b+split1),split2);
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

                for(int extraCount = 1; extraCount < extraValues + 1;++extraCount){
                    miniSums[outerCount] += a[outerCount * COARSENESS + innerCount + extraCount] * b[outerCount * COARSENESS + innerCount + extraCount];
                }

            }
        }

    }


// if coarseness > n
    if(outerCountMax == 0){
        double sum = 0;
        for (int innerCount = 0; innerCount < n; ++innerCount){
            sum += a[innerCount] * b[innerCount];

            // This loop is only for extra values due to n/COARSENESS rounding down with division.
            if(extraValues > 0 && innerCount == n - 1){
                for(int extraCount = 1; extraCount < extraValues + 1;++extraCount){
                    sum  += a[innerCount + extraCount] * b[innerCount + extraCount];
                }
            }
        }
        return sum;
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

    cilk::reducer< cilk::op_add<double> > sum1;
    cilk_for (int outerCount = 0; outerCount < outerCountMax; ++outerCount){

        for (int innerCount = 0; innerCount < COARSENESS; ++innerCount){
            *sum1 += a[outerCount * COARSENESS + innerCount] * b[outerCount * COARSENESS + innerCount];

            // This loop is only for extra values due to n/COARSENESS rounding down with division.
            if(extraValues > 0 && outerCount == outerCountMax - 1 && innerCount == COARSENESS - 1){
                for(int extraCount = 1; extraCount < extraValues + 1;++extraCount){
                    *sum1  += a[outerCount * COARSENESS + innerCount + extraCount] * b[outerCount * COARSENESS + innerCount + extraCount];
                }

            }
        }

    }


    if(outerCountMax == 0){
        for (int innerCount = 0; innerCount < n; ++innerCount){
            *sum1 += a[innerCount] * b[innerCount];

            // This loop is only for extra values due to n/COARSENESS rounding down with division.
            if(extraValues > 0 && innerCount == n - 1){
                for(int extraCount = 1; extraCount < extraValues + 1;++extraCount){
                    *sum1  += a[innerCount + extraCount] * b[innerCount + extraCount];
                }

            }
        }
    }

    return sum1.get_value();


}

