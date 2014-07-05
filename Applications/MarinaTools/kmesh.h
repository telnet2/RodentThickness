//
//  kmesh.h
//  ktools
//
//  Created by Joohwi Lee on 12/5/13.
//
//

#ifndef __ktools__kmesh__
#define __ktools__kmesh__

#include <iostream>

namespace std {
    template <class T>
    class MeanStd {
    public:
        int n;
        double sum, sum2;
        MeanStd(): n(0), sum(0), sum2(0) {}
        void operator()(T x) {
            sum += x;
            sum2 += x*x;
            n++;
        }
        double mean() {
            return sum / n;
        }
        double std() {
            double m = mean();
            return (sum2) / n - m*m;
        }
    };
}

#endif /* defined(__ktools__kmesh__) */
