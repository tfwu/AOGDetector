#ifndef RGM_UTILMATH_HPP_
#define RGM_UTILMATH_HPP_

#include <math.h>
#include "opencv2/opencv.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace RGM
{

template<typename T>
class MathUtil_
{

public:
    static T                  median(std::vector<T> & vData, int num=std::numeric_limits<int>::max(), bool sortData=false);
    static std::vector<float> pdist(const std::vector<cv::Point_<T> > & pts, int num=std::numeric_limits<int>::max());
    static std::vector<T>     linspace(T s, T e, T interval);
    static std::vector<T>     hist(std::vector<T> & y, std::vector<T> & x);
    static std::vector<T>     convnSame(std::vector<T> & x, std::vector<T> & filter);

    static T calcErr(std::vector<T> & pscores, std::vector<T> & nscores, T & thr) ;
    static T calcVar(std::vector<T> & pscores, T & thr) ;
}; // MathUtil_

} // namespace RGM

#endif // RGM_UTILMATH_HPP_
