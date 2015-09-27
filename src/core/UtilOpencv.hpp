#ifndef RGM_UTILOPENCV_HPP_
#define RGM_UTILOPENCV_HPP_

#include "Common.hpp"

namespace RGM
{

class OpencvUtil
{
public:
    /// Visualize HOG
    static cv::Mat pictureHOG(cv::Mat_<Scalar> & filter, int bs);

    /// Rotates
    static cv::Mat rotate(cv::Mat m, float degree);

    /// Get subarray
    static cv::Mat subarray(cv::Mat img, cv::Rect roi, float padFactor, int padType);
    static cv::Mat subarray(cv::Mat img, cv::Rect roi, int padType);

}; // class OpencvUtil


template <typename T>
class OpencvUtil_
{
public:
    /// Resize a 3D matrix
    static cv::Mat_<T> resize(cv::Mat_<T> & m, float factor, int method);

}; // class OpencvUtil_

} // namespace RGM

#endif // RGM_UTILOPENCV_HPP_
