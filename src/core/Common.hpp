#ifndef RGM_COMMON_HPP_
#define RGM_COMMON_HPP_

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/lexical_cast.hpp>

#include <Eigen/Core>

#include <opencv2/opencv.hpp>

namespace RGM
{
/// Disable the copy and assignment operator for a class.
#define DISABLE_COPY_AND_ASSIGN(classname) \
    private:\
    classname(const classname&);\
    classname& operator=(const classname&)

/// Instantiate serialization
#define INSTANTIATE_BOOST_SERIALIZATION(classname) \
    template  void classname::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar, const unsigned int version); \
    template  void classname::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar, const unsigned int version); \
    template  void classname::serialize<boost::archive::text_iarchive>(boost::archive::text_iarchive& ar, const unsigned int version); \
    template  void classname::serialize<boost::archive::text_oarchive>(boost::archive::text_oarchive& ar, const unsigned int version);


/// Type of a Scalar
#ifdef RGM_USE_DOUBLE
typedef double Scalar;
#else
typedef float Scalar;
#endif

/// Type of a complex value.
typedef std::complex<Scalar> CScalar;

/// Type of a color image pixel
typedef cv::Vec<Scalar, 3> Pixel;

/// Type of a int matrix
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXi;

/// Type of a matrix.
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

/// Type of a complex matrix.
typedef Eigen::Matrix<CScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> CMatrix;

/// Double precision for gradient used in optimization
#ifndef RGM_USE_DOUBLE
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> dMatrix;
#else
typedef Matrix  dMatrix;
#endif

/// round
#define ROUND(x) \
    boost::math::round(x)

}  // namespace RGM

#endif  // RGM_COMMON_HPP_
