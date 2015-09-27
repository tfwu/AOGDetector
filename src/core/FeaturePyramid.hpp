/// This file is adapted from FFLDv2 (the Fast Fourier Linear Detector version 2)
/// Copyright (c) 2013 Idiap Research Institute, <http://www.idiap.ch/>
/// Written by Charles Dubout <charles.dubout@idiap.ch>

#ifndef RGM_FEATUREPYRAMID_HPP_
#define RGM_FEATUREPYRAMID_HPP_

#include "Common.hpp"
#include "UtilString.hpp"
#include "UtilLog.hpp"

namespace RGM
{

#ifndef RGM_USE_PCA_DIM
#define RGM_USE_PCA_DIM  5
#endif

class FeaturePyramid
{
public:
    /// Number of features (guaranteed to be even). Fixed at compile time for both ease of use
    /// and optimal performance.
#ifdef RGM_USE_CAFFE
    static const int NbFeatures = 96; // learned filter kernel at 1st layer
#else
#if (!defined RGM_USE_EXTRA_FEATURES) || (defined RGM_USE_FELZENSZWALB_HOG_FEATURES)
    static const int NbFeatures = 32; // HOG
#else
    static const int NbFeatures = 48; // HOG + LBP + Color
#endif
#endif

    //static const int MinLevelSz = 5;

    /// Type of a pyramid level cell (fixed-size array of length NbFeatures).
    typedef Eigen::Array<Scalar, NbFeatures, 1> Cell;

    /// Type of a pyramid level (matrix of cells).
    typedef Eigen::Matrix<Cell, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Level;

    /// Double precision for gradient used in optimization
#ifndef RGM_USE_DOUBLE
    typedef Eigen::Array<double, NbFeatures, 1>                                    dCell;
    typedef Eigen::Matrix<dCell, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>  dLevel;
#else
    typedef Cell    dCell;
    typedef Level   dLevel;
#endif

    /// Types for PCA
#if RGM_USE_PCA_DIM
    static const int NbPCAFeatures = RGM_USE_PCA_DIM + 1;
    typedef Eigen::Array<Scalar, NbPCAFeatures, 1> PCACell;
    typedef Eigen::Matrix<PCACell, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> PCALevel;
#endif

    /// Constructs an empty pyramid. An empty pyramid has no level.
    FeaturePyramid();

    /// Constructs a pyramid from the image of a Scene.
    /// @param[in] image The image of the Scene.
    /// @param[in] padx Amount of horizontal zero padding (in cells).
    /// @param[in] pady Amount of vertical zero padding (in cells).
    /// @param[in] interval Number of levels per octave in the pyramid.
    /// @note The amount of padding and the interval should be at least 1.
    FeaturePyramid(const cv::Mat & image, int cellSize, int padx, int pady, int octave = 0, int interval = 5, bool extraOctave = false);

    FeaturePyramid(const cv::Mat & image, int cellSize, int padx, int pady, Cell & bgmu,
                   int octave = 0, int interval = 5, bool extraOctave = false);

    /// Constructs a pyramid from parameters and a list of levels and scales [and PCAlevels].
    /// @param[in] cellSize Size of a cell (in pixels)
    /// @param[in] padx Amount of horizontal zero padding (in cells).
    /// @param[in] pady Amount of vertical zero padding (in cells).
    /// @param[in] octave Number of octaves in the pyramid
    /// @param[in] interval Number of levels per octave in the pyramid.
    /// @param[in] extraOctave If an extra octave is used
    /// @param[in] levels List of pyramid levels.
    /// @param[in] PCAlevels (optional) List of pyramid PCA projected levels
    /// @param[in] scales List of scales of pyramid levels
    /// @note The amount of padding and the interval must both be at least 1.
    /// @note The input levels are swapped with empty ones on exit.
#if RGM_USE_PCA_DIM
    FeaturePyramid(int cellSize, int padx, int pady, int octave, int interval, bool extraOctave, int imgWd, int imgHt,
                   std::vector<Level> & levels, std::vector<PCALevel> & PCAlevels, std::vector<Scalar> & scales);
#endif

    FeaturePyramid(int cellSize, int padx, int pady, int octave, int interval, bool extraOctave, int imgWd, int imgHt,
                   std::vector<Level> & levels, std::vector<Scalar> & scales);


    /// Returns whether the pyramid is empty. An empty pyramid has no level.
    bool empty() const;

#if RGM_USE_PCA_DIM
    bool emptyPCA() const;
#endif

    /// Returns the cell size
    int cellSize() const;

    /// Returns the amount of horizontal zero padding (in cells).
    int padx() const;

    /// Returns the amount of vertical zero padding (in cells).
    int pady() const;

    /// Returns the number of octaves in the pyramid
    int octave() const;

    /// Returns the number of levels per octave in the pyramid.
    int interval() const;

    /// Returns if an extra octave is used
    bool extraOctave() const;

    /// Returns the pyramid levels.
    /// @note Scales are given by the following formula: 2^(1 - @c index / @c interval).
    const std::vector<Level> & levels() const;

#if RGM_USE_PCA_DIM
    /// Return the PCA projected pyramid levels
    const std::vector<PCALevel> & PCAlevels() const;

    /// ProjectS the pyramid using PCA coef @p coef
    void project(const Matrix & coef);
#endif

    /// Returns the scales
    const std::vector<Scalar> & scales() const;

    /// Returns the status of levels
    const std::vector<bool> & validLevels() const;
    std::vector<bool> & getValidLevels();

    /// Returns the indice of valid levels
    const std::vector<int> & idxValidLevels();

    /// Returns the number of levels
    int nbLevels() const;

    /// Returns the number of valid levels
    int nbValidLevels() const;

    /// Returns the original image size
    int imgWd() const;
    int imgHt() const;

    /// Returns the convolutions of the pyramid with a filter.
    /// @param[in] filter Filter.
    /// @param[out] convolutions Convolution of each level.
    void convolve(const Level & filter, std::vector<Matrix> & convolutions) const;

    /// Returns the flipped version (horizontally) of a level.
    static Level Flip(const FeaturePyramid::Level & level);

    /// Maps a pyramid level to a simple matrix (useful to apply standard matrix operations to it).
    /// @note The size of the matrix will be rows x (cols * NbFeatures).
    static Eigen::Map<Matrix, Eigen::Aligned> Map(Level & level);

    static Eigen::Map<dMatrix, Eigen::Aligned> dMap(dLevel & level);

    /// Maps a const pyramid level to a simple const matrix (useful to apply standard matrix
    /// operations to it).
    /// @note The size of the matrix will be rows x (cols * NbFeatures).
    static const Eigen::Map<const Matrix, Eigen::Aligned> Map(const Level & level);

    /// Converts a pyramid level ta a cv::Mat_
    /// @param[in] startDim Index of starting dimension (out of the total dimension, i.e. NbFeatures)
    /// @param[in] endDim Index of ending dimension
    static cv::Mat_<Scalar> convertToMat(const Level & level, int startDim, int endDim);

    /// Returns the contrast insensitive orientations
    static cv::Mat_<Scalar> fold(const Level & level);

    /// Resizes a level
    static void resize(const Level & in, Scalar factor,  Level & out);

    /// Visualizes a pyramid level (HOG features)
    static cv::Mat visualize(const Level & level, int bs=20);

    /// Computes the virtual padding
    static int VirtualPadding(int padding, int ds);

    /// Efficiently computes Histogram of Oriented Gradient (HOG) features [+ LBP + Color]
    /// Code to compute HOG features as described in "Object Detection with Discriminatively Trained
    /// Part Based Models" by Felzenszwalb, Girshick, McAllester and Ramanan, PAMI 2010
    static void computeHOGFeature(const cv::Mat & image, Level & level, Cell & bgmu,
                                  int padx = 1, int pady = 1, int cellSize = 8);

private:
    /// Compute the feature pyramid
    void computePyramid(const cv::Mat & image, int cellSize, int padx, int pady,
                        Cell & bgmu, int octave, int interval, bool extraOctave);

    /// Computes the 2D convolution of a pyramid level with a filter
    static void Convolve(const Level & x, const Level & y, Matrix & z);    

private:
    int cellSize_;
    int padx_;
    int pady_;
    int interval_;
    int octave_;
    bool extraOctave_;

    Cell mu_; // background mean for padding

    std::vector<Level> levels_;
#if RGM_USE_PCA_DIM
    std::vector<PCALevel> PCAlevels_;
#endif
    std::vector<Scalar> scales_;

    std::vector<bool> validLevels_;
    std::vector<int>  idxValidLevels_;

    int imgWd_;
    int imgHt_;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

    DEFINE_RGM_LOGGER;

}; // class FeaturePyramid

} // namespace RGM


/// Some compilers complain about the lack of a NumTraits for Eigen::Array<Scalar, NbFeatures, 1>
namespace Eigen
{
template <>
struct NumTraits<Array<RGM::Scalar, RGM::FeaturePyramid::NbFeatures, 1> > :
        GenericNumTraits<Array<RGM::Scalar, RGM::FeaturePyramid::NbFeatures, 1> > {
    static inline RGM::Scalar dummy_precision()
    {
        return 0; // Never actually called
    }
};

#ifndef RGM_USE_DOUBLE
template <>
struct NumTraits<Array<double, RGM::FeaturePyramid::NbFeatures, 1> > :
        GenericNumTraits<Array<double, RGM::FeaturePyramid::NbFeatures, 1> > {
    static inline double dummy_precision()
    {
        return 0; // Never actually called
    }
};
#endif

#if RGM_USE_PCA_DIM
template <>
struct NumTraits<Array<RGM::Scalar, RGM::FeaturePyramid::NbPCAFeatures, 1> > :
        GenericNumTraits<Array<RGM::Scalar, RGM::FeaturePyramid::NbPCAFeatures, 1> > {
    static inline RGM::Scalar dummy_precision()
    {
        return 0; // Never actually called
    }
};
#endif

} // namespace Eigen

/// Serialize FeaturePyramid::Level
namespace boost
{
namespace serialization
{
template <class Archive>
void save(Archive & ar, const RGM::FeaturePyramid::Level & m, const unsigned int version)
{
    int rows=m.rows(),cols=m.cols();
    ar & BOOST_SERIALIZATION_NVP(rows);
    ar & BOOST_SERIALIZATION_NVP(cols);

    int dim = RGM::FeaturePyramid::NbFeatures;
    ar & BOOST_SERIALIZATION_NVP(dim);

    for ( int i = 0; i < rows; ++i ) {
        for ( int j = 0; j < cols; ++j ) {
            ar & make_array(m(i,j).data(), dim);
        }
    }

//    for ( int i = 0; i < rows; ++i ) {
//        std::string str = "row"+ RGM::NumToString_<int>(i);
//        ar & make_nvp(str.c_str(),
//                      make_array(reinterpret_cast<const char *>(m.row(i).data()), cols * sizeof(RGM::FeaturePyramid::Cell)));
//    }

}

template <class Archive>
void load(Archive & ar, RGM::FeaturePyramid::Level & m, const unsigned int version)
{
    int rows,cols, dim;

    ar & BOOST_SERIALIZATION_NVP(rows);
    ar & BOOST_SERIALIZATION_NVP(cols);
    ar & BOOST_SERIALIZATION_NVP(dim);

    assert(dim == RGM::FeaturePyramid::NbFeatures);

    m = RGM::FeaturePyramid::Level::Constant(rows, cols, RGM::FeaturePyramid::Cell::Zero());

    for ( int i = 0; i < rows; ++i ) {
        for ( int j = 0; j < cols; ++j ) {
            ar & make_array(m(i,j).data(), dim);
        }
    }

//    for ( int i = 0; i < rows; ++i ) {
//        std::string str = "row"+ RGM::NumToString_<int>(i);
//        ar & make_nvp(str.c_str(),
//                      make_array(reinterpret_cast<char *>(m.row(i).data()), cols * sizeof(RGM::FeaturePyramid::Cell)));
//    }
}

template <class Archive>
void serialize(Archive & ar, RGM::FeaturePyramid::Level & m, const unsigned int version)
{
    split_free(ar,m,version);
}


} // namespace serialization
} // namespace boost

#endif // RGM_FEATUREPYRAMID_HPP_
