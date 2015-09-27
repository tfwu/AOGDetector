#ifndef RGM_PARAMETERS_HPP_
#define RGM_PARAMETERS_HPP_

#include <string>
#include <vector>

#include "FeaturePyramid.hpp"
#include "Rectangle.hpp"


namespace RGM
{
/// The ParamUtil class consists of the common utility parameters used in learning the grammar
class ParamUtil
{
public:
    /// Default constructor
    ParamUtil();

    /// Default destructor
    virtual ~ParamUtil();

    /// Returns the learning status of the offset parameter
    Scalar   learningStatus() const;
    Scalar & getLearningStatus();

    /// Returns the regularization cost of the offset parameter
    Scalar   regularizationCost() const;
    Scalar & getRegularizationCost();

protected:
    Scalar learningStatus_; // 0: parameters are not learned; >0: they are learned
    Scalar regularizationCost_;

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};


/// The Appearance class represents the appearance parameters of a TERMINAL-node in the grammar or
/// the features in the feature pyramid of a deformed TERMINAL-node instantiated in the parse tree
class Appearance : public ParamUtil
{
public:
    /// Type of parameter vector
    typedef FeaturePyramid::Level Param;

    typedef FeaturePyramid::dLevel dParam;

    /// Type of cell
    typedef FeaturePyramid::Cell Cell;

    typedef FeaturePyramid::dCell dCell;

#if RGM_USE_PCA_DIM
    /// Type of PCA-projected parameter vector
    typedef FeaturePyramid::PCALevel PCAParam;

    /// Type of PCA-projected cell
    typedef FeaturePyramid::PCACell PCACell;
#endif

    /// Default constuctor
    Appearance();

    /// Copy constructor
    Appearance(const Appearance & app);

    /// Default destructor
    virtual ~Appearance();

    /// Init
    void init(int wd, int ht);

    /// Returns the param
    const Param & w() const;
    Param & getW();

    /// Returns the lower bound
    const Param & lowerBound() const;
    Param & getLowerBound();

    /// Returns the gradient
    const dParam & gradient() const;
    dParam & getGradient();

private:
    Param  w_;
    Param  lowerBound_;

    dParam  gradient_; // gradient used in learning

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

}; // class Appearance


/// The Offset class represents the bias term
class Offset : public ParamUtil
{
public:
    /// Default constructor
    Offset();

    /// Copy constructor
    explicit Offset(const Offset & off);

    /// Default destructor
    virtual ~Offset();

    /// Returns the offset parameter
    Scalar   w() const;
    Scalar & getW();

    /// Returns the parameter lower bound
    Scalar   lowerBound() const;
    Scalar & getLowerBound();

    /// Returns the gradient
    const double & gradient() const;
    double & getGradient();

private:
    Scalar  w_;
    Scalar  lowerBound_;

    double  gradient_;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

}; // class Offset



/// The ScalePrior class represents the scale prior parameter
class Scaleprior : public ParamUtil
{
public:
    /// Type of parameter vector
    typedef Eigen::Matrix<Scalar, 1, 3, Eigen::RowMajor> Param;

    /// Type of parameter gradient
#ifndef RGM_USE_DOUBLE
    typedef Eigen::Matrix<double, 1, 3, Eigen::RowMajor> dParam;
#else
    typedef Param  dParam;
#endif

    /// Type of a row vector
    typedef Eigen::Matrix<Scalar, 1, Eigen::Dynamic, Eigen::RowMajor>  Vector;

    /// Default constructor
    Scaleprior();

    /// Copy constructor
    explicit Scaleprior(const Scaleprior & prior);

    /// Default destructor
    virtual ~Scaleprior();

    /// Returns the scale prior parameter
    const Param &  w() const;
    Param &  getW();

    /// Returns the parameter lower bound
    const Param & lowerBound() const;
    Param & getLowerBound();

    /// Returns the gradient
    const dParam & gradient() const;
    dParam & getGradient();

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
    Param  w_;
    Param  lowerBound_;

    dParam  gradient_;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

}; // class Scaleprior



/// The Deformation class represents a 2d quadratic deformation (dx^2 dx dy^2 dy).
class Deformation : public ParamUtil
{
public:
    /// Type of a 2d quadratic deformation
    typedef Eigen::Matrix<Scalar, 4, 1, Eigen::ColMajor, 4, 1> Param;

#ifndef RGM_USE_DOUBLE
    typedef Eigen::Matrix<double, 4, 1, Eigen::ColMajor, 4, 1> dParam;
#else
    typedef Param   dParam;
#endif

    /// Bounded shift in DT
    static const int BoundedShiftInDT = 4;

    /// Default constructor
    Deformation();

    /// Copy constructor
    explicit Deformation(const Deformation & def);

    /// Constructs a deformation with given arguments
    explicit Deformation(Scalar dx, Scalar dy);

    /// Default destructor
    virtual ~Deformation();

    /// Returns the scale prior parameter
    const Param &  w() const;
    Param &  getW();

    /// Returns the parameter lower bound
    const Param & lowerBound() const;
    Param & getLowerBound();

    /// Returns the gradient
    const dParam & gradient() const;
    dParam & getGradient();

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
    Param  w_;
    Param  lowerBound_;

    dParam  gradient_;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

}; // class Deformation



/// The ParseInfo class
class ParseInfo : public Rectangle_<Scalar>
{
public:
    /// Default constructor
    ParseInfo();

    /// Constructs a parse info with given arguments
    explicit ParseInfo(int c, int l, int x, int y, int ds, int dx, int dy, Scalar score,
              const Rectangle_<Scalar> & bbox, Scalar loss = 0);

    /// Copy constructor
    explicit ParseInfo(const ParseInfo & info);

    /// clip the bbox
    bool clipBbox(int wd, int ht);

    int c_;
    int l_;
    int x_;
    int y_;
    int ds_;
    int dx_;
    int dy_;

    Scalar score_;
    Scalar loss_;

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

}; // class ParseNode


} // namespace RGM


#endif // RGM_PARAMETERS_HPP_




