#include <boost/serialization/base_object.hpp>

#include "Parameters.hpp"
#include "UtilSerialization.hpp"

namespace RGM
{
// ------ ParamUtil ------

ParamUtil::ParamUtil() :
    learningStatus_(0), regularizationCost_(0)
{
}

ParamUtil::~ParamUtil()
{
}

Scalar ParamUtil::learningStatus() const
{
    return learningStatus_;
}

Scalar & ParamUtil::getLearningStatus()
{
    return learningStatus_;
}

Scalar ParamUtil::regularizationCost() const
{
    return regularizationCost_;
}

Scalar & ParamUtil::getRegularizationCost()
{
    return regularizationCost_;
}

template<class Archive>
void ParamUtil::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(learningStatus_);
    ar & BOOST_SERIALIZATION_NVP(regularizationCost_);
}

INSTANTIATE_BOOST_SERIALIZATION(ParamUtil);



// ------ Appearance ------

Appearance::Appearance()
{
}

Appearance::Appearance(const Appearance & app)
{
    getW()                  = app.w();
    getLowerBound()         = app.lowerBound();
    getLearningStatus()     = app.learningStatus();
    getRegularizationCost() = app.regularizationCost();
}

Appearance::~Appearance()
{
}

void Appearance::init(int wd, int ht)
{
    Scalar Inf              = std::numeric_limits<Scalar>::infinity();

    getW()                  = Param::Constant(ht, wd, Cell::Zero());
    getLowerBound()         = Param::Constant(ht, wd, Cell::Constant(-Inf));

    getRegularizationCost() = 1.0F;
    getLearningStatus()     = 1.0F;
}

const Appearance::Param & Appearance::w() const
{
    return w_;
}

Appearance::Param & Appearance::getW()
{
    return w_;
}

const Appearance::Param & Appearance::lowerBound() const
{
    return lowerBound_;
}

Appearance::Param & Appearance::getLowerBound()
{
    return lowerBound_;
}

const Appearance::dParam & Appearance::gradient() const
{
    return gradient_;
}

Appearance::dParam & Appearance::getGradient()
{
    return gradient_;
}

template<class Archive>
void Appearance::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ParamUtil);
    ar & BOOST_SERIALIZATION_NVP(w_);
    ar & BOOST_SERIALIZATION_NVP(lowerBound_);
}

INSTANTIATE_BOOST_SERIALIZATION(Appearance);



// ------ Offset ------

Offset::Offset() :
    w_(0), lowerBound_(0), gradient_(0)
{
}

Offset::Offset(const Offset & off)
{
    getW()                  = off.w();
    getLowerBound()         = off.lowerBound();
    getLearningStatus()     = off.learningStatus();
    getRegularizationCost() = off.regularizationCost();
}

Offset::~Offset()
{
}

Scalar Offset::w() const
{
    return w_;
}

Scalar & Offset::getW()
{
    return w_;
}

Scalar Offset::lowerBound() const
{
    return lowerBound_;
}

Scalar & Offset::getLowerBound()
{
    return lowerBound_;
}

const double & Offset::gradient() const
{
    return gradient_;
}

double & Offset::getGradient()
{
    return gradient_;
}

template<class Archive>
void Offset::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ParamUtil);
    ar & BOOST_SERIALIZATION_NVP(w_);
    ar & BOOST_SERIALIZATION_NVP(lowerBound_);
}

INSTANTIATE_BOOST_SERIALIZATION(Offset);



// ------ Scaleprior ------

Scaleprior::Scaleprior() :
    w_(Param::Zero()), lowerBound_(Param::Zero())
{
}

Scaleprior::Scaleprior(const Scaleprior & prior)
{
    getW()                  = prior.w();
    getLowerBound()         = prior.lowerBound();
    getLearningStatus()     = prior.learningStatus();
    getRegularizationCost() = prior.regularizationCost();
}

Scaleprior::~Scaleprior()
{
}

const Scaleprior::Param & Scaleprior::w() const
{
    return w_;
}

Scaleprior::Param & Scaleprior::getW()
{
    return w_;
}

const Scaleprior::Param & Scaleprior::lowerBound() const
{
    return lowerBound_;
}

Scaleprior::Param & Scaleprior::getLowerBound()
{
    return lowerBound_;
}

const Scaleprior::dParam & Scaleprior::gradient() const
{
    return gradient_;
}

Scaleprior::dParam & Scaleprior::getGradient()
{
    return gradient_;
}

template<class Archive>
void Scaleprior::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ParamUtil);
    ar & BOOST_SERIALIZATION_NVP(w_);
    ar & BOOST_SERIALIZATION_NVP(lowerBound_);
}

INSTANTIATE_BOOST_SERIALIZATION(Scaleprior);



// ------ Deformation ------

Deformation::Deformation() :
    w_(Param::Zero()), lowerBound_(Param::Zero())
{
}

Deformation::Deformation(const Deformation & def)
{
    getW()                  = def.w();
    getLowerBound()         = def.lowerBound();
    getLearningStatus()     = def.learningStatus();
    getRegularizationCost() = def.regularizationCost();
}

Deformation::Deformation(Scalar dx, Scalar dy) :
    lowerBound_(Param::Zero())
{
    getW() << dx * dx, dx, dy * dy, dy;
}

Deformation::~Deformation()
{
}

const Deformation::Param & Deformation::w() const
{
    return w_;
}

Deformation::Param & Deformation::getW()
{
    return w_;
}

const Deformation::Param & Deformation::lowerBound() const
{
    return lowerBound_;
}

Deformation::Param & Deformation::getLowerBound()
{
    return lowerBound_;
}

const Deformation::dParam & Deformation::gradient() const
{
    return gradient_;
}

Deformation::dParam & Deformation::getGradient()
{
    return gradient_;
}

template<class Archive>
void Deformation::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ParamUtil);
    ar & BOOST_SERIALIZATION_NVP(w_);
    ar & BOOST_SERIALIZATION_NVP(lowerBound_);
}

INSTANTIATE_BOOST_SERIALIZATION(Deformation);



// ------ ParseInfo ------

ParseInfo::ParseInfo() :
    c_(-1), l_(-1), x_(-1), y_(-1), ds_(-1), dx_(-1), dy_(-1), score_(-10.0F), loss_(0.0F)
{
}

ParseInfo::ParseInfo(int c, int l, int x, int y, int ds, int dx, int dy, Scalar score,
                     const Rectangle_<Scalar> & bbox, Scalar loss) :
    c_(c), l_(l), x_(x), y_(y), ds_(ds), dx_(dx), dy_(dy), score_(score), loss_(loss),
    Rectangle_<Scalar>(bbox)
{
}

ParseInfo::ParseInfo(const ParseInfo & info) :
    c_(info.c_), l_(info.l_), x_(info.x_), y_(info.y_),
    ds_(info.ds_), dx_(info.dx_), dy_(info.dy_), score_(info.score_), loss_(info.loss_),
    Rectangle_<Scalar>(info.x(), info.y(), info.width(), info.height())
{
}

bool ParseInfo::clipBbox(int wd, int ht)
{
    float x1                = std::max<float>(0, x());
    float y1                = std::max<float>(0, y());
    float x2                = std::min<float>(wd-1, right());
    float y2                = std::min<float>(ht-1, bottom());
    float width             = x2 - x1 + 1;
    float height            = y2 - y1 + 1;

    if (width <=0 || height <=0) {
        return false;
    }

    setX(x1);
    setY(y1);
    setWidth(width);
    setHeight(height);

    return true;
}

template<class Archive>
void ParseInfo::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Rectangle_<Scalar>);
    ar & BOOST_SERIALIZATION_NVP(c_);
    ar & BOOST_SERIALIZATION_NVP(l_);
    ar & BOOST_SERIALIZATION_NVP(x_);
    ar & BOOST_SERIALIZATION_NVP(y_);
    ar & BOOST_SERIALIZATION_NVP(ds_);
    ar & BOOST_SERIALIZATION_NVP(dx_);
    ar & BOOST_SERIALIZATION_NVP(dy_);
    ar & BOOST_SERIALIZATION_NVP(score_);
    ar & BOOST_SERIALIZATION_NVP(loss_);
}

INSTANTIATE_BOOST_SERIALIZATION(ParseInfo);



} // namespace RGM



