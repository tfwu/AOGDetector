#include <algorithm>
#include <iostream>

#include "Rectangle.hpp"

namespace RGM
{

using namespace std;

// ------- Rectangle_ ------

template <typename T>
Rectangle_<T>::Rectangle_() : x_(0), y_(0), width_(0), height_(0)
{
}

template <typename T>
Rectangle_<T>::Rectangle_(T width, T height) : x_(0), y_(0), width_(width), height_(height)
{
}

template <typename T>
Rectangle_<T>::Rectangle_(T x, T y, T width, T height) : x_(x), y_(y), width_(width), height_(height)
{
}

template <typename T>
Rectangle_<T>::Rectangle_(const Rectangle_<T> & rect)
{
    setX(rect.x());
    setY(rect.y());
    setWidth(rect.width());
    setHeight(rect.height());
}

template <typename T>
Rectangle_<T>::Rectangle_(const cv::Rect_<T> & rect)
{
    setX(rect.x);
    setY(rect.y);
    setWidth(rect.width);
    setHeight(rect.height);
}

template <typename T>
T Rectangle_<T>::x() const
{
    return x_;
}

template <typename T>
void Rectangle_<T>::setX(T x)
{
    x_ = x;
}

template <typename T>
T Rectangle_<T>::y() const
{
    return y_;
}

template <typename T>
void Rectangle_<T>::setY(T y)
{
    y_ = y;
}

template <typename T>
T Rectangle_<T>::width() const
{
    return width_;
}

template <typename T>
void Rectangle_<T>::setWidth(T width)
{
    width_ = width;
}

template <typename T>
T Rectangle_<T>::height() const
{
    return height_;
}

template <typename T>
void Rectangle_<T>::setHeight(T heigth)
{
    height_ = heigth;
}

template <typename T>
T Rectangle_<T>::left() const
{
    return x();
}

template <typename T>
void Rectangle_<T>::setLeft(T left)
{
    setWidth(right() - left + 1);
    setX(left);
}

template <typename T>
T Rectangle_<T>::top() const
{
    return y();
}

template <typename T>
void Rectangle_<T>::setTop(T top)
{
    setHeight(bottom() - top + 1);
    setY(top);
}


template <typename T>
T Rectangle_<T>::right() const
{
    return x() + width() - 1;
}

template <typename T>
void Rectangle_<T>::setRight(T right)
{
    setWidth(right - left() + 1);
}

template <typename T>
T Rectangle_<T>::bottom() const
{
    return y() + height() - 1;
}

template <typename T>
void Rectangle_<T>::setBottom(T bottom)
{
    setHeight(bottom - top() + 1);
}

template <typename T>
T Rectangle_<T>::xcenter() const
{
    return x() + width() / 2;
}

template <typename T>
T Rectangle_<T>::ycenter() const
{
    return y() + height() / 2;
}

template <typename T>
void Rectangle_<T>::setMaxWidth(T wd)
{
    setWidth(std::max(wd, width()));
}

template <typename T>
void Rectangle_<T>::setMaxHeight(T ht)
{
    setHeight(std::max(ht, height()));
}

template <typename T>
void Rectangle_<T>::setMax(const Rectangle_<T> & rect)
{
    setMaxWidth(rect.width());
    setMaxHeight(rect.height());
}

template <typename T>
void Rectangle_<T>::setMinWidth(T wd)
{
    setWidth(std::min(wd, width()));
}

template <typename T>
void Rectangle_<T>::setMinHeight(T ht)
{
    setHeight(std::min(ht, height()));
}

template <typename T>
void Rectangle_<T>::setMin(const Rectangle_<T> & rect)
{
    setMinWidth(rect.width());
    setMinHeight(rect.height());
}

template <typename T>
bool Rectangle_<T>::empty() const
{
    return (width() <= 0) || (height() <= 0) || x()==std::numeric_limits<T>::quiet_NaN() ||
           y()==std::numeric_limits<T>::quiet_NaN() ||
           width()==std::numeric_limits<T>::quiet_NaN() ||
           height()==std::numeric_limits<T>::quiet_NaN();
}

template <typename T>
T Rectangle_<T>::area() const
{
    return max<T>(width(), 0) * max<T>(height(), 0);
}

template <typename T>
bool Rectangle_<T>::operator==(const Rectangle_<T> &rect) const
{
    return (x()==rect.x() && y()==rect.y() &&
            width()==rect.width() && height()==rect.height());
}

template <typename T>
bool Rectangle_<T>::isSameType(const Rectangle_<T> &rect) const
{
    return (width()==rect.width() && height()==rect.height());
}

template <typename T>
std::vector<Rectangle_<T> > Rectangle_<T>::partition(Rectangle_<T> &rect)
{
    std::vector<Rectangle_<T> > bb;

    // Not overlap with rect
    T x1 = std::max<T>(left(), rect.left());
    T x2 = std::min<T>(right(), rect.right());
    if (x1>x2) {
        bb.push_back(*this);
        bb.push_back(rect);
        return bb;
    }

    T y1 = std::max<T>(top(), rect.top());
    T y2 = std::min<T>(bottom(), rect.bottom());
    if (y1 > y2) {
        bb.push_back(*this);
        bb.push_back(rect);
        return bb;
    }

    // Overlap with rect
    Rectangle_<T> overlap(x1, y1, x2-x1+1, y2-y1+1);

    bb.push_back(overlap);
    bb.push_back(overlap);

    if (left()==rect.left() && right()==rect.right()) {
        // Hor cut
        if (top() < rect.top()) {
            bb.push_back( Rectangle_<T>(x1, top(),    x2-x1+1, y1-top()+1) );
            bb.push_back( Rectangle_<T>(x1, bottom(), x2-x1+1, rect.bottom()-bottom()+1) );
        } else {
            bb.push_back( Rectangle_<T>(x1, rect.top(),    x2-x1+1, y1-rect.top()+1) );
            bb.push_back( Rectangle_<T>(x1, rect.bottom(), x2-x1+1, bottom()-rect.bottom()+1) );
        }
    } else {
        // Ver cut
        if (left() < rect.left() ) {
            bb.push_back( Rectangle_<T>(left(),  y1, x1-left()+1,            y2-y1+1) );
            bb.push_back( Rectangle_<T>(right(), y1, rect.right()-right()+1, y2-y1+1) );
        } else {
            bb.push_back( Rectangle_<T>(rect.left(),  y1, x1-rect.left()+1,       y2-y1+1) );
            bb.push_back( Rectangle_<T>(rect.right(), y1, right()-rect.right()+1, y2-y1+1) );
        }
    }

    return bb;
}

template <typename T>
cv::Rect_<T> Rectangle_<T>::cvRect() const
{
    return cv::Rect_<T>(x(), y(), width(), height());
}

template <typename T>
template <class Archive>
void Rectangle_<T>::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(x_);
    ar & BOOST_SERIALIZATION_NVP(y_);
    ar & BOOST_SERIALIZATION_NVP(width_);
    ar & BOOST_SERIALIZATION_NVP(height_);
}

/// Specification
template class Rectangle_<int>;
template class Rectangle_<Scalar>;

INSTANTIATE_BOOST_SERIALIZATION(Rectangle_<int>);
INSTANTIATE_BOOST_SERIALIZATION(Rectangle_<Scalar>);



std::vector<Rectangle2i> getOverlappedRectangles(std::vector<Rectangle2i> &rects1, std::vector<Rectangle2i> &rects2)
{
    std::vector<Rectangle2i> bb;

    for ( int i=0; i<rects2.size(); ++i ) {
        Rectangle2i &rect2( rects2[i] );
        for ( int j=0; j<rects1.size(); ++j ) {
            Rectangle2i &rect1( rects1[j] );

            int x1 = std::max<int>(rect1.left(), rect2.left());
            int y1 = std::max<int>(rect1.top(), rect2.top());
            int x2 = std::min<int>(rect1.right(), rect2.right());
            int y2 = std::min<int>(rect1.bottom(), rect2.bottom());

            if (x1>x2 || y1>y2) {
                continue;
            }

            bb.push_back( Rectangle2i(x1, y1, x2-x1+1, y2-y1+1) );
        }
    }

    return bb;
}


// ------- Detection --------

template<typename T>
Detection_<T>::Detection_() :
    c_(-1), l_(-1), x_(-1), y_(-1), score_(-10.0F), ptIdx_(-1), ptNodeIdx_(-1)
{
}

template<typename T>
Detection_<T>::Detection_(Scalar score, const Rectangle_<T> & bndbox) :
    c_(-1), l_(-1), x_(-1), y_(-1), score_(score), Rectangle_<T>(bndbox), ptIdx_(-1), ptNodeIdx_(-1)
{
}

template<typename T>
Detection_<T>::Detection_(int l, int x, int y, Scalar score) :
    c_(-1), l_(l), x_(x), y_(y), score_(score), ptIdx_(-1), ptNodeIdx_(-1)
{
}

template<typename T>
Detection_<T>::Detection_(int c, Scalar score, const Rectangle_<T> & bndbox) :
    c_(c), l_(-1), x_(-1), y_(-1), score_(score), Rectangle_<T>(bndbox), ptIdx_(-1), ptNodeIdx_(-1)
{
}

template<typename T>
Detection_<T>::Detection_(const Detection_<T> & detection) :
    c_(detection.c_), l_(detection.l_), x_(detection.x_), y_(detection.y_), score_(detection.score_),
    Rectangle_<T>(detection.x(), detection.y(), detection.width(), detection.height()),
    ptIdx_(detection.ptIdx_), ptNodeIdx_(detection.ptNodeIdx_)
{
}

template<typename T>
Detection_<T>::Detection_(int c, int l, int x, int y, Scalar score, const Rectangle_<T> & bndbox) :
    c_(c), l_(l), x_(x), y_(y), score_(score), Rectangle_<T>(bndbox), ptIdx_(-1), ptNodeIdx_(-1)
{
}

template<typename T>
Detection_<T>::Detection_(int c, int l, int x, int y, Scalar score, const Rectangle_<T> & bndbox, int ptIdx, int ptNodeIdx) :
    c_(c), l_(l), x_(x), y_(y), score_(score), Rectangle_<T>(bndbox), ptIdx_(ptIdx), ptNodeIdx_(ptNodeIdx)
{
}

template<typename T>
bool Detection_<T>::operator<(const Detection_<T> & detection) const
{
    return score_ > detection.score_; // for decreasing sort
}

template<typename T>
bool Detection_<T>::clipBbox(int wd, int ht)
{
    T x1                = std::max<T>(0, this->x());
    T y1                = std::max<T>(0, this->y());
    T x2                = std::min<T>(wd-1, this->right());
    T y2                = std::min<T>(ht-1, this->bottom());
    T width             = x2 - x1 + 1;
    T height            = y2 - y1 + 1;

    if (width <=0 || height <=0) {
        return false;
    }

    this->setX(x1);
    this->setY(y1);
    this->setWidth(width);
    this->setHeight(height);

    return true;
}

template<typename T>
void Detection_<T>::show(cv::Mat img, bool display)
{
    cv::rectangle(img, this->cvRect(), cv::Scalar::all(255), 5);
    cv::rectangle(img, this->cvRect(), cv::Scalar(0, 0, 255), 3);

    if ( display ) {
        cv::String winName("AOGDetection");
        cv::imshow(winName, img);
        cv::waitKey(0);
    }
}

template<typename T>
template<class Archive>
void Detection_<T>::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Rectangle_<T>);
    ar & BOOST_SERIALIZATION_NVP(c_);
    ar & BOOST_SERIALIZATION_NVP(l_);
    ar & BOOST_SERIALIZATION_NVP(x_);
    ar & BOOST_SERIALIZATION_NVP(y_);
    ar & BOOST_SERIALIZATION_NVP(score_);
    ar & BOOST_SERIALIZATION_NVP(ptIdx_);
    ar & BOOST_SERIALIZATION_NVP(ptNodeIdx_);
}

/// Specification
template class Detection_<int>;
template class Detection_<Scalar>;

INSTANTIATE_BOOST_SERIALIZATION(Detection_<int>);
INSTANTIATE_BOOST_SERIALIZATION(Detection_<Scalar>);


// ------- Intersector ------

template <typename T>
Intersector_<T>::Intersector_(const Rectangle_<T> & reference, float threshold, bool dividedByUnion) :
    reference_(reference), threshold_(threshold), dividedByUnion_(dividedByUnion)
{
}

template <typename T>
bool Intersector_<T>::operator()(const Rectangle_<T> & rect, float * score) const
{
    if (score) {
        *score = 0.0;
    }

    const int left = std::max<int>(reference_.left(), rect.left());
    const int right = std::min<int>(reference_.right(), rect.right());

    if (right < left) {
        return false;
    }

    const int top = std::max<int>(reference_.top(), rect.top());
    const int bottom = std::min<int>(reference_.bottom(), rect.bottom());

    if (bottom < top) {
        return false;
    }

    const int intersectionArea = (right - left + 1) * (bottom - top + 1);
    const int rectArea = rect.area();

    if (dividedByUnion_) {
        const int referenceArea = reference_.area();
        const int unionArea = referenceArea + rectArea - intersectionArea;

        if (intersectionArea >= unionArea * threshold_) {
            if (score) {
                *score = static_cast<float>(intersectionArea) / unionArea;
            }

            return true;
        }
    } else {
        if (intersectionArea >= rectArea * threshold_) {
            if (score) {
                *score = static_cast<float>(intersectionArea) / rectArea;
            }

            return true;
        }
    }

    return false;
}

/// Specification
template class Intersector_<int>;
template class Intersector_<Scalar>;

} // namespace RGM
