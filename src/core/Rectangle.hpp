// This file is adapted from FFLDv2 (the Fast Fourier Linear Detector version 2)
// Copyright (c) 2013 Idiap Research Institute, <http://www.idiap.ch/>
// Written by Charles Dubout <charles.dubout@idiap.ch>

#ifndef RGM_RECTANGLE_HPP_
#define RGM_RECTANGLE_HPP_

#include <stdio.h>
#include <vector>
#include <algorithm>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "opencv2/opencv.hpp"

#include "Common.hpp"

namespace RGM
{
/// The Rectangle_ class defines a rectangle in the plane using T precision. If the coordinates
/// of the top left corner of the rectangle are (x, y), the coordinates of the bottom right corner
/// are (x + width - 1, y + height - 1), where width and height are the dimensions of the rectangle.
/// The corners are thus understood as the extremal points still inside the rectangle.
template <typename T>
class Rectangle_
{
public:
    /// Constructs an empty rectangle. An empty rectangle has no area.
    Rectangle_();

    /// Constructs a rectangle with the given @p width and @p height.
    Rectangle_(T width, T height);

    /// Constructs a rectangle with coordinates (@p x, @p y) and the given @p width and @p height.
    Rectangle_(T x, T y, T width, T height);

    /// Constructs a rectangle by copying from @p rect
    Rectangle_(const Rectangle_<T> & rect);

    /// Constructs a rectangle with cv::Rect_ @p rect
    Rectangle_(const cv::Rect_<T> & rect);

    /// Returns the x-coordinate of the rectangle.
    T x() const;

    /// Sets the x coordinate of the rectangle to @p x.
    void setX(T x);

    /// Returns the y-coordinate of the rectangle.
    T y() const;

    /// Sets the y coordinate of the rectangle to @p y.
    void setY(T y);

    /// Returns the width of the rectangle.
    T width() const;

    /// Sets the height of the rectangle to the given @p width.
    void setWidth(T width);

    /// Returns the height of the rectangle.
    T height() const;

    /// Sets the height of the rectangle to the given @p height.
    void setHeight(T height);

    /// Returns the left side of the rectangle.
    /// @note Equivalent to x().
    T left() const;

    /// Sets the left side of the rectangle to @p left.
    /// @note The right side of the rectangle is not modified.
    void setLeft(T left);

    /// Returns the top side of the rectangle.
    /// @note Equivalent to y().
    T top() const;

    /// Sets the top side of the rectangle to @p top.
    /// @note The bottom side of the rectangle is not modified.
    void setTop(T top);

    /// Returns the right side of the rectangle.
    /// @note Equivalent to x() + width() - 1.
    T right() const;

    /// Sets the right side of the rectangle to @p right.
    /// @note The left side of the rectangle is not modified.
    void setRight(T right);

    /// Returns the bottom side of the rectangle.
    /// @note Equivalent to y() + height() - 1.
    T bottom() const;

    /// Sets the bottom side of the rectangle to @p bottom.
    /// @note The top side of the rectangle is not modified.
    void setBottom(T bottom);

    /// Returns the center
    T xcenter() const;
    T ycenter() const;

    /// Sets the width and height to the maximum
    void setMaxWidth(T wd);
    void setMaxHeight(T ht);
    void setMax(const Rectangle_<T> & rect);

    /// Sets the width and height to the minimum
    void setMinWidth(T wd);
    void setMinHeight(T ht);
    void setMin(const Rectangle_<T> & rect);

    /// Returns whether the rectangle is empty. An empty rectangle has no area.
    bool empty() const;

    /// Returns the area of the rectangle.
    /// @note Equivalent to max(width(), 0) * max(height(), 0).
    T area() const;

    /// Return if it is the same as rect
    bool operator== (const Rectangle_<T> & rect) const;

    /// Returns if it is the same type as rect
    /// @note The type of a rectangle is defined by @p width and @p height
    bool isSameType(const Rectangle_<T> & rect) const;

    /// Returns the partitioned rectangles with @p rect
    /// @note If it does not overlap with @p rect Returns both of them directly
    /// @note If it overlaps with @p rect Returns the partitioned rectangles
    std::vector<Rectangle_<T> > partition(Rectangle_<T> &rect);

    /// Converts to cv::Rect
    cv::Rect_<T> cvRect() const;

private:
    T x_;
    T y_;
    T width_;
    T height_;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

}; // class Rectangle_

/// Typedef
typedef Rectangle_<int> Rectangle2i;


/// Returns all the overlapped rectangles between any pair of rectangles in @p rects1 and @p rects2
std::vector<Rectangle2i> getOverlappedRectangles(std::vector<Rectangle2i> &rects1, std::vector<Rectangle2i> &rects2);


/// The Detection class
template <typename T>
class Detection_ : public Rectangle_<T>
{
public:
    /// Constructs an empty Detection
    Detection_();

    /// Constructs a Detection with @p score and @p bndbox
    Detection_(Scalar score, const Rectangle_<T> & bndbox);

    /// Constructs a Detection with @p l, @p x, @p y and @p score
    Detection_(int l, int x, int y, Scalar score);

    /// Constructs a Detection with @p c, @p score and @p bndbox
    Detection_(int c, Scalar score, const Rectangle_<T> & bndbox);

    /// Constructs a Detection with full specification except for pt idx
    Detection_(int c, int l, int x, int y, Scalar score, const Rectangle_<T> & bndbox);

    /// Constructs a Detection with full specification
    Detection_(int c, int l, int x, int y, Scalar score, const Rectangle_<T> & bndbox, int ptIdx, int ptNodeIdx);

    /// Copy constructor
    Detection_(const Detection_<T> & detection);

    /// Compares the scores in decreasing order
    bool operator<(const Detection_<T> & detection) const;

    /// clip the bbox
    bool clipBbox(int wd, int ht);

    /// show
    void show(cv::Mat img, bool display);

    int    c_;
    int    l_;
    int    x_;
    int    y_;
    Scalar score_;

    int ptIdx_;
    int ptNodeIdx_; // usually for single object And PtNode

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);


}; // class Detection


/// Functor used to test for the intersection of two rectangles
/// according to the Pascal criterion (area of intersection over area of union).
template <typename T>
class Intersector_
{
public:
    /// Constructor.
    /// @param[in] reference The reference rectangle.
    /// @param[in] threshold The threshold of the criterion.
    /// @param[in] dividedByUnion Use Felzenszwalb's criterion instead (area of intersection over area
    /// of second rectangle). Useful to remove small detections inside bigger ones.
    Intersector_(const Rectangle_<T> & reference, float threshold = 0.5, bool dividedByUnion = false);

    /// Tests for the intersection between a given rectangle and the reference.
    /// @param[in] rect The rectangle to intersect with the reference.
    /// @param[out] score The score of the intersection.
    bool operator()(const Rectangle_<T> & rect, float * score = 0) const;

private:
    const Rectangle_<T> & reference_;
    float threshold_;
    bool   dividedByUnion_;

}; // class Intersector_


} // namespace RGM


#endif // RGM_RECTANGLE_HPP_



