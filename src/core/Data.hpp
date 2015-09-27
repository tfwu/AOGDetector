#ifndef RGM_DATA_HPP_
#define RGM_DATA_HPP_

#include "Rectangle.hpp"
#include "XMLReader.hpp"
#include "UtilLog.hpp"


namespace RGM
{
/// The Object class represents an object in a Scene. It stores all the information present
/// in between <object> tags in a Pascal VOC .xml annotation file, although the bounding box is
/// represented slightly differently (top left coordinates of (0, 0) instead of (1, 1)).
/// It is adapted from FFLDv2 (the Fast Fourier Linear Detector version 2)
/// Copyright (c) 2013 Idiap Research Institute, <http://www.idiap.ch/>
/// Written by Charles Dubout <charles.dubout@idiap.ch>
class Object
{
public:
    /// The possible object views.
    enum Pose {
        FRONTAL, LEFT, REAR, RIGHT, UNSPECIFIED
    };

    /// Constructs an empty object. An empty object has name 'unknown', pose 'unspecified', and all
    /// other parameters set to their default values.
    Object();

    /// Constructs an object from a name, a pose, annotation flags and a bounding box.
    /// @param[in] pose View of the object.
    /// @param[in] truncated Whether the object is annotated as being truncated.
    /// @param[in] difficult Whether the object is annotated as being difficult.
    /// @param[in] bndbox Bounding box of the object.
    Object(Pose pose, bool truncated, bool difficult, bool flip, Rectangle2i & bndbox, int dataId);

    /// Returns the pose (view) of the object.
    Pose   pose() const;
    Pose & getPose();

    /// Returns whether the object is annotated as being truncated.
    bool   truncated() const;
    bool & getTruncated();

    /// Returns whether the object is annotated as being difficult.
    bool    difficult() const;
    bool &  getDifficult();

    /// Returns the bounding box of the object.
    const Rectangle2i & bndbox() const;
    Rectangle2i & getBndbox();

    /// Returns the data ID
    int   dataId() const;
    int & getDataId();

    /// Returns whether the object is empty.
    /// An empty object has name 'unknown', pose 'unspecified', and all other parameters set to
    /// their default values.
    bool empty() const;

private:
    Pose pose_;
    bool truncated_;
    bool difficult_;
    Rectangle2i bndbox_;
    int dataId_;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

}; // class Object


/// Positive data represents a positive image with a list of foreground objects
class PosData
{
public:
    /// Default constructor
    PosData();

    /// Returns the image file
    const std::string & imgFile() const;
    std::string & getImgFile();

    /// Returns the object category name
    const std::string & objName() const;
    std::string & getObjName();

    /// Returns the objects
    const std::vector<Object> & objects() const;
    std::vector<Object> & getObjects();

    /// Returns the left-right flipping status
    bool   flip() const;
    bool & getFlip();

    /// Returns the maximum area of object bounding boxes
    int   maxArea() const;
    int & getMaxArea();

    /// Crops the objects and returns the new bboxes
    cv::Mat getCropImg(std::vector<Rectangle2i> & cropBoxes) const;
    cv::Mat getCropImg(int idx, Rectangle2i & cropBox) const;

    /// Crops the object image patch
    /// @note @p wd and @p ht are in the unit of @p cellSize
    std::pair<cv::Mat, cv::Mat> cropObjectPatch(int objIdx, int cellSize, int wd, int ht, bool padding=true) const;

    /// Get a set of data from a given set of (imgIdx, bboxIdx)'s
    static void getPosData(const std::vector<PosData> & pos, std::vector<std::pair<int, int> > & idx,
                           int minSz, std::vector<PosData> & usedPos);

    /// Get resized object patches w.r.t. a given protoarea
    void getResizedObjectPatch(std::vector<cv::Mat> &imgs, float protoarea, int cellSize, int factor) const;

private:
    std::string              imgFile_;
    std::string              objName_;
    std::vector<Object>      objects_;
    bool                     flip_;
    int                      maxArea_;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);


}; // Positive data


/// Negative data
class NegData
{
public:
    /// Default constructor
    NegData();

    /// Get an image
    cv::Mat getImg() const;

    /// Return image file
    const std::string & imgFile() const;
    std::string & getImgFile();

    /// Return status of left-right flipping
    bool   flip() const;
    bool & getFlip();

    /// Return the data id
    int   dataId() const;
    int & getDataId();

private:
    std::string imgFile_;
    bool        flip_;
    int         dataId_;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

}; // class NegData


} // namespce RGM

#endif // RGM_DATA_HPP_


