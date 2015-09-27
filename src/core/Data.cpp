#include "Data.hpp"
#include "UtilOpencv.hpp"


namespace RGM
{

using namespace std;


// ----- Object -----

Object::Object() :
    pose_(UNSPECIFIED), truncated_(false), difficult_(false), dataId_(-1)
{
}

Object::Object(Pose pose, bool truncated, bool difficult, bool flip, Rectangle2i & bndbox, int dataId) :
    pose_(pose), truncated_(truncated), difficult_(difficult), bndbox_(bndbox), dataId_(dataId)
{
}

Object::Pose Object::pose() const
{
    return pose_;
}

Object::Pose & Object::getPose()
{
    return pose_;
}

bool Object::truncated() const
{
    return truncated_;
}

bool & Object::getTruncated()
{
    return truncated_;
}

bool Object::difficult() const
{
    return difficult_;
}

bool & Object::getDifficult()
{
    return difficult_;
}

const Rectangle2i & Object::bndbox() const
{
    return bndbox_;
}

Rectangle2i & Object::getBndbox()
{
    return bndbox_;
}

int Object::dataId() const
{
    return dataId_;
}

int & Object::getDataId()
{
    return dataId_;
}

bool Object::empty() const
{
    return (pose() == UNSPECIFIED) && !truncated() && !difficult() && bndbox().empty() && (dataId()==-1);
}

template<class Archive>
void Object::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(pose_);
    ar & BOOST_SERIALIZATION_NVP(truncated_);
    ar & BOOST_SERIALIZATION_NVP(difficult_);
    ar & BOOST_SERIALIZATION_NVP(bndbox_);
    ar & BOOST_SERIALIZATION_NVP(dataId_);
}

INSTANTIATE_BOOST_SERIALIZATION(Object);



// ----- PosData -----

PosData::PosData() :
    imgFile_("UNKNOWN"), objName_("UNKNOWN"), flip_(false), maxArea_(0)
{
}

const std::string & PosData::imgFile() const
{
    return imgFile_;
}

std::string & PosData::getImgFile()
{
    return imgFile_;
}
const std::string & PosData::objName() const
{
    return objName_;
}

std::string & PosData::getObjName()
{
    return objName_;
}

const std::vector<Object> & PosData::objects() const
{
    return objects_;
}

std::vector<Object> & PosData::getObjects()
{
    return objects_;
}

bool    PosData::flip() const
{
    return flip_;
}

bool &  PosData::getFlip()
{
    return flip_;
}

int   PosData::maxArea() const
{
    return maxArea_;
}

int & PosData::getMaxArea()
{
    return maxArea_;
}

std::pair<cv::Mat, cv::Mat> PosData::cropObjectPatch(int objIdx, int cellSize, int wd, int ht, bool padding) const
{
    DEFINE_RGM_LOGGER;

    cv::Mat img;

    if (objIdx < 0 || objIdx >= objects().size() ) {
        RGM_LOG(error, "Wrong object index");
        return std::make_pair(img, img);
    }

    int pixelWd = cellSize * wd;
    int pixelHt = cellSize * ht;

    int cropHt =  padding ? (cellSize * (ht+2)) : pixelHt;
    int cropWd =  padding ? (cellSize * (wd+2)) : pixelWd;
    cv::Size cropSz(cropWd, cropHt);

#pragma omp critical
    {
        img = cv::imread(imgFile(), cv::IMREAD_COLOR);
    }

    if (flip()) {
        cv::flip(img, img, 1);
    }

    const Object & obj(objects()[objIdx]);

    int bbHt = obj.bndbox().height();
    int bbWd = obj.bndbox().width();

    int pady = padding ? ROUND((float)cellSize * bbHt / pixelHt) : 0;
    int padx = padding ? ROUND((float)cellSize * bbWd / pixelWd) : 0;

    int x1 = obj.bndbox().x()      - padx;
    int x2 = obj.bndbox().right()  + padx;
    int y1 = obj.bndbox().y()      - pady;
    int y2 = obj.bndbox().bottom() + pady;

    cv::Rect roi(x1, y1, x2-x1+1, y2-y1+1);

    cv::Mat roiImg = OpencvUtil::subarray(img, roi, 1);

    cv::resize(roiImg, img, cropSz, 0, 0, cv::INTER_AREA);

    return std::make_pair(img, roiImg);
}

cv::Mat PosData::getCropImg(std::vector<Rectangle2i> & cropBoxes) const
{
    cv::Mat img;
#pragma omp critical
    {
        img = cv::imread(imgFile_, cv::IMREAD_COLOR);
    }

    if ( flip_ ) {
        cv::flip(img, img, 1);
    }

    cropBoxes.clear();
    cropBoxes.push_back(objects()[0].bndbox());

    int x1 = objects()[0].bndbox().x();
    int x2 = objects()[0].bndbox().right();
    int y1 = objects()[0].bndbox().y();
    int y2 = objects()[0].bndbox().bottom();
    int ht = objects()[0].bndbox().height();
    int wd = objects()[0].bndbox().width();
    for ( int i=1; i<objects().size(); ++i ) {
        x1 = std::min<int>(x1, objects()[i].bndbox().x());
        x2 = std::max<int>(x2, objects()[i].bndbox().right());
        y1 = std::min<int>(y1, objects()[i].bndbox().y());
        y2 = std::max<int>(y2, objects()[i].bndbox().bottom());

        ht = std::max<int>(ht, objects()[i].bndbox().height());
        wd = std::max<int>(wd, objects()[i].bndbox().width());

        cropBoxes.push_back(objects()[i].bndbox());
    }

    const float alpha = 0.5F;
    int padx = ROUND(wd * alpha);
    int pady = ROUND(ht * alpha);

    cv::Rect roi;

    roi.x = std::max<int>(0, x1 - padx);
    roi.y = std::max<int>(0, y1 - pady);
    roi.width = std::min<int>(img.cols, x2 + padx) - roi.x;
    roi.height = std::min<int>(img.rows, y2 + pady) - roi.y;

    cv::Mat cropImg = img(roi).clone();

    for ( int i=0; i<cropBoxes.size(); ++i ) {
        cropBoxes[i].setX(cropBoxes[i].x() - roi.x);
        cropBoxes[i].setY(cropBoxes[i].y() - roi.y);
    }

    return cropImg;
}

cv::Mat PosData::getCropImg(int idx, Rectangle2i & cropBox) const
{
    cv::Mat img;
    if ( idx < 0 || idx >= objects().size() )
        return img;


#pragma omp critical
    {
        img = cv::imread(imgFile_, cv::IMREAD_COLOR);
    }

    if ( flip_ ) {
        cv::flip(img, img, 1);
    }

    int x1 = objects()[idx].bndbox().x();
    int x2 = objects()[idx].bndbox().right();
    int y1 = objects()[idx].bndbox().y();
    int y2 = objects()[idx].bndbox().bottom();
    int ht = objects()[idx].bndbox().height();
    int wd = objects()[idx].bndbox().width();

    const float alpha = 0.5F;
    int padx = ROUND(wd * alpha);
    int pady = ROUND(ht * alpha);

    cv::Rect roi;

    roi.x = std::max<int>(0, x1 - padx);
    roi.y = std::max<int>(0, y1 - pady);
    roi.width = std::min<int>(img.cols, x2 + padx) - roi.x;
    roi.height = std::min<int>(img.rows, y2 + pady) - roi.y;

    cv::Mat cropImg = img(roi).clone();

    cropBox = objects()[idx].bndbox();
    cropBox.setX(cropBox.x() - roi.x);
    cropBox.setY(cropBox.y() - roi.y);

    return cropImg;
}

void PosData::getPosData(const std::vector<PosData> & pos, std::vector<std::pair<int, int> > & idx,
                         int minSz, std::vector<PosData> & usedPos)
{
    int num = idx.size();

    usedPos.clear();

    for ( int i = 0; i < num; ++i ) {
        int dataIdx = idx[i].first;
        int objIdx = idx[i].second;

        PosData d;
        d.getImgFile() = pos[dataIdx].imgFile();
        d.getObjName() = pos[dataIdx].objName();
        d.getFlip() = pos[dataIdx].flip();
        d.getMaxArea() = pos[dataIdx].objects()[objIdx].bndbox().area();

        if ( pos[dataIdx].objects()[objIdx].bndbox().area() < minSz ) {
            continue;
        }

        d.getObjects().push_back(pos[dataIdx].objects()[objIdx]);
        d.getMaxArea() = pos[dataIdx].objects()[objIdx].bndbox().area();

        usedPos.push_back(d);

    } // for i


    //std::sort(idx.begin(), idx.end());

    //int dataIdx = idx[0].first;
    //
    //std::vector<int> boxIdx;
    //boxIdx.push_back(idx[0].second);

    //for ( int i = 1; i < num; ++i ) {
    //	if ( idx[i].first == dataIdx ) {
    //		boxIdx.push_back(idx[i].second);
    //	} else {
    //		PosData d;
    //		d.getImgFile() = pos[dataIdx].imgFile();
    //		d.getObjName() = pos[dataIdx].objName();
    //		d.getFlip() = pos[dataIdx].flip();
    //		d.getMaxArea() = 0;
    //		for ( int j = 0; j < boxIdx.size(); ++j ) {
    //			if ( pos[dataIdx].objects()[boxIdx[j]].bndbox().area() < minSz )
    //				continue;
    //			d.getObjects().push_back(pos[dataIdx].objects()[boxIdx[j]]);
    //			d.getMaxArea() = std::max(d.maxArea(), d.getObjects().back().bndbox().area());
    //		}// for j

    //		if ( d.objects().size() > 0 )
    //			usedPos.push_back(d);

    //		dataIdx = idx[i].first;
    //		boxIdx.clear();
    //		boxIdx.push_back(idx[i].second);
    //	}
    //} // for i
}

void PosData::getResizedObjectPatch(std::vector<cv::Mat> &imgs, float protoarea,
                                    int cellSize, int factor) const
{
    DEFINE_RGM_LOGGER;

    RGM_CHECK_GE(factor, 1);

    int num = objects().size();
    imgs.resize( num );

    cv::Mat img;
#pragma omp critical
    {
        img = cv::imread(imgFile(), cv::IMREAD_COLOR);
    }

    if (flip()) {
        cv::flip(img, img, 1);
    }

    cv::Mat resizedImg;
    cv::Rect roi, roix;

    for ( int i = 0; i < num; ++i ) {
        const Object &obj(objects()[i]);
        float ratio = std::sqrt( protoarea / obj.bndbox().area() );

        // for simplicity, resize the whole image and then crop
        cv::Size sz(ROUND(img.cols * ratio * factor),
                    ROUND(img.rows * ratio * factor));
        cv::resize(img, resizedImg, sz, 0, 0, cv::INTER_AREA);

        roi.x      = ROUND(obj.bndbox().x() * ratio);
        roi.y      = ROUND(obj.bndbox().y() * ratio);
        roi.width  = ROUND(obj.bndbox().width() * ratio);
        roi.height = ROUND(obj.bndbox().height() * ratio);

        // pad width
        int tmp = roi.width % cellSize;
        if ( tmp ) {
            int l = ROUND(float(cellSize - tmp) / 2.0F);
            roi.x -= l;
            roi.width += cellSize - tmp;
        }

        // pad height
        tmp = roi.height % cellSize;
        if ( tmp ) {
            int t = ROUND(float(cellSize - tmp) / 2.0F);
            roi.y -= t;
            roi.height += cellSize - tmp;
        }

        if ( factor == 1 ) {
            roix = roi;
        } else {
            roix.x      = ROUND(obj.bndbox().x() * ratio * factor);
            roix.y      = ROUND(obj.bndbox().y() * ratio * factor);
            roix.width  = ROUND(obj.bndbox().width() * ratio * factor);
            roix.height = ROUND(obj.bndbox().height() * ratio * factor);

            // pad width
            int tmp = roix.width - roi.width * factor;
            if ( tmp > 0 ) {
                int l = ROUND(tmp / 2.0F);
                roix.x += l;
                roix.width += tmp;
            } else if ( tmp < 0 ) {
                tmp *= -1;
                int l = ROUND(tmp / 2.0F);
                roix.x -= l;
                roix.width += tmp;
            }

            // pad heigth
            tmp = roix.height - roi.height * factor;
            if ( tmp > 0 ) {
                int t = ROUND(tmp / 2.0F);
                roix.y += t;
                roix.height += tmp;
            } else if ( tmp < 0 ) {
                tmp *= -1;
                int t = ROUND(tmp / 2.0F);
                roix.y -= t;
                roix.height += tmp;
            }
        }

        roix.x -= cellSize;
        roix.width += cellSize * 2;

        roix.y -= cellSize;
        roix.height += cellSize * 2;

        imgs[i] = OpencvUtil::subarray(resizedImg, roix, 1);

#if 0
        cv::imshow("resizedObj", imgs[i]);
        cv::waitKey(0);
#endif
    }
}

template<class Archive>
void PosData::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(imgFile_);
    ar & BOOST_SERIALIZATION_NVP(objName_);
    ar & BOOST_SERIALIZATION_NVP(objects_);
    ar & BOOST_SERIALIZATION_NVP(flip_);
    ar & BOOST_SERIALIZATION_NVP(maxArea_);
}

INSTANTIATE_BOOST_SERIALIZATION(PosData);



// ----- NegData -----

NegData::NegData() :
    imgFile_("UNKNOWN"), flip_(false), dataId_(-1)
{
}

cv::Mat NegData::getImg() const
{
    cv::Mat img = cv::imread(imgFile(), cv::IMREAD_COLOR);
    if ( flip() ) {
        cv::flip(img, img, 1);
    }

    return img;
}

const std::string & NegData::imgFile() const
{
    return imgFile_;
}

std::string & NegData::getImgFile()
{
    return imgFile_;
}

bool   NegData::flip() const
{
    return flip_;
}

bool & NegData::getFlip()
{
    return flip_;
}

int   NegData::dataId() const
{
    return dataId_;
}

int & NegData::getDataId()
{
    return dataId_;
}

template<class Archive>
void NegData::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(imgFile_);
    ar & BOOST_SERIALIZATION_NVP(flip_);
    ar & BOOST_SERIALIZATION_NVP(dataId_);
}

INSTANTIATE_BOOST_SERIALIZATION(NegData);


} // namespace RGM
