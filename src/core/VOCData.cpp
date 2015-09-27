#include <algorithm>
#include <random>
#include <fstream>

#include <boost/serialization/vector.hpp>

#include "VOCData.hpp"
#include "UtilString.hpp"
#include "UtilFile.hpp"
#include "XMLReader.hpp"

namespace  RGM
{

// ------- VOCData ------

VOCData::VOCData() :
    config_(NULL), flag_(TRAINING), maxWidth_(0), maxHeight_(0), flipPositives_(false)
{
}

VOCData::VOCData(VOCConfig & config) :
    config_(&config), flag_(TRAINING), maxWidth_(0), maxHeight_(0), flipPositives_(false)
{
    RGM_CHECK_NOTNULL(config_);
}

void VOCData::setVOCConfig(VOCConfig & config)
{
    config_ = &config;
}

bool VOCData::prepareTrainData(const std::string & objName)
{
    RGM_CHECK_NOTNULL(config_);

    bool foundObj = config_->findObject(objName);
    if ( !foundObj ) {
        RGM_LOG(error,
                boost::format("Not find training data for: %s") % objName);
        return false;
    }

    flag_ = TRAINING;

    VOCConfig &voc(*config_);
    RGM_LOG(normal,
            boost::format(" Preparing the training data: %s ") % voc.ObjectCategory());

    int dataID = 0;

    std::string fileName = voc.cacheDir_ + voc.ObjectCategory() + "_" +
            voc.trainSetFg_ + "_" + voc.trainSetBg_ + ".bin";
    if (!load(fileName)) {
        // Get positives
        std::string imgSetFgFileName = voc.objectTrainSetFg();
        std::vector<std::string> posImgs;
        std::vector<std::string> negImgs;

        readTxt(imgSetFgFileName, posImgs, negImgs);

        readPositives(dataID, posImgs);

        // Get negatives
        imgSetFgFileName = voc.objectTrainSetBg();
        readTxt(imgSetFgFileName, posImgs, negImgs);

        readNegatives(dataID, negImgs);

        save(fileName);
    }

#ifdef _DEBUG
    //visualize();
#endif

    // shuffle neg
    std::random_shuffle(negData_.begin(), negData_.end());

    RGM_LOG(normal, " and shuffle negData .... Done.");

    return true;
}

const std::vector<NegData> & VOCData::negData() const
{
    return negData_;
}

const std::vector<PosData>  & VOCData::posData() const
{
    return posData_;
}

bool VOCData::prepareTestData(const std::string & objName)
{

    VOCConfig &voc(*config_);
    bool foundObj = voc.findObject(objName);
    if ( !foundObj ) {
        RGM_LOG(error, boost::format("Not find training data for: %s") % objName);
        return false;
    }

    flag_ = TESTING;

    RGM_LOG(normal,
            boost::format(" Preparing the testing data: %s ...") % voc.ObjectCategory());

    std::string fileName = voc.cacheDir_ + voc.ObjectCategory() + "_" + voc.testSet_ + ".bin";

    testImgFiles_.clear();
    maxWidth_ = 0;
    maxHeight_ = 0;

    if (!load(fileName)) {
        std::string imgSetFileName = voc.objectTestSet();

        FILE *f = fopen(imgSetFileName.c_str(), "r");
        if (f == NULL) {
            RGM_LOG(error, boost::format("can not read %s") % imgSetFileName);
            return false;
        }

        char buffer[256];
        int  flag;

        const std::string objectName = voc.ObjectCategory();
        const std::string imgFileDir = voc.imgFileDir_;

        int count = 0;
        while ( !feof(f) ) {
            if ( count % 100 ==0 ) {
                RGM_LOG(normal,
                        boost::format("%s parsing %s testing images (%d)")
                        % voc.project_ % objectName % count);
            }

            fscanf(f, "%s %i\n", buffer, &flag);

            std::string imgFile = imgFileDir + std::string(buffer) + VOCImgExtName;

            testImgFiles_.push_back(imgFile);

            cv::Mat img = cv::imread(imgFile);
            maxWidth_ = std::max<int>(maxWidth_, img.cols);
            maxHeight_ = std::max<int>(maxHeight_, img.rows);

            count++;
        }

        fclose(f);

        save(fileName);
    }

    if (testImgFiles_.size()==0 ) {
        RGM_LOG(error, "No data found.");
    } else {
        RGM_LOG(normal,
                boost::format("    .... got %d images, Done.") % testImgFiles_.size());
    }

    return true;
}

const std::vector<std::string> & VOCData::testImgFiles() const
{
    return testImgFiles_;
}

int VOCData::maxWidth() const
{
    return maxWidth_;
}

int VOCData::maxHeight() const
{
    return maxHeight_;
}

void VOCData::readTxt(const std::string& fileName, std::vector<std::string> &posImgs, std::vector<std::string> &negImgs)
{
    FILE *f = fopen(fileName.c_str(), "r");
    if (f == NULL) {
        RGM_LOG(error,
                boost::format("can not read %s") % fileName);
        return;
    }

    char buffer[256];
    int  flag;

    posImgs.clear();
    negImgs.clear();

    while ( !feof(f) ) {
        fscanf(f, "%s %i\n", buffer, &flag);
        if (flag==-1) {
            negImgs.push_back(std::string(buffer));
        } else {
            posImgs.push_back(std::string(buffer));
        }
    }

    fclose(f);
}

Object VOCData::read(XMLData &annotation, XMLData::XMLNode *pNode)
{
    Object obj;

    std::string pose = annotation.GetString("pose", pNode);

    const std::string* iter = find(VOCObjectPoses, VOCObjectPoses+VOCObjectPoseNum, pose);
    RGM_CHECK_NOTEQ(iter,  VOCObjectPoses+VOCObjectPoseNum);

    obj.getPose() = static_cast<Object::Pose>(iter - VOCObjectPoses);
    obj.getTruncated() = annotation.GetBoolean("truncated", pNode);
    obj.getDifficult() = annotation.GetBoolean("difficult", pNode);

    // top left coordinates of (0, 0) instead of (1, 1)
    obj.getBndbox().setX( annotation.GetInt("bndbox/xmin", pNode)-1 );
    obj.getBndbox().setY( annotation.GetInt("bndbox/ymin", pNode)-1 );
    obj.getBndbox().setRight( annotation.GetInt("bndbox/xmax", pNode)-1 );
    obj.getBndbox().setBottom( annotation.GetInt("bndbox/ymax", pNode)-1 );

    return obj;
}

void VOCData::readPositives(int &dataID, std::vector<std::string> &Imgs)
{
    VOCConfig &voc(*config_);

    const std::string objectName = voc.ObjectCategory();
    const std::string annotationDir = voc.annotationDir_;
    const std::string imgFileDir = voc.imgFileDir_;
    flipPositives_ = voc.flipPositives_;

    posData_.clear();

    for ( int i=0; i<Imgs.size(); ++i ) {
        if ( i % 100 ==0 ) {
            RGM_LOG(normal,
                    boost::format("%s parsing %s positives (%d/%d)")
                    % voc.project_ % objectName % i % Imgs.size());
        }

        std::string annotationFile = annotationDir + Imgs[i] + ".xml";

        // Read the xml file
        XMLData anno;
        anno.ReadFromFile(annotationFile, "annotation");

        // Get image size
        int imgWidth  = anno.GetInt("size/width");
        int imgHeight = anno.GetInt("size/height");

        maxWidth_ = std::max<int>(maxWidth_, imgWidth);
        maxHeight_ = std::max<int>(maxHeight_, imgHeight);

        // Get annotated objects
        std::vector<XMLData::XMLNode *> objectNodes = anno.GetNodes("object");

        std::vector<Object> objs;
        std::vector<Object> flipObjs;
        int maxArea = 0;

        for ( int j=0; j<objectNodes.size(); ++j ) {
            XMLData::XMLNode *objectNode = objectNodes[j];
            std::string name = anno.GetString("name", objectNode);

            if ( name.compare(objectName)==0 ) {
                //std::string imgFile = imgFileDir + Imgs[i] + VOCImgExtName;

                const std::string *iter =
                        std::find(VOCObjectNames, VOCObjectNames+VOCObjectClassNum, name);

                assert((iter != VOCObjectNames + VOCObjectClassNum));

                Object curObject = read(anno, objectNode);

                if ( voc.useDifficultPos_ || !curObject.difficult() ) {

                    curObject.getDataId() = dataID++;
                    objs.push_back(curObject);

                    maxArea = std::max<int>(maxArea, curObject.bndbox().area());

                    if ( flipPositives_ ) {
                        curObject.getBndbox().setX(imgWidth - curObject.bndbox().right() + 1);
                        curObject.getDataId() = dataID++;

                        flipObjs.push_back(curObject);
                    }
                }
            }
        }// for j

        if (objs.size() == 0) {
            continue;
        }

        PosData posImg;
        posImg.getImgFile() = imgFileDir + Imgs[i] + VOCImgExtName;
        posImg.getObjName() = objectName;
        posImg.getObjects().swap(objs);
        posImg.getMaxArea() = maxArea;
        posImg.getFlip() = false;

        posData_.push_back(posImg);

        if (flipPositives_) {
            posImg.getObjects().swap(flipObjs);
            posImg.getFlip() = true;
            posData_.push_back(posImg);
        }
    }
}

void VOCData::readNegatives(int &dataID, std::vector<std::string> &Imgs)
{
    VOCConfig &voc(*config_);

    const std::string objectName = voc.ObjectCategory();
    const std::string annotationDir = voc.annotationDir_;
    const std::string imgFileDir = voc.imgFileDir_;

    negData_.clear();

    for ( int i=0; i<Imgs.size(); ++i ) {
        if ( i % 100 ==0 ) {
            RGM_LOG(normal,
                    boost::format("%s parsing %s negatives (%d/%d)")
                    % voc.project_ % objectName % i % Imgs.size());
        }

        std::string annotationFile = annotationDir + Imgs[i] + ".xml";

        // Read the xml file
        XMLData anno;
        anno.ReadFromFile(annotationFile, "annotation");

        // Get image size
        int imgWidth  = anno.GetInt("size/width");
        int imgHeight = anno.GetInt("size/height");

        maxWidth_ = std::max<int>(maxWidth_, imgWidth);
        maxHeight_ = std::max<int>(maxHeight_, imgHeight);

        NegData neg;

        neg.getImgFile() = imgFileDir + Imgs[i] + VOCImgExtName;;
        neg.getFlip() = false;
        neg.getDataId() = dataID++;

        negData_.push_back(neg);
    }
}

void VOCData::visualize()
{
    if (posData().size()==0) {
        return;
    }

    for ( int i = 0;  i < posData().size(); ++i ) {
        cv::Mat img = cv::imread(posData()[i].imgFile(), cv::IMREAD_COLOR);
        if ( posData()[i].flip() ) {
            cv::flip(img, img, 1);
        }

        for ( int j = 0; j < posData()[i].objects().size(); ++j ) {
            cv::rectangle(img, posData()[i].objects()[j].bndbox().cvRect(), cv::Scalar(0, 0, 255), 2);
        }

        cv::imshow("VOCData: Press ESC to quit", img);
        if ( 27 == cv::waitKey(0) ) {
            return;
        }
    }
}

template<class Archive>
void VOCData::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(flag_);
    if ( flag_ == TRAINING ) {
        ar & BOOST_SERIALIZATION_NVP(posData_);
        ar & BOOST_SERIALIZATION_NVP(flipPositives_);
        ar & BOOST_SERIALIZATION_NVP(negData_);
    } else {
        ar & BOOST_SERIALIZATION_NVP(testImgFiles_);
    }

    ar & BOOST_SERIALIZATION_NVP(maxWidth_);
    ar & BOOST_SERIALIZATION_NVP(maxHeight_);
}

INSTANTIATE_BOOST_SERIALIZATION(VOCData);


void VOCData::save(const std::string& fileName)
{
    std::ofstream os(fileName.c_str(), std::ios::out);
    if (!os.is_open()) {
        RGM_LOG(error,
                boost::format("can not write data to %s") % fileName);
        return;
    }

    boost::archive::binary_oarchive oa(os);

    oa << *this;
}

bool VOCData::load(const std::string& fileName)
{
    std::ifstream is(fileName.c_str(), std::ios::in);
    if ( !is.is_open() ) {
        RGM_LOG(error,
                boost::format("can not read data from %s") % fileName);
        return false;
    }

    boost::archive::binary_iarchive ia(is);

    ia >> *this;

    return true;
}

} // namespace  RGM
