#ifndef RGM_VOCDATA_HPP_
#define RGM_VOCDATA_HPP_

#include <ostream>
#include <istream>

#include "Data.hpp"
#include "Environment.hpp"
#include "VOCConfig.hpp"

namespace RGM
{
/// PASCAL VOC Constants
static const int VOCObjectClassNum = 20;

static const std::string VOCObjectNames[VOCObjectClassNum+1] = {
    "aeroplane", "bicycle", "bird", "boat", "bottle", "bus", "car", "cat", "chair", "cow",
    "diningtable", "dog", "horse", "motorbike", "person", "pottedplant", "sheep", "sofa",
    "train", "tvmonitor", "Unknown"
};

static const int VOCObjectPoseNum = 5;

static const std::string VOCObjectPoses[5] = {
    "Frontal", "Left", "Rear", "Right", "Unspecified"
};

static const char* VOCImgExtName = ".jpg";

/// PASCAL VOC data
class VOCData
{
public:
    /// Dataset type
    enum DatasetFlag { TRAINING=0, TESTING };

    /// Default constructor
    VOCData();

    /// Constructs the VOCData with given configuration settings
    VOCData(VOCConfig & config);

    /// Sets the configuration settings
    void setVOCConfig(VOCConfig & config);

    /// Read VOC training data
    bool prepareTrainData(const std::string & objName);

    /// Returns the training datasets
    const std::vector<PosData> & posData() const;
    const std::vector<NegData> & negData() const;

    /// Read VOC testing data
    bool prepareTestData(const std::string & objName);

    /// Returns testing image files
    const std::vector<std::string> & testImgFiles() const;

    /// Returns the maximum size of images
    int maxWidth() const;
    int maxHeight() const;

    /// Visualize training data
    void visualize();

    /// Serialization
    void save(const std::string& fileName);

    /// Unserialization
    bool load(const std::string& fileName);

private:
    void readTxt(const std::string& fileName, std::vector<std::string> &posImgs, std::vector<std::string> &negImgs);
    void readPositives(int &dataID, std::vector<std::string> &Imgs);
    Object read(XMLData &annotation, XMLData::XMLNode *pNode);
    void readNegatives(int &dataID, std::vector<std::string> &Imgs);

private:
    VOCConfig * config_;

    DatasetFlag flag_;

    std::vector<PosData> posData_;
    bool flipPositives_;

    std::vector<NegData> negData_;

    std::vector<std::string> testImgFiles_; // testing data

    int maxWidth_; // used for initializing FFTW in inference
    int maxHeight_;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

    DEFINE_RGM_LOGGER;

}; // VOCData

} // namespace RGM

#endif // RGM_VOCDATA_HPP_
