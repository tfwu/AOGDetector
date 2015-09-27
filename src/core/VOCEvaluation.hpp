#ifndef RGM_VOCEVALUATION_HPP_
#define RGM_VOCEVALUATION_HPP_

#include "Environment.hpp"
#include "VOCData.hpp"
#include "AOGrammar.hpp"
#include "ParseTree.hpp"


namespace RGM
{
/// The VOCEvaluation class
class VOCEvaluation
{
public:
    /// typedef
    typedef ParseTree::Detection Detection;

    /// Constructs the evaluation with given configuration file
    /// @param testType, 0: voc, 1: single image, 2: images in a directory
    VOCEvaluation(std::string & configFile, int testType=0);

    /// Evaluation using DP detection
    void evalDP();

    /// Run detection in a testing image
    void detect(XMLData* config);

    /// Run detection for all images in a directory
    void detectBatch(XMLData* config);

private:
    /// Loads the parse trees
    void saveParseTrees(const std::string & fileName);

    /// Loads the parse trees
    bool loadParseTrees(const std::string & fileName);

    /// Loads the detections
    void saveDetections(const std::string & fileName);

    /// Loads the detections
    bool loadDetections(const std::string & fileName);    

    /// Saves the bounding boxes of policy detections
    void saveDetBboxes(const std::string & fileName, const std::vector<std::string> & imgFiles);

private:
    std::string configFile_;

    Environment env_;
    VOCConfig   voc_;
    VOCData     data_;

    AOGrammar   grammar_;
    std::vector<std::vector<ParseTree> >  allpts_;

    std::vector<std::vector<Detection> > allDets_;

    DEFINE_RGM_LOGGER;

}; // class VOCEvaluation

} // namespace RGM

#endif // RGM_VOCEVALUATION_HPP_
