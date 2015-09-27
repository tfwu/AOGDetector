#include <fstream>
#include <algorithm>
#include <numeric>
#include <map>

#include "AOGrammar.hpp"
#include "UtilOpencv.hpp"
#include "UtilFile.hpp"
#include "UtilString.hpp"
#include "UtilGeneric.hpp"

namespace RGM
{

// ------ AOGrammar's Edge -------

Edge::Edge()
{
    init();
}

Edge::Edge(const Edge & e)
{
    init();

    edgeType_ = e.type();
    isLRFlip_ = e.isLRFlip();
    idx_ = e.idx();
}

Edge::Edge(edgeType t, Node * fromNode, Node * toNode)
{    
    RGM_CHECK_NOTEQ(t, UNKNOWN_EDGE);
    RGM_CHECK_NOTNULL(toNode);
    RGM_CHECK_NOTNULL(fromNode);
    RGM_CHECK_NOTEQ(fromNode, toNode);

    init();

    edgeType_ = t;
    fromNode_ = fromNode;
    toNode_ = toNode;
}

Edge::~Edge()
{
    init();
}

void Edge::init()
{
    edgeType_ = UNKNOWN_EDGE;
    fromNode_ = NULL;
    toNode_ = NULL;
    isLRFlip_ = false;
    LRMirrorEdge_ = NULL;
    idx_.fill(-1);
}

Edge::edgeType Edge::type() const
{
    return edgeType_;
}

Edge::edgeType & Edge::getType()
{
    return edgeType_;
}

const Node * Edge::fromNode() const
{
    return fromNode_;
}

Node *& Edge::getFromNode()
{
    return fromNode_;
}

const Node * Edge::toNode() const
{
    return toNode_;
}

Node *& Edge::getToNode()
{
    return toNode_;
}

bool Edge::isLRFlip() const
{
    return isLRFlip_;
}

bool & Edge::getIsLRFlip()
{
    return isLRFlip_;
}

const Edge * Edge::lrMirrorEdge() const
{
    return LRMirrorEdge_;
}

Edge *& Edge::getLRMirrorEdge()
{
    return LRMirrorEdge_;
}

const Edge::Index & Edge::idx() const
{
    return idx_;
}

Edge::Index & Edge::getIdx()
{
    return idx_;
}

void Edge::assignIdx(AOGrammar * g)
{
    RGM_CHECK_NOTNULL(g);

    getIdx()(IDX_FROM) = g->idxNode(fromNode());
    getIdx()(IDX_TO) = g->idxNode(toNode());
    getIdx()(IDX_MIRROR) = g->idxEdge(lrMirrorEdge());

    RGM_CHECK_GE(idx()(IDX_FROM), 0);
    RGM_CHECK_GE(idx()(IDX_TO), 0);
}

void Edge::assignConnections(AOGrammar * g)
{
    RGM_CHECK_NOTNULL(g);

    // from Node
    getFromNode() = g->findNode( idx()(IDX_FROM) );

    // to Node
    getToNode() = g->findNode( idx()(IDX_TO) );

    // flip
    getLRMirrorEdge() = g->findEdge( idx()(IDX_MIRROR) );
}

template <class Archive>
void Edge::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(edgeType_);
    ar & BOOST_SERIALIZATION_NVP(isLRFlip_);
    ar & BOOST_SERIALIZATION_NVP(idx_);
}

INSTANTIATE_BOOST_SERIALIZATION(Edge);




// ------ AOGrammar_Node ------

Node::Node()
{
    init();
}

Node::Node(nodeType t)
{
    init();
    nodeType_ = t;
}

Node::Node(const Node & n)
{
    init();

    nodeType_ = n.type();
    isLRFlip_ = n.isLRFlip();
    detectWindow_ = n.detectWindow();
    anchor_ = n.anchor();

    idxInEdge_ = n.idxInEdge();

    idxOutEdge_ = n.idxOutEdge();

    idx_ = n.idx();
}

Node::~Node()
{
    init();
}

void Node::init()
{
    nodeType_ = UNKNOWN_NODE;
    isLRFlip_ = false;
    LRMirrorNode_ = NULL;
    anchor_.setZero();
    scaleprior_ = NULL;
    offset_ = NULL;
    deformation_ = NULL;
    appearance_ = NULL;
    cachedFFTFilter_ = NULL;
#if RGM_USE_PCA_DIM
    wpca_ = NULL;
    cachedFFTPCAFilter_ = NULL;
#endif

    tag_ = boost::uuids::random_generator()();

    idxInEdge_.clear();
    idxOutEdge_.clear();
    idx_.fill(-1);
}

Node::nodeType Node::type() const
{
    return nodeType_;
}

Node::nodeType & Node::getType()
{
    return nodeType_;
}

const std::vector<Edge *> & Node::inEdges() const
{
    return inEdges_;
}

std::vector<Edge *> & Node::getInEdges()
{
    return inEdges_;
}

const std::vector<Edge *> & Node::outEdges() const
{
    return outEdges_;
}

std::vector<Edge *> & Node::getOutEdges()
{
    return outEdges_;
}

bool Node::isLRFlip() const
{
    return isLRFlip_;
}

bool & Node::getIsLRFlip()
{
    return isLRFlip_;
}

const Node * Node::lrMirrorNode() const
{
    return LRMirrorNode_;
}

Node *& Node::getLRMirrorNode()
{
    return LRMirrorNode_;
}

const Rectangle2i & Node::detectWindow() const
{
    return detectWindow_;
}

Rectangle2i & Node::getDetectWindow()
{
    return detectWindow_;
}

const Scaleprior * Node::scaleprior() const
{
    return scaleprior_;
}

Scaleprior *& Node::getScaleprior()
{
    return scaleprior_;
}

const Offset * Node::offset() const
{
    return offset_;
}

Offset *& Node::getOffset()
{
    return offset_;
}

const Node::Anchor &  Node::anchor() const
{
    return anchor_;
}

Node::Anchor & Node::getAnchor()
{
    return anchor_;
}

const Deformation *  Node::deformation() const
{
    return deformation_;
}

Deformation *& Node::getDeformation()
{
    return deformation_;
}

Deformation::Param Node::deformationParam() const
{
    if ( !isLRFlip() ) {
        return deformation()->w();
    }

    Deformation::Param w = deformation()->w();
    w(1) *= -1.0F;

    return w;
}

const Appearance * Node::appearance() const
{
    return appearance_;
}

Appearance *& Node::getAppearance()
{
    return appearance_;
}

Appearance::Param Node::appearanceParam() const
{
    if ( !isLRFlip() ) {
        return appearance()->w();
    }

    return FeaturePyramid::Flip( appearance()->w() );
}

#if RGM_USE_PCA_DIM
const Appearance::PCAParam * Node::wpca() const
{
    return wpca_;
}

Appearance::PCAParam *& Node::getWpca()
{
    return wpca_;
}

const Node::FFTPCAFilter *  Node::cachedFFTPCAFilter() const
{
    return cachedFFTPCAFilter_;
}

Node::FFTPCAFilter *& Node::getCachedFFTPCAFilter()
{
    return cachedFFTPCAFilter_;
}
#endif

const Node::FFTFilter * Node::cachedFFTFilter() const
{
    return cachedFFTFilter_;
}

Node::FFTFilter *& Node::getCachedFFTFilter()
{
    return cachedFFTFilter_;
}

const boost::uuids::uuid& Node::tag() const
{
    return tag_;
}

const std::vector<int> & Node::idxInEdge() const
{
    return idxInEdge_;
}

std::vector<int> & Node::getIdxInEdge()
{
    return idxInEdge_;
}

const std::vector<int> & Node::idxOutEdge() const
{
    return idxOutEdge_;
}

std::vector<int> & Node::getIdxOutEdge()
{
    return idxOutEdge_;
}

const Node::Index &  Node::idx() const
{
    return idx_;
}

Node::Index &  Node::getIdx()
{
    return idx_;
}

void Node::assignIdx(AOGrammar * g)
{
    RGM_CHECK_NOTNULL(g);

    int num = inEdges().size();
    getIdxInEdge().resize(num);
    for ( int i = 0; i < num; ++i ) {
        getIdxInEdge()[i] = g->idxEdge(inEdges()[i]);
    }

    num = outEdges().size();
    getIdxOutEdge().resize(num);
    for ( int i = 0; i < num; ++i ) {
        getIdxOutEdge()[i] = g->idxEdge(outEdges()[i]);
    }

    getIdx()(IDX_MIRROR) = g->idxNode(lrMirrorNode());
    getIdx()(IDX_SCALEPRIOR) = g->idxScaleprior(scaleprior());
    getIdx()(IDX_BIAS) = g->idxOffset(offset());
    getIdx()(IDX_DEF) = g->idxDeformation(deformation());
    getIdx()(IDX_APP) = g->idxAppearance(appearance());
    getIdx()(IDX_MIRROR) = g->idxNode(lrMirrorNode());
    getIdx()(IDX_FILTER) = g->idxFFTFilter(cachedFFTFilter());
#if RGM_USE_PCA_DIM
    getIdx()(IDX_WPCA) = g->idxWpca(wpca());
    getIdx()(IDX_PCAFILTER) = g->idxFFTPCAFilter(cachedFFTPCAFilter());
#endif
}

void Node::assignConnections(AOGrammar * g)
{
    RGM_CHECK_NOTNULL(g);

    // Assigning inEdges
    getInEdges().resize(idxInEdge_.size(), NULL);
    for ( int i = 0; i < idxInEdge_.size(); ++i ) {
        getInEdges()[i] = g->findEdge( idxInEdge_[i] );
    }

    // Assigning outEdges
    getOutEdges().resize(idxOutEdge_.size(), NULL);
    for ( int i = 0; i < idxOutEdge_.size(); ++i ) {
        getOutEdges()[i] = g->findEdge( idxOutEdge_[i] );
    }

    // Assigning LR mirror node
    getLRMirrorNode() = g->findNode(idx()(IDX_MIRROR));

    // scale prior
    getScaleprior() = g->findScaleprior( idx()(IDX_SCALEPRIOR) );

    // offset
    getOffset() = g->findOffset( idx()(IDX_BIAS) );

    // deformation
    getDeformation() = g->findDeformation(idx()(IDX_DEF));

    // appearance
    getAppearance() = g->findAppearance(idx()(IDX_APP));

    // fft
    getCachedFFTFilter() = g->findCachedFFTFilter(idx()(IDX_FILTER));

#if RGM_USE_PCA_DIM
    getWpca() = g->findWpca(idx()(IDX_WPCA));
    getCachedFFTPCAFilter() = g->findCachedFFTPCAFilter(idx()(IDX_PCAFILTER));
#endif
}

template<class Archive>
void Node::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(nodeType_);
    ar & BOOST_SERIALIZATION_NVP(isLRFlip_);
    ar & BOOST_SERIALIZATION_NVP(detectWindow_);
    ar & BOOST_SERIALIZATION_NVP(anchor_);
    ar & BOOST_SERIALIZATION_NVP(idxInEdge_);
    ar & BOOST_SERIALIZATION_NVP(idxOutEdge_);
    ar & BOOST_SERIALIZATION_NVP(idx_);
}

INSTANTIATE_BOOST_SERIALIZATION(Node);



// ------ AOGrammar::FeatureExtraction ------

AOGrammar::FeatureExtraction::FeatureExtraction() :
    cellSize_(8), extraOctave_(false), featureBias_(10.0F), interval_(10)
{
}

void AOGrammar::FeatureExtraction::init()
{
    cellSize_ = 8;
    extraOctave_ = false;
    featureBias_ = 10.0F;
    interval_ = 10;
}

template<class Archive>
void AOGrammar::FeatureExtraction::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(cellSize_);
    ar & BOOST_SERIALIZATION_NVP(extraOctave_);
    ar & BOOST_SERIALIZATION_NVP(featureBias_);
    ar & BOOST_SERIALIZATION_NVP(interval_);
}

INSTANTIATE_BOOST_SERIALIZATION(AOGrammar::FeatureExtraction);



// ------ AOGrammar ------

AOGrammar::AOGrammar()
{
    init();
}

AOGrammar::~AOGrammar()
{
    clear();
}

AOGrammar::AOGrammar(const std::string & modelFile)
{
    init();
    read(modelFile);
}

void AOGrammar::init()
{
    gType_ = UNKNOWN_GRAMMAR;
    regMethod_ = REG_L2;
    name_ = "UNKNOWN";
    note_ = "UNKNOWN";
    year_ = "UNKNOWN";
    isSingleObjModel_ = true;
    rootNode_ = NULL;
    isLRFlip_ = false;
    maxDetectWindow_ = Rectangle2i(0, 0);
    minDetectWindow_ = Rectangle2i(0, 0);
    featureExtraction_.init();
    thresh_ = 0.0F;
    cached_ = false;
#if RGM_USE_PCA_DIM
    pcaCoef_ = NULL;
#endif
    idxRootNode_ = -1;

}

void AOGrammar::clear()
{
    gType_ = UNKNOWN_GRAMMAR;
    regMethod_ = REG_L2;
    name_  = "UNKNOWN";
    note_  = "UNKNOWN";
    year_  = "UNKNOWN";
    isSingleObjModel_ = true;

    // node set
    for ( int i = 0; i < nodeSet_.size(); ++i ) {
        delete nodeSet_[i];
    }
    nodeSet_.clear();

    nodeDFS_.clear();
    nodeBFS_.clear();

    compNodeDFS_.clear();
    compNodeBFS_.clear();

    // root
    rootNode_ = NULL;

    isLRFlip_ = false;

    // edge set
    for ( int i = 0; i < edgeSet_.size(); ++i ) {
        delete edgeSet_[i];
    }
    edgeSet_.clear();

    // appearance
    for ( int i = 0; i < appearanceSet_.size(); ++i ) {
        delete appearanceSet_[i];
    }
    appearanceSet_.clear();

    // offset
    for ( int i = 0; i < biasSet_.size(); ++i ) {
        delete biasSet_[i];
    }
    biasSet_.clear();

    // deformation
    for ( int i = 0; i < deformationSet_.size(); ++i ) {
        delete deformationSet_[i];
    }
    deformationSet_.clear();

    // scale prior
    for ( int i = 0; i < scalepriorSet_.size(); ++i ) {
        delete scalepriorSet_[i];
    }
    scalepriorSet_.clear();

    maxDetectWindow_ = Rectangle2i();
    minDetectWindow_ = Rectangle2i();

    cached_ = false;

    for ( int i = 0; i < cachedFFTFilters_.size(); ++i ) {
        delete cachedFFTFilters_[i];
    }
    cachedFFTFilters_.clear();
#if RGM_USE_PCA_DIM
    if (pcaCoef_ != NULL ) {
        delete pcaCoef_;
    }

    for ( int i = 0; i < pcaAppearanceParam_.size(); ++i ) {
        if (pcaAppearanceParam_[i] != NULL) {
            delete pcaAppearanceParam_[i];
        }
    }
    pcaAppearanceParam_.clear();

    for ( int i = 0; i < cachedFFTPCAFilters_.size(); ++i ) {
        delete cachedFFTPCAFilters_[i];
    }
    cachedFFTPCAFilters_.clear();
#endif

    bboxPred_.clear();

    idxRootNode_ = -1;

}

bool AOGrammar::empty() const
{
    return nodeSet().size() == 0;
}

AOGrammar::grammarType AOGrammar::type() const
{
    return gType_;
}

AOGrammar::grammarType & AOGrammar::getType()
{
    return gType_;
}

const std::string & AOGrammar::name() const
{
    return name_;
}

std::string & AOGrammar::getName()
{
    return name_;
}

const std::string & AOGrammar::note() const
{
    return note_;
}

std::string & AOGrammar::getNote()
{
    return note_;
}

const std::string & AOGrammar::year() const
{
    return year_;
}

std::string & AOGrammar::getYear()
{
    return year_;
}

bool AOGrammar::isSingleObjModel() const
{
    return isSingleObjModel_;
}

bool & AOGrammar::getIsSingleObjModel()
{
    return isSingleObjModel_;
}

const std::vector<Node *> & AOGrammar::nodeSet() const
{
    return nodeSet_;
}

std::vector<Node *> & AOGrammar::getNodeSet()
{
    return nodeSet_;
}

void AOGrammar::traceNodeDFS(Node * curNode, std::vector<int> & visited, std::vector<Node *> & nodeDFS)
{
    // assume visited.size() == nodeSet.size()

    int idx = idxNode(curNode);
    if (visited[idx] == 1 ) {
        RGM_LOG(error, "Cycle detected in grammar!");
        return;
    }

    visited[idx] = 1;

    int numChild = curNode->outEdges().size();

    for ( int i=0; i<numChild; ++i ) {
        Node * ch = curNode->getOutEdges()[i]->getToNode();
        int chIdx = idxNode( ch );
        if ( visited[chIdx] < 2 ) {
            traceNodeDFS(ch, visited, nodeDFS);
        }
    }

    nodeDFS.push_back(curNode);
    visited[idx] = 2;
}

const std::vector<Node *> & AOGrammar::nodeDFS() const
{
    return nodeDFS_;
}

std::vector<Node *> & AOGrammar::getNodeDFS()
{
    return nodeDFS_;
}

const std::vector<std::vector<Node *> > & AOGrammar::compNodeDFS() const
{
    return compNodeDFS_;
}

std::vector<std::vector<Node *> > & AOGrammar::getCompNodeDFS()
{
    return compNodeDFS_;
}

void AOGrammar::traceNodeBFS(Node * curNode, std::vector<int> & visited, std::vector<Node *> & nodeBFS)
{
    // assume visited.size() == nodeSet.size()

    int idx = idxNode(curNode);
    if (visited[idx] == 1 ) {
        RGM_LOG(error, "Cycle detected in grammar!");
        return;
    }

    nodeBFS.push_back(curNode);

    visited[idx] = 1;

    int numChild = curNode->outEdges().size();

    for ( int i=0; i<numChild; ++i ) {
        Node * ch = curNode->getOutEdges()[i]->getToNode();
        int chIdx = idxNode( ch );
        if ( visited[chIdx] < 2 ) {
            traceNodeBFS(ch, visited, nodeBFS);
        }
    }

    visited[idx] = 2;
}

const std::vector<Node *> & AOGrammar::nodeBFS() const
{
    return nodeBFS_;
}

std::vector<Node *> & AOGrammar::getNodeBFS()
{
    return nodeBFS_;
}

const std::vector<std::vector<Node *> > & AOGrammar::compNodeBFS() const
{
    return compNodeBFS_;
}

std::vector<std::vector<Node *> > & AOGrammar::getCompNodeBFS()
{
    return compNodeBFS_;
}


void AOGrammar::traceCompNodeDFSandBFS()
{
    int numComp = rootNode()->outEdges().size();
    if ( isLRFlip() ) {
        numComp /= 2;
    }

    getCompNodeDFS().resize(numComp);
    getCompNodeBFS().resize(numComp);

    for ( int i = 0, j = 0; i < getRootNode()->getOutEdges().size(); ++i ) {
        Node * node = getRootNode()->getOutEdges()[i]->getToNode();
        if ( node->isLRFlip() ) {
            continue;
        }

        std::vector<int> visited(nodeSet().size(), 0);

        getCompNodeDFS()[j].clear();
        traceNodeDFS(node, visited, getCompNodeDFS()[j]);

        visited.assign(visited.size(), 0);

        getCompNodeBFS()[j].clear();
        traceNodeBFS(node, visited, getCompNodeBFS()[j]);

        j++;
    }
}


const std::vector<Edge *> & AOGrammar::edgeSet() const
{
    return edgeSet_;
}

std::vector<Edge *> & AOGrammar::getEdgeSet()
{
    return edgeSet_;
}

const Node *  AOGrammar::rootNode() const
{
    return rootNode_;
}

Node *& AOGrammar::getRootNode()
{
    return rootNode_;
}

bool AOGrammar::isLRFlip() const
{
    return isLRFlip_;
}

bool & AOGrammar::getIsLRFlip()
{
    return isLRFlip_;
}

AOGrammar::regType   AOGrammar::regMethod() const
{
    return regMethod_;
}

AOGrammar::regType & AOGrammar::getRegMethod()
{
    return regMethod_;
}

const std::vector<Appearance *>  & AOGrammar::appearanceSet() const
{
    return appearanceSet_;
}

std::vector<Appearance *>  & AOGrammar::getAppearanceSet()
{
    return appearanceSet_;
}

const std::vector<Offset *>  & AOGrammar::biasSet() const
{
    return biasSet_;
}

std::vector<Offset *>  & AOGrammar::getBiasSet()
{
    return biasSet_;
}

const std::vector<Deformation *>  & AOGrammar::deformationSet() const
{
    return deformationSet_;
}

std::vector<Deformation *>  & AOGrammar::getDeformationSet()
{
    return deformationSet_;
}

const std::vector<Scaleprior *> & AOGrammar::scalepriorSet() const
{
    return scalepriorSet_;
}

std::vector<Scaleprior *> & AOGrammar::getScalepriorSet()
{
    return scalepriorSet_;
}

const Rectangle2i & AOGrammar::maxDetectWindow() const
{
    return maxDetectWindow_;
}

Rectangle2i & AOGrammar::getMaxDetectWindow()
{
    return maxDetectWindow_;
}

const Rectangle2i & AOGrammar::minDetectWindow() const
{
    return minDetectWindow_;
}

Rectangle2i & AOGrammar::getMinDetectWindow()
{
    return minDetectWindow_;
}

const std::vector<Node::FFTFilter *> &  AOGrammar::cachedFFTFilters() const
{
    return cachedFFTFilters_;
}

std::vector<Node::FFTFilter *> &  AOGrammar::getCachedFFTFilters()
{
    return cachedFFTFilters_;
}

bool AOGrammar::cachedFFTStatus() const
{
    return cached_;
}

bool & AOGrammar::getCachedFFTStatus()
{
    return cached_;
}

void AOGrammar::cachingFFTFilters(bool withPCA)
{
    if ( cachedFFTStatus() ) {
        return;
    }

    // release old ones
    for ( int i = 0; i < cachedFFTFilters().size(); ++i ) {
        delete getCachedFFTFilters()[i];
    }
    getCachedFFTFilters().clear();

    std::vector<Appearance::Param> w;
    for ( int i = 0; i < nodeSet().size(); ++i ) {
        if (nodeSet()[i]->type() == Node::T_NODE) {
            w.push_back( nodeSet()[i]->appearanceParam() );
        }
    }

    getCachedFFTFilters().resize(w.size(), NULL);

#pragma omp parallel for
    for ( int i = 0; i < w.size(); ++i ) {
        getCachedFFTFilters()[i] = new Node::FFTFilter();
        Patchwork::TransformFilter(w[i], *(getCachedFFTFilters()[i]));
    }

#if RGM_USE_PCA_DIM
    if (withPCA) {
        cachingFFTPCAFilters();
    }
#endif

    getCachedFFTStatus() = true;
}

#if RGM_USE_PCA_DIM
bool AOGrammar::project(const std::string & pcaCoefFile)
{
    FILE * pf = fopen(pcaCoefFile.c_str(), "rb");
    if (pf == NULL) {
        RGM_LOG(error, "Can not read pca data from " + pcaCoefFile);
        return false;
    }

    int rows, cols;
    int count = fread(&rows, sizeof(int), 1, pf);
    RGM_CHECK_EQ(count, 1);

    count = fread(&cols, sizeof(int), 1, pf);
    RGM_CHECK_EQ(count, 1);

    const int dim = FeaturePyramid::NbFeatures;
    const int notrucatDim = dim - 1;
    if (rows != notrucatDim || cols != notrucatDim) {
        RGM_LOG(error, " The input pca data is not valid.\n");
        fclose(pf);
        return false;
    }

    Matrix coef(notrucatDim, notrucatDim);
    int num = notrucatDim * notrucatDim;
    count = fread(coef.data(), sizeof(Scalar), num, pf);
    RGM_CHECK_EQ(count, num);

    fclose(pf);

    if (pcaCoef_ == NULL) {
        pcaCoef_ = new Matrix();
    }

    const int pcaDim = FeaturePyramid::NbPCAFeatures;

    *pcaCoef_ = Matrix::Zero(dim, pcaDim);
    (*pcaCoef_)(dim-1, pcaDim-1) = 1;

    (*pcaCoef_).block(0, 0, coef.rows(), RGM_USE_PCA_DIM) =
            coef.block(0, 0, coef.rows(), RGM_USE_PCA_DIM);

    pcaCoef_->transposeInPlace();

    //std::cout << _coef << std::endl;

    for ( int i = 0; i < pcaAppearanceParam().size(); ++i ) {
        if (pcaAppearanceParam()[i] != NULL) {
            delete getPcaAppearanceParam()[i];
        }
    }
    getPcaAppearanceParam().clear();

    for ( int i = 0; i < nodeSet().size(); ++i ) {
        if (nodeSet()[i]->type() == Node::T_NODE) {
            Appearance::Param w = nodeSet()[i]->appearanceParam();
            getPcaAppearanceParam().push_back(new Appearance::PCAParam());
            Appearance::PCAParam & wpca(*(getPcaAppearanceParam().back()));
            wpca.resize(w.rows(), w.cols());
            for ( int r = 0; r < w.rows(); ++r ) {
                for ( int c = 0; c < w.cols(); ++c ) {
                    wpca(r, c) = (*pcaCoef_) * w(r, c).matrix();
                } // for c
            } // for r

            getNodeSet()[i]->getWpca() = getPcaAppearanceParam().back();
        }
    } // for i

    return true;
}


const Matrix * AOGrammar::pcaCoef() const
{
    return pcaCoef_;
}

Matrix *& AOGrammar::getPcaCoef()
{
    return pcaCoef_;
}

const std::vector<Appearance::PCAParam *>  & AOGrammar::pcaAppearanceParam() const
{
    return pcaAppearanceParam_;
}

std::vector<Appearance::PCAParam *>  & AOGrammar::getPcaAppearanceParam()
{
    return pcaAppearanceParam_;
}

const std::vector<Node::FFTPCAFilter *> &  AOGrammar::cachedFFTPCAFilters() const
{
    return cachedFFTPCAFilters_;
}

std::vector<Node::FFTPCAFilter *> &  AOGrammar::getCachedFFTPCAFilters()
{
    return cachedFFTPCAFilters_;
}

void AOGrammar::cachingFFTPCAFilters()
{
    int num = pcaAppearanceParam().size();
    if (num == 0) {
        return;
    }

    // release old ones
    for ( int i = 0; i < cachedFFTPCAFilters().size(); ++i ) {
        delete getCachedFFTPCAFilters()[i];
    }
    getCachedFFTPCAFilters().clear();

    getCachedFFTPCAFilters().resize(num, NULL);

#pragma omp parallel for
    for ( int i = 0; i < num; ++i ) {
        getCachedFFTPCAFilters()[i] = new Node::FFTPCAFilter();
        Patchwork::TransformPCAFilter(*(pcaAppearanceParam()[i]), *(getCachedFFTPCAFilters()[i]));
    }
}
#endif

const AOGrammar::FeatureExtraction & AOGrammar::featureExtraction() const
{
    return featureExtraction_;
}

int   AOGrammar::cellSize() const
{
    return featureExtraction_.cellSize_;
}

int & AOGrammar::getCellSize()
{
    return featureExtraction_.cellSize_;
}

int AOGrammar::minCellSize()
{
    return extraOctave() ? cellSize() / 4 : cellSize() / 2;
}

bool   AOGrammar::extraOctave() const
{
    return featureExtraction_.extraOctave_;
}

bool & AOGrammar::getExtraOctave()
{
    return featureExtraction_.extraOctave_;
}

Scalar   AOGrammar::featureBias() const
{
    return featureExtraction_.featureBias_;
}

Scalar & AOGrammar::getFeatureBias()
{
    return featureExtraction_.featureBias_;
}

int   AOGrammar::interval() const
{
    return featureExtraction_.interval_;
}

int & AOGrammar::getInterval()
{
    return featureExtraction_.interval_;
}

int   AOGrammar::padx() const
{
    return maxDetectWindow().width();
}

int   AOGrammar::pady() const
{
    return maxDetectWindow().height();
}

Scalar   AOGrammar::thresh() const
{
    return thresh_;
}

Scalar & AOGrammar::getThresh()
{
    return thresh_;
}

const std::map<int, Matrix> & AOGrammar::bboxPred() const
{
    return bboxPred_;
}

std::map<int, Matrix> & AOGrammar::getBboxPred()
{
    return bboxPred_;
}

int AOGrammar::dim() const
{
    int d = 0;

    // Count all the parameters using DFS
    const std::vector<Node * > & DFS = nodeDFS();

    for ( int i = 0; i < DFS.size(); ++i ) {
        const Node * n = DFS[i];
        if ( n->isLRFlip() ) {
            continue;
        }

        Node::nodeType t = n->type();
        switch (t) {
        case Node::T_NODE: {
            const Appearance::Param & w = n->appearance()->w();
            d += static_cast<int>(w.size()) * FeaturePyramid::NbFeatures;
            break;
        }
        case Node::AND_NODE: {
            if ( n->deformation() != NULL ) {
                d += 4;
            }

            if ( n->scaleprior() != NULL ) {
                d += 3;
            }

            if ( n->offset() != NULL ) {
                d++;
            }

            break;
        }
        } // switch
    } // for i

    return d;
}

int AOGrammar::idxNode(const Node * node) const
{
    if ( node == NULL ) {
        return -1;
    }

    for ( int i = 0; i < nodeSet().size(); ++i ) {
        if ( node == nodeSet()[i] ) {
            return i;
        }
    }

    return -1;
}

int AOGrammar::idxObjAndNodeOfTermNode(const Node *tnode) const
{
    if ( tnode == NULL || tnode->type() != Node::T_NODE ) {
        return -1;
    }

    // goes up 3 layes: tnode -> and-node -> or-node -> obj and-node
    const Node * obj = tnode->inEdges()[0]->fromNode()->inEdges()[0]->fromNode()->inEdges()[0]->fromNode();

    return idxNode(obj);
}

int AOGrammar::idxEdge(const Edge * edge) const
{
    if ( edge == NULL ) {
        return -1;
    }

    for ( int i = 0; i < edgeSet().size(); ++i ) {
        if ( edge == edgeSet()[i] ) {
            return i;
        }
    }

    return -1;
}

int AOGrammar::idxAppearance(const Appearance * app) const
{
    if ( app == NULL ) {
        return -1;
    }

    for ( int i = 0; i < appearanceSet().size(); ++i ) {
        if ( app == appearanceSet()[i] ) {
            return i;
        }
    }

    return -1;
}

int AOGrammar::idxOffset(const Offset * off) const
{
    if ( off == NULL ) {
        return -1;
    }

    for ( int i = 0; i < biasSet().size(); ++i ) {
        if ( off == biasSet()[i] ) {
            return i;
        }
    }

    return -1;
}

int AOGrammar::idxDeformation(const Deformation * def) const
{
    if ( def == NULL ) {
        return -1;
    }

    for ( int i = 0; i < deformationSet().size(); ++i ) {
        if ( def == deformationSet()[i] ) {
            return i;
        }
    }

    return -1;
}

int AOGrammar::idxScaleprior(const Scaleprior * scale) const
{
    if ( scale == NULL ) {
        return -1;
    }

    for ( int i = 0; i < scalepriorSet().size(); ++i ) {
        if ( scale == scalepriorSet()[i] ) {
            return i;
        }
    }

    return -1;
}

int AOGrammar::idxFFTFilter(const Node::FFTFilter *  filter) const
{
    if ( filter == NULL ) {
        return -1;
    }

    for ( int i = 0; i < cachedFFTFilters().size(); ++i ) {
        if ( filter == cachedFFTFilters()[i] ) {
            return i;
        }
    }

    return -1;
}

#if RGM_USE_PCA_DIM
int AOGrammar::idxWpca(const Appearance::PCAParam *   w) const
{
    if ( w == NULL ) {
        return -1;
    }

    for ( int i = 0; i < pcaAppearanceParam().size(); ++i ) {
        if ( w == pcaAppearanceParam()[i] ) {
            return i;
        }
    }

    return -1;
}

int AOGrammar::idxFFTPCAFilter(const Node::FFTPCAFilter *   filter) const
{
    if ( filter == NULL ) {
        return -1;
    }

    for ( int i = 0; i < cachedFFTPCAFilters().size(); ++i ) {
        if ( filter == cachedFFTPCAFilters()[i] ) {
            return i;
        }
    }

    return -1;
}
#endif

Node * AOGrammar::findNode(int idx)
{
    if (idx < 0 || idx >= nodeSet().size()) {
        return NULL;
    }

    return getNodeSet()[idx];
}

Edge * AOGrammar::findEdge(int idx)
{
    if (idx < 0 || idx >= edgeSet().size()) {
        return NULL;
    }

    return getEdgeSet()[idx];
}

Appearance * AOGrammar::findAppearance(int idx)
{
    if (idx < 0 || idx >= appearanceSet().size()) {
        return NULL;
    }

    return getAppearanceSet()[idx];
}

Offset * AOGrammar::findOffset(int idx)
{
    if (idx < 0 || idx >= biasSet().size()) {
        return NULL;
    }

    return getBiasSet()[idx];
}

Deformation * AOGrammar::findDeformation(int idx)
{
    if (idx < 0 || idx >= deformationSet().size()) {
        return NULL;
    }

    return getDeformationSet()[idx];
}

Scaleprior * AOGrammar::findScaleprior(int idx)
{
    if (idx < 0 || idx >= scalepriorSet().size()) {
        return NULL;
    }

    return getScalepriorSet()[idx];
}


Node::FFTFilter * AOGrammar::findCachedFFTFilter(int idx)
{
    if (idx < 0 || idx >= cachedFFTFilters().size()) {
        return NULL;
    }

    return getCachedFFTFilters()[idx];
}

#if RGM_USE_PCA_DIM
Appearance::PCAParam *  AOGrammar::findWpca(int idx)
{
    if (idx < 0 || idx >= pcaAppearanceParam().size()) {
        return NULL;
    }

    return getPcaAppearanceParam()[idx];
}

Node::FFTPCAFilter * AOGrammar::findCachedFFTPCAFilter(int idx)
{
    if (idx < 0 || idx >= cachedFFTPCAFilters().size()) {
        return NULL;
    }

    return getCachedFFTPCAFilters()[idx];
}
#endif

template<class Archive>
void AOGrammar::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(gType_);
    ar & BOOST_SERIALIZATION_NVP(name_);
    ar & BOOST_SERIALIZATION_NVP(note_);
    ar & BOOST_SERIALIZATION_NVP(year_);
    ar & BOOST_SERIALIZATION_NVP(isSingleObjModel_);

    ar.register_type(static_cast<Node *>(NULL));
    ar.register_type(static_cast<Edge *>(NULL));
    ar.register_type(static_cast<ParamUtil *>(NULL));
    ar.register_type(static_cast<Appearance *>(NULL));
    ar.register_type(static_cast<Offset *>(NULL));
    ar.register_type(static_cast<Deformation *>(NULL));
    ar.register_type(static_cast<Scaleprior *>(NULL));

    ar.template register_type<Node>();
    ar.template register_type<Edge>();
    ar.template register_type<ParamUtil>();
    ar.template register_type<Appearance>();
    ar.template register_type<Offset>();
    ar.template register_type<Deformation>();
    ar.template register_type<Scaleprior>();

    ar & BOOST_SERIALIZATION_NVP(nodeSet_);
    ar & BOOST_SERIALIZATION_NVP(edgeSet_);
    ar & BOOST_SERIALIZATION_NVP(isLRFlip_);

    ar & BOOST_SERIALIZATION_NVP(regMethod_);
    ar & BOOST_SERIALIZATION_NVP(appearanceSet_);
    ar & BOOST_SERIALIZATION_NVP(biasSet_);
    ar & BOOST_SERIALIZATION_NVP(deformationSet_);
    ar & BOOST_SERIALIZATION_NVP(scalepriorSet_);

    ar & BOOST_SERIALIZATION_NVP(maxDetectWindow_);
    ar & BOOST_SERIALIZATION_NVP(minDetectWindow_);

    ar & BOOST_SERIALIZATION_NVP(featureExtraction_);

    ar & BOOST_SERIALIZATION_NVP(idxRootNode_);

    ar & BOOST_SERIALIZATION_NVP(bboxPred_);
}

INSTANTIATE_BOOST_SERIALIZATION(AOGrammar);



void AOGrammar::save(const std::string & modelFile, int archiveType)
{
    std::ofstream out;
    out.open(modelFile.c_str(), std::ios::out);

    if ( !out.is_open() ) {
        RGM_LOG(error, ("Failed to write to file " + modelFile) );
        return;
    }

    finalize(false);

    switch (archiveType) {
    case 1: {
        boost::archive::text_oarchive oa(out);
        oa << *this;
        break;
    }
    case 0:
    default: {
        boost::archive::binary_oarchive oa(out);
        oa << *this;
        break;
    }

    } // switch

    out.close();
}

bool AOGrammar::read(const std::string & modelFile, int archiveType)
{
    std::ifstream in;
    in.open(modelFile.c_str(), std::ios::in);

    if ( !in.is_open() ) {
        //RGM_LOG(error, ("Failed to read file " + modelFile) );
        return false;
    }

    clear();

    switch (archiveType) {
    case 1: {
        boost::archive::text_iarchive ia(in);
        ia >> *this;
        break;
    }
    case 0:
    default: {
        boost::archive::binary_iarchive ia(in);
        ia >> *this;
        break;
    }

    } // switch

    in.close();

    finalize(true);

    return true;
}

void AOGrammar::visualize(const std::string & saveDir)
{
    std::string modelName = name() + "_" + year();
    std::string strDir = saveDir + FileUtil::FILESEP + modelName + "_vis"+FileUtil::FILESEP;
    FileUtil::VerifyDirectoryExists(strDir);

    std::string dotfile = strDir+modelName+".dot";
    FILE *f = fopen(dotfile.c_str(), "w");
    if ( f == NULL ) {
        RGM_LOG(error, ("Can not write file " + dotfile) );
        return;
    }

    // visualize all T-nodes
    pictureTNodes(strDir);

    // deformations
    pictureDeformation(strDir);

    // write .dot file for graphviz
    fprintf(f, "digraph %s {\n label=\"%s\";\n", modelName.c_str(), modelName.c_str());

    fprintf(f, "pack=true;\n overlap=false;\n labelloc=t;\n center=true;\n");

    // write the legend
    /*float wd = 2.0F;

    fprintf(f,  "subgraph legend {\n");

    //fprintf(f, "style=filled; color=lightgrey;\n");

    fprintf(f, "ORNode [shape=ellipse, style=bold, color=green, label=\"\"];\n");
    fprintf(f, "OR [shape=plaintext, style=solid, label=\"%s\"\r, width=%.1f];\n", "OR-node", wd);

    fprintf(f, "ANDNode [shape=ellipse, style=filled, color=blue, label=\"\"];\n");
    fprintf(f, "AND [shape=plaintext, style=solid, label=\"%s\"\r, width=%.1f];\n", "AND-node", wd);

    fprintf(f, "TNode [shape=box, style=bold, color=red, label=\"\"];\n");
    fprintf(f, "T [shape=plaintext, style=solid, label=\"%s\"\r, width=%.1f];\n", "TERMINAL-node", wd);

    fprintf(f, "sFromNode [shape=ellipse, style=bold, color=green, label=\"\"];\n");
    fprintf(f, "sToNode [shape=ellipse, style=filled, color=blue, label=\"\"];\n");
    fprintf(f, "edge [style=bold, color=green];\n");
    fprintf(f, "sFromNode -> sToNode;\n");
    fprintf(f, "Switching [shape=plaintext, style=solid, label=\"%s\"\r, width=%.1f];\n", "SWITCHING-edge", wd);

    fprintf(f, "cFromNode [shape=ellipse,  style=filled, color=blue, label=\"\"];\n");
    fprintf(f, "cToNode [shape=ellipse,style=bold, color=green, label=\"\"];\n");
    fprintf(f, "edge [style=bold, color=blue];\n");
    fprintf(f, "cFromNode -> cToNode;\n");
    fprintf(f, "Composition [shape=plaintext, style=solid, label=\"%s\"\r, width=%.1f];\n", "COMPOSITION-edge", wd);

    fprintf(f, "dFromNode [shape=ellipse,  style=filled, color=blue, label=\"\"];\n");
    fprintf(f, "dToNode [shape=box,style=bold, color=red, label=\"\"];\n");
    fprintf(f, "edge [style=bold, color=red];\n");
    fprintf(f, "dFromNode -> dToNode;\n");
    fprintf(f, "Deformation [shape=plaintext, style=solid, label=\"%s\"\r, width=%.1f];\n", "DEFORMATION-edge", wd);

    fprintf(f, "{ rank=source; rankdir=LR; OR AND T}\n");
    fprintf(f, "{ rank=source; rankdir=LR; Switching Composition Deformation}\n");

    fprintf( f, "};\n"); */

    const std::string imgExt(".png");
    const std::string textlabel("\"\"");
    const int bs = 20;
    const int padding = 2;
    const int zeroPaddingNum = 5;

    std::string saveName;

    // Write nodes using DFS
    for ( int i = 0; i < nodeDFS().size(); ++i ) {
        const Node * curNode = nodeDFS()[i];
        Node::nodeType t = curNode->type();
        int nodeIdx = idxNode(curNode);

        int flip = static_cast<int>(curNode->isLRFlip());
        std::string strFlip = NumToString_<int>(flip);

        switch (t) {
        case Node::T_NODE: {
            int appIdx = idxAppearance(curNode->appearance());
            std::string strApp = NumToString_<int>(appIdx, zeroPaddingNum);

            fprintf(f, "node%d [shape=box, style=bold, color=red, label=\"%s_%s\", labelloc=b, image=\"%sAppTemplate_%s_%s%s\"];\n",
                    nodeIdx, strApp.c_str(), strFlip.c_str(), strDir.c_str(), strApp.c_str(), strFlip.c_str(), imgExt.c_str());
            break;
        }
        case Node::AND_NODE: {
            if ( curNode->deformation() != NULL ) {
                int defIdx = idxDeformation(curNode->deformation());
                std::string strDef = NumToString_<int>(defIdx, zeroPaddingNum);

                fprintf(f, "node%d [shape=ellipse, style=filled, color=blue, label=%s, image=\"%sDeformation_%s_%s%s\"];\n",
                        nodeIdx, textlabel.c_str(), strDir.c_str(), strDef.c_str(), strFlip.c_str(), imgExt.c_str());
            } else if (curNode->offset() != NULL) {

                cv::Mat rootApp;
                std::vector<std::pair<cv::Mat, Node::Anchor> > parts;

                for ( int j = 0; j < curNode->outEdges().size(); ++j ) {
                    const Node * to = curNode->outEdges()[j]->toNode();

                    const Node * T = to;
                    while ( T->type() != Node::T_NODE) {
                        T = T->outEdges()[0]->toNode();
                    }

                    int appIdx = idxAppearance(T->appearance());
                    std::string strApp = NumToString_<int>(appIdx, zeroPaddingNum);
                    std::string imgFile = strDir + "AppTemplate_" + strApp + "_" + strFlip + imgExt;

                    if ( to->anchor()(2) == 0 ) {
                        rootApp = cv::imread(imgFile, cv::IMREAD_UNCHANGED);
                        assert(!rootApp.empty());
                    } else {
                        parts.push_back( std::make_pair(cv::imread(imgFile, cv::IMREAD_UNCHANGED), to->anchor()) );
                    }
                } // for j

                if ( parts.size() > 0 ) {
                    cv::Mat rootxApp(rootApp.rows * 2, rootApp.cols*2, rootApp.type());
                    cv::resize(rootApp, rootxApp, rootxApp.size(), 0, 0, cv::INTER_CUBIC);

                    for ( int j = 0; j < parts.size(); ++j ) {
                        int x = parts[j].second(0) * bs + padding;
                        int y = parts[j].second(1) * bs + padding;
                        parts[j].first.copyTo(rootxApp(cv::Rect(x, y, parts[j].first.cols, parts[j].first.rows)));

                        /*cv::imshow("debug", rootxApp);
                    cv::waitKey(0);*/
                    }

                    int nodeIdx = idxNode(curNode);
                    std::string strNode = NumToString_<int>(nodeIdx, zeroPaddingNum);

                    std::string saveName = strDir + "Component_" + strNode + "_" + strFlip  + imgExt;
                    cv::imwrite(saveName, rootxApp);

                    fprintf(f, "node%d [shape=ellipse, style=filled, color=blue, label=%s, image=\"%sComponent_%s_%s%s\"];\n",
                            nodeIdx, textlabel.c_str(), strDir.c_str(), strNode.c_str(), strFlip.c_str(), imgExt.c_str());
                } else {
                    fprintf(f, "node%d [shape=ellipse, style=filled, color=blue, label=%s];\n",
                            nodeIdx, textlabel.c_str());
                }
            } else {
                fprintf(f, "node%d [shape=ellipse, style=filled, color=blue, label=%s];\n",
                        nodeIdx, textlabel.c_str());
            }
            break;
        }
        case Node::OR_NODE: {
            fprintf(f, "node%d [shape=ellipse, style=bold, color=green, label=%s];\n",
                    nodeIdx, textlabel.c_str());
            break;
        }
        } // switch
    } // for i

    //fprintf(f, "nodeOR [shape=ellipse, style=bold, color=green, label=\"OR-node\"];\n");

    // Write edges using BFS
    for ( int i = 0; i < nodeBFS().size(); ++i ) {
        const Node * fromNode = nodeBFS()[i];
        int idxFrom = idxNode(fromNode);

        const std::vector<Edge *> & outEdges = fromNode->outEdges();
        for ( int j = 0; j < outEdges.size(); ++j ) {
            const Edge * curEdge = outEdges[j];
            Edge::edgeType t = curEdge->type();
            const Node * toNode = curEdge->toNode();
            int idxTo = idxNode(toNode);

            switch ( t ) {
            case Edge::SWITCHING: {
                fprintf(f, "edge [style=bold, color=green];\n");
                fprintf(f, "node%d -> node%d;\n", idxFrom, idxTo);
                break;
            }
            case Edge::COMPOSITION: {
                fprintf(f, "edge [style=bold, color=blue];\n");
                fprintf(f, "node%d -> node%d;\n", idxFrom, idxTo);
                break;
            }
            case Edge::DEFORMATION: {
                fprintf(f, "edge [style=bold, color=red];\n");
                fprintf(f, "node%d -> node%d;\n", idxFrom, idxTo);
                break;
            }
            case Edge::TERMINATION: {
                fprintf(f, "edge [style=bold, color=black];\n");
                fprintf(f, "node%d -> node%d;\n", idxFrom, idxTo);
                break;
            }
            }
        } // for j
    } // for i

    fprintf(f, "}");

    fclose(f);

    /// Use GraphViz
    std::string cmd = "dot -Tpdf " + dotfile + " -o " + strDir + modelName + ".pdf";
    std::system(cmd.c_str());

    cmd = "dot -Tpng " + dotfile + " -o " + strDir + modelName + ".png";
    std::system(cmd.c_str());
}

void AOGrammar::pictureTNodes(const std::string & saveDir)
{
    const std::string imgExt(".png");

    const int bs = 20;
    const int padding = 2;
    const int zeroPaddingNum = 5;

    for ( int i = 0; i < nodeSet().size(); ++i ) {
        const Node * curNode = nodeSet()[i];
        if (curNode->type() != Node::T_NODE) {
            continue;
        }

        int appIdx = idxAppearance(curNode->appearance());
        std::string strApp = NumToString_<int>(appIdx, zeroPaddingNum);

        int flip = static_cast<int>(curNode->isLRFlip());
        std::string strFlip = NumToString_<int>(flip);

        // Return the contrast insensitive orientations
        cv::Mat_<Scalar> wFolded = FeaturePyramid::fold( curNode->appearanceParam() );

        cv::Mat img = OpencvUtil::pictureHOG(wFolded, bs);

        cv::Mat imgPadding(img.rows+2*padding, img.cols+2*padding, img.type(), cv::Scalar::all(128));
        img.copyTo(imgPadding(cv::Rect(padding, padding, img.cols, img.rows)));

        cv::Mat imgShow;
        imgPadding.convertTo(imgShow, CV_8UC1);
        //cv::normalize(img, imgShow(cv::Rect(padding, padding, img.cols, img.rows)), 255, 0.0, CV_MINMAX, CV_8UC1);

        cv::imshow("AppTemplate", imgShow);
        cv::waitKey(2);

        std::string saveName = saveDir + std::string("AppTemplate_") + strApp + "_" + strFlip  + imgExt;
        cv::imwrite(saveName, imgShow);
    } // for i

    cv::destroyWindow("AppTemplate");
}

void AOGrammar::pictureDeformation(const std::string & saveDir)
{
    const std::string imgExt(".png");

    const int bs = 20;
    const int padding = 2;
    const int zeroPaddingNum = 5;
    const Scalar defScale = 500;

    for ( int i = 0; i < nodeSet().size(); ++i ) {
        const Node * curNode = nodeSet()[i];
        if (curNode->deformation() == NULL) {
            continue;
        }

        int defIdx = idxDeformation(curNode->deformation());
        std::string strDef = NumToString_<int>(defIdx, zeroPaddingNum);

        int flip = static_cast<int>(curNode->isLRFlip());
        std::string strFlip = NumToString_<int>(flip);

        const Deformation::Param w = curNode->deformationParam();

        int partHt = curNode->detectWindow().height() * bs;
        int partWd = curNode->detectWindow().width() * bs;

        cv::Mat_<Scalar> def(partHt, partWd, Scalar(0));

        int probex = partWd/2;
        int probey = partHt/2;

        Deformation::Param displacement;

        for ( int y=0; y<partHt; ++y ) {
            Scalar py = Scalar(probey - y) / bs;
            displacement(2) = py*py;
            displacement(3) = py;

            for ( int x=0; x<partWd; ++x ) {
                Scalar px = Scalar(probex - x) / bs;
                displacement(0) = px*px;
                displacement(1) = px;

                Scalar penalty = w.dot(displacement) * defScale;

                def(y, x) = penalty;
            }
        }

        cv::Mat imgPadding(def.rows+2*padding, def.cols+2*padding, def.type(), cv::Scalar::all(128));
        def.copyTo(imgPadding(cv::Rect(padding, padding, def.cols, def.rows)));

        cv::Mat imgShow;
        imgPadding.convertTo(imgShow, CV_8UC1);
        //cv::normalize(img, imgShow(cv::Rect(padding, padding, img.cols, img.rows)), 255, 0.0, CV_MINMAX, CV_8UC1);

        cv::imshow("Deformation", imgShow);
        cv::waitKey(2);

        std::string saveName = saveDir + std::string("Deformation_") + strDef + "_" + strFlip  + imgExt;
        cv::imwrite(saveName, imgShow);

    } // for i
}


Node * AOGrammar::addNode(Node::nodeType t)
{
    getNodeSet().push_back(new Node(t));
    return getNodeSet().back();
}

Edge * AOGrammar::addEdge(Node * from, Node * to, Edge::edgeType t)
{
    RGM_CHECK_NOTNULL(from);
    RGM_CHECK_NOTNULL(to);

    bool dup = false;
    for ( int i = 0; i < edgeSet().size(); ++i ) {
        if ( edgeSet()[i]->fromNode() == from && edgeSet()[i]->toNode() == to ) {
            dup = true;
            break;
        }
    }

    if ( dup ) {
        RGM_LOG(error, "add duplicated edge");
    }

    getEdgeSet().push_back(new Edge(t, from, to));
    Edge * e = getEdgeSet().back();

    from->getOutEdges().push_back(e);
    to->getInEdges().push_back(e);

    return e;
}

Offset * AOGrammar::addOffset(const Offset & bias)
{
    getBiasSet().push_back(new Offset(bias));
    return getBiasSet().back();
}

Scaleprior * AOGrammar::addScaleprior(const Scaleprior & prior)
{
    getScalepriorSet().push_back(new Scaleprior(prior));
    return getScalepriorSet().back();
}

Deformation * AOGrammar::addDeformation(const Deformation & def)
{
    getDeformationSet().push_back(new Deformation(def));
    return getDeformationSet().back();
}

std::pair<Node *, Edge *> AOGrammar::addChild(Node * parent, Node::nodeType chType, Edge::edgeType edgeType)
{
    RGM_CHECK_NOTNULL(parent);

    Node * ch = addNode(chType);

    Edge * e = addEdge(parent, ch, edgeType);

    return std::make_pair(ch, e);
}

std::pair<Node *, Edge *> AOGrammar::addParent(Node * ch, Node::nodeType paType, Edge::edgeType edgeType)
{
    RGM_CHECK_NOTNULL(ch);

    Node * pa = addNode(paType);

    Edge * e = addEdge(pa, ch, edgeType);

    return std::make_pair(pa, e);
}

Node * AOGrammar::addLRMirror(Node * root)
{
    RGM_CHECK_NOTNULL(root);

    // using BFS starting from the input root
    std::vector<Node *> bfs;
    std::vector<Node *> mirror;

    bfs.push_back(root);

    Node * mroot = addNode(root->type());
    mirror.push_back( mroot );

    setNodeLRFlip(root, mroot);

    int head = 0;
    while ( head < bfs.size() ) {
        Node * n  = bfs[head];
        Node * mn = mirror[head];
        head++;

        mn->getDetectWindow() = n->detectWindow();

        switch ( n->type() ) {
        case Node::OR_NODE: {
            for ( int i = 0; i < n->outEdges().size(); ++i ) {
                std::pair<Node *, Edge *> ch = addChild(mn, n->outEdges()[i]->toNode()->type(), n->outEdges()[i]->type());
                setNodeLRFlip(n->getOutEdges()[i]->getToNode(), ch.first);
                setEdgeLRFlip(n->getOutEdges()[i], ch.second);

                bfs.push_back(n->getOutEdges()[i]->getToNode());
                mirror.push_back(ch.first);
            }

            break;
        }
        case Node::AND_NODE: {
            mn->getOffset() = n->getOffset();
            mn->getDeformation() = n->getDeformation();
            mn->getScaleprior() = n->getScaleprior();

            for ( int i = 0; i < n->outEdges().size(); ++i ) {
                std::pair<Node *, Edge *> ch = addChild(mn, n->outEdges()[i]->toNode()->type(), n->outEdges()[i]->type());
                setNodeLRFlip(n->getOutEdges()[i]->getToNode(), ch.first);
                setEdgeLRFlip(n->getOutEdges()[i], ch.second);

                bfs.push_back(n->getOutEdges()[i]->getToNode());
                mirror.push_back(ch.first);
            }

            break;
        }
        case Node::T_NODE: {
            mn->getAppearance() = n->getAppearance();

            break;
        }
        } // switch
    } // while

    return getNodeSet()[idxNode(mroot)];
}

void AOGrammar::setNodeLRFlip(Node * n, Node * nFlip)
{
    RGM_CHECK_NOTNULL(n);
    RGM_CHECK_NOTNULL(nFlip);

    n->getIsLRFlip() = false;
    n->getLRMirrorNode() = nFlip;

    nFlip->getIsLRFlip() = true;
    nFlip->getLRMirrorNode() = n;
}

void AOGrammar::setEdgeLRFlip(Edge * e, Edge * eFlip)
{
    RGM_CHECK_NOTNULL(e);
    RGM_CHECK_NOTNULL(eFlip);

    e->getIsLRFlip() = false;
    e->getLRMirrorEdge() = eFlip;

    eFlip->getIsLRFlip() = true;
    eFlip->getLRMirrorEdge() = e;
}

void AOGrammar::finalize(bool hasIdx)
{
    if ( hasIdx ) {
        for ( int i = 0; i < nodeSet().size(); ++i ) {
            getNodeSet()[i]->assignConnections(this);
        }

        for ( int i = 0; i < edgeSet().size(); ++i ) {
            getEdgeSet()[i]->assignConnections(this);
        }

        if (rootNode() == NULL) {
            getRootNode() = findNode(idxRootNode_);
        }

    } else {
        for ( int i = 0; i < nodeSet().size(); ++i ) {
            getNodeSet()[i]->assignIdx(this);
        }

        for ( int i = 0; i < edgeSet().size(); ++i ) {
            getEdgeSet()[i]->assignIdx(this);
        }

        idxRootNode_ = idxNode(rootNode());
    }

    std::vector<int> visitedDFS(nodeSet().size(), 0);
    getNodeDFS().clear();
    traceNodeDFS(getRootNode(), visitedDFS, getNodeDFS());

    std::vector<int> visitedBFS(nodeSet().size(), 0);
    getNodeBFS().clear();
    traceNodeBFS(getRootNode(), visitedBFS, getNodeBFS());

    traceCompNodeDFSandBFS();
}

void save(const std::string & modelFile, const std::vector<AOGrammar> & models, int archiveType)
{
    std::ofstream out;
    out.open(modelFile.c_str(), std::ios::out);

    DEFINE_RGM_LOGGER;

    if ( !out.is_open() ) {
        RGM_LOG(error, ("Can not write models to " + modelFile));
        return;
    }


    switch (archiveType) {
    case 1: {
        boost::archive::text_oarchive oa(out);
        oa << models;
        break;
    }
    case 0:
    default: {
        boost::archive::binary_oarchive oa(out);
        oa << models;
        break;
    }
    } // switch

    out.close();

}

bool load(const std::string & modelFile, std::vector<AOGrammar> & models, int archiveType)
{
    std::ifstream in;
    in.open(modelFile.c_str(), std::ios::in);

    if ( !in.is_open() ) {
        return false;
    }

    models.clear();

    switch (archiveType) {
    case 1: {
        boost::archive::text_iarchive ia(in);
        ia >> models;
        break;
    }
    case 0:
    default: {
        boost::archive::binary_iarchive ia(in);
        ia >> models;
        break;
    }
    } // switch

    in.close();

    for ( int i = 0; i < models.size(); ++i ) {
        models[i].finalize(true);
    }

    return true;
}

#ifdef RGM_USE_MATLAB
namespace detail
{
/// Struct of a parameter block in voc-r5 matlab model file
struct matBlock
{
    char  t;
    double * pW;
    double * pLB;
    int shapeDim; // rows, cols, dim
    double * pShape;
    int len;
    Scalar learn;
    Scalar regCost;
};  // struct matBlock

void readBlock(matvar_t *blocks, int b,  matBlock &block)
{
    block.t      = *reinterpret_cast<char *>(Mat_VarGetStructFieldByName(blocks, "type", b)->data);
    block.pW     = reinterpret_cast<double *>(Mat_VarGetStructFieldByName(blocks, "w", b)->data);
    block.pLB    = reinterpret_cast<double *>(Mat_VarGetStructFieldByName(blocks, "lb", b)->data);
    block.pShape = reinterpret_cast<double *>(Mat_VarGetStructFieldByName(blocks, "shape", b)->data); // rows, cols, dim
    block.shapeDim = (int)block.pShape[1];
    block.len    = (int)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(blocks, "dim", b)->data);
    block.learn    = (Scalar)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(blocks, "learn", b)->data);
    block.regCost  = (Scalar)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(blocks, "reg_mult", b)->data);
}

} // namespace detail

void AOGrammar::cvtVOCRel5Model(const std::string &modelFile)
{
    using namespace detail;

    clear();

    mat_t *matfp = Mat_Open(modelFile.c_str(), MAT_ACC_RDONLY);
    RGM_CHECK_NOTNULL(matfp);

    matvar_t *model = Mat_VarRead(matfp, "model");
    RGM_CHECK_NOTNULL(model);

    // assume the DFS of the symbols has been sorted (using model_sort function in voc-r5)
    matvar_t *DFS = Mat_VarGetStructFieldByName(model, "DFS", 0);
    RGM_CHECK_NOTNULL(DFS);

    int numDFS = std::max<int>(DFS->dims[0], DFS->dims[1]);
    const double *pDFS = reinterpret_cast<const double *>( DFS->data );

    // Grammar type
    matvar_t *gType = Mat_VarGetStructFieldByName(model, "type", 0);
    char t;
    if ( gType == NULL ) {
        t = 'M';
    } else {
        t = *reinterpret_cast<char *>(gType->data);
    }

    switch (t) {
    case 'M':
        getType() = STARMIXTURE;
        break;
    case  'G':
        getType() = GRAMMAR;
        break;
    default:
        RGM_LOG(error, " Wrong grammar type.");
        return;
    }

    // name, note, year
    char * tmp = reinterpret_cast<char *>(Mat_VarGetStructFieldByName(model, "class", 0)->data);
    getName() = tmp==NULL ? "UNKNOWN" : std::string(tmp);
    tmp = reinterpret_cast<char *>(Mat_VarGetStructFieldByName(model, "note", 0)->data);
    getNote()= tmp==NULL ? "UNKNOWN" : std::string(tmp);
    tmp = reinterpret_cast<char *>(Mat_VarGetStructFieldByName(model, "year", 0)->data);
    getYear() = tmp==NULL ? "UNKNOWN" : std::string(tmp);

    getIsSingleObjModel() = true;

    // parameters
    matvar_t *blocks = Mat_VarGetStructFieldByName(model, "blocks", 0);

    // filters
    int numFilters = (int)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(model, "numfilters", 0)->data);
    matvar_t *filters = Mat_VarGetStructFieldByName(model, "filters", 0);

    // symbols
    matvar_t *symbols = Mat_VarGetStructFieldByName(model, "symbols", 0);

    int startSymbolIdx = *reinterpret_cast<double *>(Mat_VarGetStructFieldByName(model, "start", 0)->data) - 1;

    // rules
    matvar_t *rules = Mat_VarGetStructFieldByName(model, "rules", 0);

    // Correspondence
    std::vector<std::pair<int, Appearance *> >  correspFilter;
    std::vector<std::pair<int, Offset *> >  correspOffset;
    std::vector<std::pair<int, Deformation *> >  correspDef;
    std::vector<std::pair<int, Scaleprior *> >  correspScaleprior;

    matBlock block;

    // Add all appearance
    for ( int i = 0; i < numFilters; ++i ) {
        bool flip     = (bool)( *reinterpret_cast<char *>(Mat_VarGetStructFieldByName(filters, "flip", i)->data) );
        if ( flip )
            continue;

        int b  = (int)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(filters, "blocklabel", i)->data) - 1;

        readBlock(blocks, b, block);
        RGM_CHECK_EQ(block.t, 'F');
        RGM_CHECK_EQ(block.pShape[2], 32);

        getAppearanceSet().push_back( new Appearance() );

        getAppearanceSet().back()->getW() = Appearance::Param::Constant(block.pShape[0], block.pShape[1], Appearance::Cell::Zero());
        getAppearanceSet().back()->getLowerBound() = getAppearanceSet().back()->w();

        for ( int d = 0; d < block.pShape[2]; ++d ) {
            for ( int c = 0; c < block.pShape[1]; ++c) {
                for ( int r = 0; r < block.pShape[0]; ++r ) {
                    getAppearanceSet().back()->getW()(r, c)(d) = *(block.pW++);
                    getAppearanceSet().back()->getLowerBound()(r, c)(d) = *(block.pLB++);
                }
            }
        }

        getAppearanceSet().back()->getLearningStatus() = block.learn;
        getAppearanceSet().back()->getRegularizationCost() = block.regCost;

        correspFilter.push_back( std::pair<int, Appearance *>(b, getAppearanceSet().back()) );
    } // for i

    // convert to AOGrammar using DFS
    std::vector<std::pair<int, Node *> > symToNode;

    for ( int i = 0; i < numDFS; ++i ) {
        int symIdx = pDFS[i] - 1;
        matvar_t * rule = 0;
        rule = Mat_VarGetCell(rules, symIdx);
        if ( rule->class_type == MAT_C_EMPTY ) {
            // 'T' symbol
            t = *reinterpret_cast<char *>(Mat_VarGetStructFieldByName(symbols, "type", symIdx)->data);
            RGM_CHECK_EQ(t, 'T');

            Node * tNode = addNode(Node::T_NODE);

            // record it
            symToNode.push_back(std::pair<int, Node *>(symIdx, tNode));

            int f = (int)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(symbols, "filter", symIdx)->data) - 1;
            bool flip  = (bool)*reinterpret_cast<char *>(Mat_VarGetStructFieldByName(filters, "flip", f)->data);
            int b  = (int)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(filters, "blocklabel", f)->data) - 1;

            tNode->getIsLRFlip() = flip;

            for ( int k = 0; k < correspFilter.size(); ++k ) {
                if (correspFilter[k].first == b) {
                    tNode->getAppearance() = correspFilter[k].second;
                    break;
                }
            } // for k

            if (tNode->appearance() == NULL ) {
                RGM_LOG(error, "Not found the appearance for a T-node.");
                return;
            }

            tNode->getDetectWindow() = Rectangle2i(tNode->appearance()->w().cols(), tNode->appearance()->w().rows());

            continue;
        }

        const int numR = std::max<int>(rule->dims[0], rule->dims[1]);

        // Add an OR-node as the root for a rule ('D' or 'S')
        Node * orNode = addNode(Node::OR_NODE);

        std::vector<Node *> sChildren;

        if (symIdx == startSymbolIdx) {
            getRootNode() = orNode;
        }

        int lhs0 = (int)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(rule, "lhs", 0)->data) - 1;

        for ( int j = 0; j < numR; ++j ) {
            t = *reinterpret_cast<char *>(Mat_VarGetStructFieldByName(rule, "type", j)->data);

            // lhs
            int lhs = (int)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(rule, "lhs", j)->data) - 1;
            RGM_CHECK_EQ(lhs, lhs0);

            // rhs
            matvar_t *mx_rhs = Mat_VarGetStructFieldByName(rule, "rhs", j);
            const double *rhs = reinterpret_cast<const double *>(mx_rhs->data);
            const int rhs_len = std::max<int>(mx_rhs->dims[0], mx_rhs->dims[1]);

            // check rhs
            int foundRhs = 0;
            std::vector<Node *> children;
            for ( int k = 0; k < rhs_len; ++k ) {
                int curRhs = rhs[k] - 1;
                for (int l = 0; l < symToNode.size(); ++l ) {
                    if ( symToNode[l].first == curRhs ) {
                        foundRhs++;
                        children.push_back(symToNode[l].second);
                        break;
                    }
                } // for l
            } // for k
            if (foundRhs != rhs_len) {
                RGM_LOG(error, "The rhs haven't been converted yet.");
                return;
            }

            if (t == 'D') {
                // deformation rule
                RGM_CHECK_EQ(rhs_len, 1);
                RGM_CHECK_EQ(children[0]->inEdges().size(), 0);

                // Add an And-node
                Node * andNode = addNode(Node::AND_NODE);

                /// Connect to children

                // DEFORMATION edge
                Edge * defEdge = addEdge(andNode, children[0], Edge::DEFORMATION);

                // get the deformation parameter
                matvar_t * mxDef = Mat_VarGetStructFieldByName(rule, "def", j);
                int b = (int)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(mxDef, "blocklabel", 0)->data) - 1;
                bool flip = (bool)*reinterpret_cast<char *>(Mat_VarGetStructFieldByName(mxDef, "flip", 0)->data);

                defEdge->getIsLRFlip() = flip;
                andNode->getIsLRFlip() = flip;

                // add def to andNode
                Deformation * def = NULL;
                for ( int k = 0; k < correspDef.size(); ++k ) {
                    if (correspDef[k].first == b) {
                        def = correspDef[k].second;
                        break;
                    }
                } // for k

                if ( def == NULL ) {
                    getDeformationSet().push_back(new Deformation());
                    def = getDeformationSet().back();

                    readBlock(blocks, b, block);
                    RGM_CHECK_EQ(block.len, 4);

                    def->getW() << block.pW[0], block.pW[1], block.pW[2], block.pW[3]; // dx^2 dx dy^2 dy
                    def->getLowerBound() << block.pLB[0], block.pLB[1], block.pLB[2], block.pLB[3];
                    def->getLearningStatus() = block.learn;
                    def->getRegularizationCost() = block.regCost;

                    // record it
                    correspDef.push_back( std::pair<int, Deformation *>(b, def) );
                }

                andNode->getDeformation() = def;
                andNode->getDetectWindow() = children[0]->detectWindow();

                sChildren.push_back(andNode);

            } else {
                // structural rule

                // Add an And-node
                Node * andNode = addNode(Node::AND_NODE);

                matvar_t * mxAnchors = Mat_VarGetStructFieldByName(rule, "anchor", j);

                // det & shift window
                const double * detWind = reinterpret_cast<const double *>(Mat_VarGetStructFieldByName(rule, "detwindow", j)->data);
                Rectangle2i  detectWindow(detWind[1], detWind[0]);

                const double * shfitWind = reinterpret_cast<const double *>(Mat_VarGetStructFieldByName(rule, "shiftwindow", j)->data);
                Rectangle2i  shiftWindow(shfitWind[1], shfitWind[0]);
                if (shiftWindow.area() > 0 ) {
                    int area = shiftWindow.area();
                    RGM_LOG(error, boost::format("shiftWind: area=%d") % area);
                }

                andNode->getDetectWindow() = detectWindow;

                // offset
                matvar_t * mxOffset = Mat_VarGetStructFieldByName(rule, "offset", j);
                int b = (int)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(mxOffset, "blocklabel", 0)->data) - 1;

                Offset * bias = NULL;
                for ( int k = 0; k < correspOffset.size(); ++k ) {
                    if (correspOffset[k].first == b) {
                        bias = correspOffset[k].second;
                        break;
                    }
                } // for k

                if (bias == NULL) {
                    readBlock(blocks, b, block);
                    RGM_CHECK_EQ(block.len, 1);

                    if (block.pW[0] != 0.0F ) {
                        getBiasSet().push_back(new Offset());
                        bias = getBiasSet().back();

                        bias->getW() = block.pW[0];
                        bias->getLowerBound() = block.pLB[0];
                        bias->getLearningStatus() = block.learn;
                        bias->getRegularizationCost() = block.regCost;

                        // record it
                        correspOffset.push_back(std::pair<int, Offset *>(b, bias));
                    }
                }

                andNode->getOffset() = bias;

                // loc
                matvar_t * mxLoc = Mat_VarGetStructFieldByName(rule, "loc", j);
                b = (int)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(mxLoc, "blocklabel", 0)->data) -1;

                Scaleprior * scaleprior = NULL;
                for ( int k = 0; k < correspScaleprior.size(); ++k ) {
                    if (correspScaleprior[k].first == b) {
                        scaleprior = correspScaleprior[k].second;
                        break;
                    }
                } // for k

                if (scaleprior == NULL) {
                    readBlock(blocks, b, block);
                    RGM_CHECK_EQ(block.len, 3);

                    if (block.pW[0] != 0.0F || block.pW[1] != 0.0F || block.pW[2] != 0.0F) {
                        getScalepriorSet().push_back(new Scaleprior());
                        scaleprior = getScalepriorSet().back();

                        scaleprior->getW() << block.pW[0], block.pW[1], block.pW[2];
                        scaleprior->getLowerBound() << block.pLB[0], block.pLB[1], block.pLB[2];
                        scaleprior->getLearningStatus() = block.learn;
                        scaleprior->getRegularizationCost() = block.regCost;

                        // record it
                        correspScaleprior.push_back(std::pair<int, Scaleprior*>(b, scaleprior));
                    }
                }

                andNode->getScaleprior() = scaleprior;

                // Connect to children

                for ( int k = 0; k < children.size(); ++k ) {
                    int curRhs = rhs[0] - 1;
                    t = *reinterpret_cast<char *>(Mat_VarGetStructFieldByName(symbols, "type", curRhs)->data);

                    if (t == 'T') {
                        Node * andNode1 = addNode(Node::AND_NODE);

                        Edge * tEdge = addEdge(andNode1, children[k], Edge::TERMINATION);

                        andNode1->getDetectWindow() = children[k]->detectWindow();
                        andNode1->getIsLRFlip() = children[k]->isLRFlip();

                        children[k] = andNode1;

                    } else {
                        if ( children[k]->inEdges().size() > 0 ) {
                            Node * ch = addNode(children[k]->type());

                            for ( int l = 0; l < children[k]->outEdges().size(); ++l ) {
                                Edge * chEdge = addEdge(ch, children[k]->getOutEdges()[l]->getToNode(), children[k]->outEdges()[l]->type());
                            } // for l

                            ch->getDetectWindow() = children[k]->detectWindow();

                            children[k] = ch;
                        }
                    }

                    const double *anchor = reinterpret_cast<const double *>(Mat_VarGetCell(mxAnchors, k)->data);

                    children[k]->getAnchor() << anchor[0], anchor[1], anchor[2];

                    Edge * cEdge = addEdge(andNode, children[k], Edge::COMPOSITION);

                } // for k

                if ( children.size() ) {
                    andNode->getIsLRFlip() = children[0]->isLRFlip();
                    if (andNode->detectWindow().area()==0)
                        andNode->getDetectWindow() = children[0]->detectWindow();
                }

                sChildren.push_back(andNode);

            } // if
        } // for j

        // Connect orNode to its children
        for ( int j = 0; j < sChildren.size(); ++j ) {
            Edge * sEdge = addEdge(orNode, sChildren[j], Edge::SWITCHING);
        }

        if ( orNode != getRootNode() ) {
            orNode->getDetectWindow() = sChildren[0]->detectWindow();
            orNode->getIsLRFlip() = sChildren[0]->isLRFlip();
        }

        // lhs
        symToNode.push_back(std::pair<int, Node *>(lhs0, orNode));

    } // for i

    matvar_t * mxFeature = Mat_VarGetStructFieldByName(model, "features", 0);
    int featDim     = (int)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(mxFeature, "dim", 0)->data);
    RGM_CHECK_EQ(featDim, FeaturePyramid::NbFeatures);

    getFeatureBias()    = (Scalar)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(mxFeature, "bias", 0)->data);
    getCellSize()       = (int)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(mxFeature, "sbin", 0)->data);
    getExtraOctave()    = (bool)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(mxFeature, "extra_octave", 0)->data);

    getInterval()       = (int)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(model, "interval", 0)->data);
    getThresh()         = (Scalar)*reinterpret_cast<double *>(Mat_VarGetStructFieldByName(model, "thresh", 0)->data);

    // max and min detect window
    const double * maxDetWind = reinterpret_cast<const double *>(Mat_VarGetStructFieldByName(model, "maxsize", 0)->data);
    getMaxDetectWindow() = Rectangle2i(maxDetWind[1], maxDetWind[0]);

    const double * minDetWind = reinterpret_cast<const double *>(Mat_VarGetStructFieldByName(model, "minsize", 0)->data);
    getMinDetectWindow() = Rectangle2i(minDetWind[1], minDetWind[0]); 

    // Finalize
    finalize(false);

    // bbox prediction model
    getBboxPred().clear();
    matvar_t *bboxpred = Mat_VarGetStructFieldByName(model, "bboxpred", 0);
    if (bboxpred != NULL) {
        int num = std::max<int>(bboxpred->dims[0], bboxpred->dims[1]);
        RGM_CHECK_EQ(num, rootNode()->outEdges().size());

        for (int i = 0; i < num; ++i) {
            matvar_t * curPred = Mat_VarGetCell(bboxpred, i);
            RGM_CHECK_NOTNULL(curPred);

            // x1
            matvar_t * x1 = Mat_VarGetStructFieldByName(curPred, "x1", 0);
            int dimX1 = std::max<int>(x1->dims[0], x1->dims[1]);

            int idx = idxNode(rootNode()->outEdges()[i]->toNode()); // object component And node idx
            RGM_CHECK_GE(idx, 0);
            getBboxPred().insert(std::make_pair(idx, Matrix::Zero(dimX1, 4)));
            Matrix &pred(getBboxPred()[idx]);

            const double *pX1 = reinterpret_cast<const double *>(x1->data);
            for (int j = 0; j < dimX1; ++j){
                pred(j, 0) = pX1[j];
            }
            // y1
            matvar_t * y1 = Mat_VarGetStructFieldByName(curPred, "y1", 0);
            int dimY1 = std::max<int>(y1->dims[0], y1->dims[1]);
            RGM_CHECK_EQ(dimY1, dimX1);
            const double *pY1 = reinterpret_cast<const double *>(y1->data);
            for (int j = 0; j < dimY1; ++j){
                pred(j, 1) = pY1[j];
            }
            // x2
            matvar_t * x2 = Mat_VarGetStructFieldByName(curPred, "x2", 0);
            int dimX2 = std::max<int>(x2->dims[0], x2->dims[1]);
            RGM_CHECK_EQ(dimX2, dimX1);
            const double *pX2 = reinterpret_cast<const double *>(x2->data);
            for (int j = 0; j < dimX2; ++j){
                pred(j, 2) = pX2[j];
            }
            // y2
            matvar_t * y2 = Mat_VarGetStructFieldByName(curPred, "y2", 0);
            int dimY2 = std::max<int>(y2->dims[0], y2->dims[1]);
            RGM_CHECK_EQ(dimY2, dimX1);
            const double *pY2 = reinterpret_cast<const double *>(y2->data);
            for (int j = 0; j < dimY2; ++j){
                pred(j, 3) = pY2[j];
            }
        }
    }

    Mat_VarFree(model);
    Mat_Close(matfp);
}

#endif

} // namespace RGM
