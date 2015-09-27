#include <algorithm>
#include <numeric>

#include "ParseTree.hpp"
#include "AOGrammar.hpp"


namespace RGM
{

// ------ PtEdge ------

PtEdge::PtEdge()
{
    idx_.fill(-1);
}

PtEdge::~PtEdge()
{
}

PtEdge::PtEdge(const PtEdge & e)
{
    idx_ = e.idx();
}

PtEdge::PtEdge(int fromNode, int toNode)
{
    idx_ << -1, fromNode, toNode, -1;
}

PtEdge::PtEdge(int fromNode, int toNode, int gEdge)
{
    idx_ << -1, fromNode, toNode, gEdge;
}

const PtEdge::Index & PtEdge::idx() const
{
    return idx_;
}

PtEdge::Index & PtEdge::getIdx()
{
    return idx_;
}

const PtNode * PtEdge::fromNode(const ParseTree &pt) const
{
    int i = idx()[IDX_FROM];
    if ( i < 0 || i >= pt.nodeSet().size() )
        return NULL;

    return pt.nodeSet()[i];
}

PtNode * PtEdge::getFromNode(ParseTree &pt)
{
    int i = idx()[IDX_FROM];
    if ( i < 0 || i >= pt.nodeSet().size() )
        return NULL;

    return pt.getNodeSet()[i];
}

const PtNode * PtEdge::toNode(const ParseTree &pt) const
{
    int i = idx()[IDX_TO];
    if ( i < 0 || i >= pt.nodeSet().size() )
        return NULL;

    return pt.nodeSet()[i];
}

PtNode * PtEdge::getToNode(ParseTree &pt)
{
    int i = idx()[IDX_TO];
    if ( i < 0 || i >= pt.nodeSet().size() )
        return NULL;

    return pt.getNodeSet()[i];
}

template<class Archive>
void PtEdge::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(idx_);
}

INSTANTIATE_BOOST_SERIALIZATION(PtEdge);



// ------ PtNode ------

PtNode::PtNode()
{
    idx_.fill(-1);
    idx_(IDX_VALID) = 1;
}

PtNode::~PtNode()
{
}

PtNode::PtNode(const PtNode & n)
{
    idxInEdges_ = n.idxInEdges();
    idxOutEdges_ = n.idxOutEdges();
    idx_ = n.idx();
}

PtNode::PtNode(int gNode)
{
    idx_.fill(-1);
    idx_(IDX_VALID) = 1;
    idx_[IDX_G] = gNode;
}

const std::vector<int> & PtNode::idxInEdges() const
{
    return idxInEdges_;
}

std::vector<int> & PtNode::getIdxInEdges()
{
    return idxInEdges_;
}

const std::vector<int> & PtNode::idxOutEdges() const
{
    return idxOutEdges_;
}

std::vector<int> & PtNode::getIdxOutEdges()
{
    return idxOutEdges_;
}

const PtEdge * PtNode::inEdge(int i, const ParseTree &pt) const
{
    if ( i < 0  || i >= idxInEdges().size() )
        return NULL;

    return pt.edgeSet()[idxInEdges()[i]];
}

PtEdge * PtNode::getInEdge(int i, ParseTree &pt)
{
    if ( i < 0  || i >= idxInEdges().size() )
        return NULL;

    return pt.getEdgeSet()[idxInEdges()[i]];
}

const PtEdge * PtNode::outEdge(int i, const ParseTree &pt) const
{
    if ( i < 0  || i >= idxOutEdges().size() )
        return NULL;

    return pt.edgeSet()[idxOutEdges()[i]];
}

PtEdge * PtNode::getOutEdge(int i, ParseTree &pt)
{
    if ( i < 0  || i >= idxOutEdges().size() )
        return NULL;

    return pt.getEdgeSet()[idxOutEdges()[i]];
}

const PtNode::Index & PtNode::idx() const
{
    return idx_;
}

PtNode::Index & PtNode::getIdx()
{
    return idx_;
}

const ParseInfo * PtNode::parseInfo(const ParseTree * pt) const
{
    DEFINE_RGM_LOGGER;

    RGM_CHECK_NOTNULL(pt);

    return pt->parseInfoSet()[idx()[IDX_PARSEINFO]];
}

ParseInfo *& PtNode::getParseInfo(ParseTree * pt)
{
    DEFINE_RGM_LOGGER;

    RGM_CHECK_NOTNULL(pt);

    return pt->getParseInfoSet()[idx()[IDX_PARSEINFO]];
}

template<class Archive>
void PtNode::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(idxInEdges_);
    ar & BOOST_SERIALIZATION_NVP(idxOutEdges_);
    ar & BOOST_SERIALIZATION_NVP(idx_);
}

INSTANTIATE_BOOST_SERIALIZATION(PtNode);




// ------ ParseTree::States ------

ParseTree::States::States() :
    isBelief_(false), score_(0), loss_(1), margin_(0), norm_(0)
{
}

ParseTree::States::States(const States & s) :
    isBelief_(s.isBelief_), score_(s.score_), loss_(s.loss_), margin_(s.margin_), norm_(s.norm_)
{
}

ParseTree::States::States(bool isBelief, Scalar loss, Scalar norm) :
    isBelief_(isBelief), score_(0), loss_(loss), margin_(0), norm_(norm)
{
}

ParseTree::States::States(bool isBelief, Scalar score, Scalar loss, Scalar norm) :
    isBelief_(isBelief), score_(score), loss_(loss), margin_(0), norm_(norm)
{
}

template<class Archive>
void ParseTree::States::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(isBelief_);
    ar & BOOST_SERIALIZATION_NVP(score_);
    ar & BOOST_SERIALIZATION_NVP(loss_);
    ar & BOOST_SERIALIZATION_NVP(margin_);
    ar & BOOST_SERIALIZATION_NVP(norm_);
}

INSTANTIATE_BOOST_SERIALIZATION(ParseTree::States);




// ------ ParseTree ------

ParseTree::ParseTree() :
    g_(NULL), idxRootNode_(-1), states_(NULL), dataId_(-1),
    appearanceX_(NULL), imgWd_(0), imgHt_(0)
{
}

ParseTree::ParseTree(const ParseTree & pt) :
    g_(NULL), idxRootNode_(-1), states_(NULL), dataId_(pt.dataId()),
    appearanceX_(NULL), imgWd_(pt.imgWd()), imgHt_(pt.imgHt())
{
    nodeSet_.resize(pt.nodeSet().size(), NULL);
    for ( int i = 0; i < nodeSet_.size(); ++i ) {
        nodeSet_[i] = new PtNode(*pt.nodeSet()[i]);
    }

    edgeSet_.resize(pt.edgeSet().size(), NULL);
    for ( int i = 0; i < edgeSet_.size(); ++i ) {
        edgeSet_[i] = new PtEdge(*pt.edgeSet()[i]);
    }

    idxRootNode_ = pt.idxRootNode();

    g_ = pt.grammar();

    appearanceSet_.resize(pt.appearanceSet().size(), NULL);
    for ( int i = 0; i < appearanceSet_.size(); ++i ) {
        appearanceSet_[i] = new Appearance::Param(*pt.appearanceSet()[i]);
    }

    biasSet_ = pt.biasSet();

    deformationSet_.resize(pt.deformationSet().size(), NULL);
    for ( int i = 0; i < deformationSet_.size(); ++i ) {
        deformationSet_[i] = new Deformation::Param(*pt.deformationSet()[i]);
    }

    scalepriorSet_.resize(pt.scalepriorSet().size(), NULL);
    for ( int i = 0; i < scalepriorSet_.size(); ++i ) {
        scalepriorSet_[i] = new Scaleprior::Param(*pt.scalepriorSet()[i]);
    }

    parseInfoSet_.resize(pt.parseInfoSet().size(), NULL);
    for ( int i = 0; i < parseInfoSet_.size(); ++i ) {
        parseInfoSet_[i] = new ParseInfo(*pt.parseInfoSet()[i]);
    }

    if (pt.states() != NULL ) {
        states_ = new States(*pt.states());
    }

    appUsage_ = pt.appUsage();

    if ( pt.appearaceX() != NULL ) {
        appearanceX_ = new Appearance::Param( *(pt.appearaceX()) );
    }

}

ParseTree::~ParseTree()
{
    clear();
}

ParseTree::ParseTree(const AOGrammar & g) :
    g_(&g), idxRootNode_(-1), states_(NULL)
{
}

void ParseTree::setGrammar(const AOGrammar & g)
{
    g_ = &g;
}

ParseTree & ParseTree::operator=(const ParseTree & pt)
{
    if (this == &pt) {
        return *this;
    }

    clear();

    nodeSet_.resize(pt.nodeSet().size(), NULL);
    for ( int i = 0; i < nodeSet_.size(); ++i ) {
        nodeSet_[i] = new PtNode(*pt.nodeSet()[i]);
    }

    edgeSet_.resize(pt.edgeSet().size(), NULL);
    for ( int i = 0; i < edgeSet_.size(); ++i ) {
        edgeSet_[i] = new PtEdge(*pt.edgeSet()[i]);
    }

    idxRootNode_ = pt.idxRootNode();

    g_ = pt.grammar();

    appearanceSet_.resize(pt.appearanceSet().size(), NULL);
    for ( int i = 0; i < appearanceSet_.size(); ++i ) {
        appearanceSet_[i] = new Appearance::Param(*pt.appearanceSet()[i]);
    }

    biasSet_ = pt.biasSet();

    deformationSet_.resize(pt.deformationSet().size(), NULL);
    for ( int i = 0; i < deformationSet_.size(); ++i ) {
        deformationSet_[i] = new Deformation::Param(*pt.deformationSet()[i]);
    }

    scalepriorSet_.resize(pt.scalepriorSet().size(), NULL);
    for ( int i = 0; i < scalepriorSet_.size(); ++i ) {
        scalepriorSet_[i] = new Scaleprior::Param(*pt.scalepriorSet()[i]);
    }

    parseInfoSet_.resize(pt.parseInfoSet().size(), NULL);
    for ( int i = 0; i < parseInfoSet_.size(); ++i ) {
        parseInfoSet_[i] = new ParseInfo(*pt.parseInfoSet()[i]);
    }

    if (pt.states() != NULL ) {
        states_ = new States(*pt.states());
    }

    appUsage_ = pt.appUsage();

    dataId_ = pt.dataId();

    if ( pt.appearaceX() != NULL ) {
        appearanceX_ = new Appearance::Param( *(pt.appearaceX()) );
    }

    imgWd_ = pt.imgWd();
    imgHt_ = pt.imgHt();

    return *this;
}

bool ParseTree::operator<(const ParseTree & pt) const
{
    return rootNode()->parseInfo(this)->score_ > pt.rootNode()->parseInfo(&pt)->score_; // for decreasing sort
}

void ParseTree::swap(ParseTree & pt)
{
    if ( this == &pt ) {
        return;
    }

    clear();

    nodeSet_.swap(pt.getNodeSet());

    edgeSet_.swap(pt.getEdgeSet());

    /*nodeSet_.resize(pt.nodeSet().size(), NULL);
    for ( int i = 0; i < nodeSet_.size(); ++i ) {
    	nodeSet_[i] = new PtNode(*pt.nodeSet()[i]);
    }

    edgeSet_.resize(pt.edgeSet().size(), NULL);
    for ( int i = 0; i < edgeSet_.size(); ++i ) {
    	edgeSet_[i] = new PtEdge(*pt.edgeSet()[i]);
    }*/

    idxRootNode_ = pt.idxRootNode();

    g_ = pt.grammar();

    appearanceSet_.swap(pt.getAppearanceSet());
    biasSet_.swap(pt.getBiasSet());
    deformationSet_.swap(pt.getDeformationSet());
    scalepriorSet_.swap(pt.getScalepriorSet());
    parseInfoSet_.swap(pt.getParseInfoSet());
    appUsage_.swap(pt.getAppUsage());

    /*appearanceSet_.resize(pt.appearanceSet().size(), NULL);
    for ( int i = 0; i < appearanceSet_.size(); ++i ) {
    	appearanceSet_[i] = new Appearance::Param(*pt.appearanceSet()[i]);
    }

    biasSet_ = pt.biasSet();

    deformationSet_.resize(pt.deformationSet().size(), NULL);
    for ( int i = 0; i < deformationSet_.size(); ++i ) {
    	deformationSet_[i] = new Deformation::Param(*pt.deformationSet()[i]);
    }

    scalepriorSet_.resize(pt.scalepriorSet().size(), NULL);
    for ( int i = 0; i < scalepriorSet_.size(); ++i ) {
    	scalepriorSet_[i] = new Scaleprior::Param(*pt.scalepriorSet()[i]);
    }

    parseInfoSet_.resize(pt.parseInfoSet().size(), NULL);
    for ( int i = 0; i < parseInfoSet_.size(); ++i ) {
    	parseInfoSet_[i] = new ParseInfo(*pt.parseInfoSet()[i]);
    }

    appUsage_ = pt.appUsage();*/

    if (pt.states() != NULL ) {
        //states_ = new States(*pt.states());
        std::swap(states_, pt.getStates());
    }

    dataId_ = pt.dataId();

    if ( pt.appearaceX() != NULL ) {
        std::swap(appearanceX_, pt.getAppearaceX());
    }

    imgWd_ = pt.imgWd();
    imgHt_ = pt.imgHt();
}

void ParseTree::clear()
{
    for ( int i = 0; i < nodeSet_.size(); ++i ) {
        delete nodeSet_[i];
    }
    nodeSet_.clear();

    for ( int i = 0; i < edgeSet_.size(); ++i ) {
        delete edgeSet_[i];
    }
    edgeSet_.clear();

    idxRootNode_ = -1;

    g_ = NULL;

    for ( int i = 0; i < appearanceSet_.size(); ++i ) {
        delete appearanceSet_[i];
    }
    appearanceSet_.clear();

    biasSet_.clear();

    for ( int i = 0; i < deformationSet_.size(); ++i ) {
        delete deformationSet_[i];
    }
    deformationSet_.clear();

    for ( int i = 0; i < scalepriorSet_.size(); ++i ) {
        delete scalepriorSet_[i];
    }
    scalepriorSet_.clear();

    for ( int i = 0; i < parseInfoSet_.size(); ++i ) {
        delete parseInfoSet_[i];
    }
    parseInfoSet_.clear();

    if ( states_ != NULL ) {
        delete states_;
    }

    dataId_ = -1;

    if ( appearanceX_ != NULL ) {
        delete appearanceX_;
    }

    imgWd_ = 0;
    imgHt_ = 0;

}

bool ParseTree::empty() const
{
    return nodeSet().size() == 0;
}

const std::vector<PtNode *> & ParseTree::nodeSet() const
{
    return nodeSet_;
}

std::vector<PtNode *> & ParseTree::getNodeSet()
{
    return nodeSet_;
}

const std::vector<PtEdge *> & ParseTree::edgeSet() const
{
    return edgeSet_;
}

std::vector<PtEdge *> & ParseTree::getEdgeSet()
{
    return edgeSet_;
}

int ParseTree::idxRootNode() const
{
    return idxRootNode_;
}

int &  ParseTree::getIdxRootNode()
{
    return idxRootNode_;
}

const PtNode * ParseTree::rootNode() const
{
    return nodeSet()[idxRootNode()];
}

PtNode *& ParseTree::getRootNode()
{
    return getNodeSet()[idxRootNode()];
}

const AOGrammar *  ParseTree::grammar() const
{
    return g_;
}

const std::vector<Appearance::Param *> & ParseTree::appearanceSet() const
{
    return appearanceSet_;
}

std::vector<Appearance::Param *> & ParseTree::getAppearanceSet()
{
    return appearanceSet_;
}

const std::vector<Scalar> &  ParseTree::biasSet() const
{
    return biasSet_;
}

std::vector<Scalar> & ParseTree::getBiasSet()
{
    return biasSet_;
}

const std::vector<Deformation::Param *> & ParseTree::deformationSet() const
{
    return deformationSet_;
}

std::vector<Deformation::Param *> & ParseTree::getDeformationSet()
{
    return deformationSet_;
}

const std::vector<Scaleprior::Param *> & ParseTree::scalepriorSet() const
{
    return scalepriorSet_;
}

std::vector<Scaleprior::Param *> & ParseTree::getScalepriorSet()
{
    return scalepriorSet_;
}

const std::vector<ParseInfo *> & ParseTree::parseInfoSet() const
{
    return parseInfoSet_;
}

std::vector<ParseInfo *> & ParseTree::getParseInfoSet()
{
    return parseInfoSet_;
}

const ParseTree::States * ParseTree::states() const
{
    return states_;
}

ParseTree::States *& ParseTree::getStates()
{
    return states_;
}

int  ParseTree::dataId() const
{
    return dataId_;
}

int & ParseTree::getDataId()
{
    return dataId_;
}

const std::vector<int> & ParseTree::appUsage() const
{
    return appUsage_;
}

std::vector<int> & ParseTree::getAppUsage()
{
    return appUsage_;
}

const Appearance::Param *  ParseTree::appearaceX() const
{
    return appearanceX_;
}

Appearance::Param *& ParseTree::getAppearaceX()
{
    return appearanceX_;
}

int ParseTree::imgWd() const
{
    return imgWd_;
}

int & ParseTree::getImgWd()
{
    return imgWd_;
}

int ParseTree::imgHt() const
{
    return imgHt_;
}

int & ParseTree::getImgHt()
{
    return imgHt_;
}

int ParseTree::idxObjComp() const
{
    assert(!empty());

    return rootNode()->parseInfo(this)->c_;
}

int ParseTree::addNode(int gNode, int type)
{
    int idx = nodeSet().size();

//    for ( int i = 0; i < idx; ++i ) {
//        if ( nodeSet()[i]->idx()[PtNode::IDX_G] == gNode ) {
//            RGM_LOG(error, "duplicated pt node" );
//            return -1;
//        }
//    }

    getNodeSet().push_back( new PtNode(gNode) );

    getNodeSet().back()->getIdx()[PtNode::IDX_MYSELF] = idx;
    getNodeSet().back()->getIdx()[PtNode::IDX_TYPE] = type;

    return idx;
}

int ParseTree::addEdge(int fromNode, int toNode, int gEdge)
{
    int idx = edgeSet().size();

    for ( int i = 0; i < idx; ++i ) {
        if ( edgeSet()[i]->idx()[PtEdge::IDX_FROM] == fromNode &&
                edgeSet()[i]->idx()[PtEdge::IDX_TO] == toNode ) {
            RGM_LOG(error, "duplicated pt edge" );
            return -1;
        }
    }

    getEdgeSet().push_back( new PtEdge(fromNode, toNode, gEdge));
    getEdgeSet().back()->getIdx()[PtEdge::IDX_MYSELF] = idx;

    getNodeSet()[fromNode]->getIdxOutEdges().push_back(idx);
    getNodeSet()[toNode]->getIdxInEdges().push_back(idx);

    return idx;
}

int ParseTree::AddBias(Scalar w)
{
    int idx = biasSet().size();

    getBiasSet().push_back(w);

    return idx;
}

int ParseTree::addScaleprior(Scaleprior::Param & w)
{
    int idx = scalepriorSet().size();

    getScalepriorSet().push_back(new Scaleprior::Param(w) );

    return idx;
}

int ParseTree::addDeformation(Scalar dx, Scalar dy, bool flip)
{
    int idx = deformationSet().size();

    Deformation::Param w;
    w << dx * dx, dx, dy * dy, dy;
    w *= -1.0F;

    if ( flip ) {
        w(1) *= -1;
    }

    getDeformationSet().push_back(new Deformation::Param(w) );

    return idx;
}

int ParseTree::addAppearance(Appearance::Param & w, bool flip)
{
    int idx = appearanceSet().size();

    getAppearanceSet().push_back( new Appearance::Param() );

    if ( flip ) {
        getAppearanceSet().back()->swap(FeaturePyramid::Flip( w ));
    } else {
        getAppearanceSet().back()->swap(w);
    }

    return idx;
}

int ParseTree::addParseInfo(ParseInfo & info)
{
    int idx = parseInfoSet().size();

    getParseInfoSet().push_back( new ParseInfo(info) );

    return idx;
}

void ParseTree::showDetection(cv::Mat img, bool display)
{
    for ( int i = 0; i < nodeSet().size(); ++i ) {
        if ( nodeSet()[i]->idx()(PtNode::IDX_TYPE) != static_cast<int>(Node::T_NODE) ) {
            continue;
        }
        const ParseInfo * info = nodeSet()[i]->parseInfo(this);
        cv::rectangle(img, info->cvRect(), cv::Scalar::all(255), 3);
        cv::rectangle(img, info->cvRect(), cv::Scalar(255, 0, 0), 2);
        /*if ( display ) {
        	cv::String winName("AOGDetection");
        	cv::imshow(winName, img);
        	cv::waitKey(0);
        }*/
    }

    cv::rectangle(img, rootNode()->parseInfo(this)->cvRect(), cv::Scalar::all(255), 5);
    cv::rectangle(img, rootNode()->parseInfo(this)->cvRect(), cv::Scalar(0, 0, 255), 3);

    if ( display ) {
        cv::String winName("AOGDetection");
        cv::imshow(winName, img);
        cv::waitKey(0);
    }
}

int ParseTree::dim() const
{
    int d = 0;

    for ( int i = 0; i < appearanceSet().size(); ++i ) {
        d += appearanceSet()[i]->size() * FeaturePyramid::NbFeatures;
    }

    d += biasSet().size();

    d += deformationSet().size() * 4;

    d += scalepriorSet().size() * 3;

    return d;
}

int ParseTree::compareFeatures(const ParseTree & pt) const
{
    assert(grammar() != NULL && grammar() == pt.grammar() );

    for ( int i = 0;  i < nodeSet().size(); ++i ) {
        int idxG = nodeSet()[i]->idx()(PtNode::IDX_G);

        int ii = 0;
        for ( ; ii < pt.nodeSet().size(); ++ii ) {
            int idxG1 = pt.nodeSet()[ii]->idx()(PtNode::IDX_G);
            if ( idxG1 == idxG ) {
                break;
            }
        } // for ii

        if ( ii == pt.nodeSet().size() ) {
            return 1;
        }

        for ( int j = PtNode::IDX_BIAS; j < PtNode::IDX_APP+1; ++j ) {
            int fidx  = nodeSet()[i]->idx()(j);
            int fidx1 = pt.nodeSet()[ii]->idx()(j);
            if (fidx != -1 && fidx1 != -1 ) {
                switch (j) {
                case PtNode::IDX_BIAS: {
                    if ( biasSet()[fidx] > pt.biasSet()[fidx1] ) {
                        return 1;
                    } else if ( biasSet()[fidx] < pt.biasSet()[fidx1] ) {
                        return -1;
                    }

                    break;
                }
                case PtNode::IDX_DEF: {
                    const Deformation::Param & p(*deformationSet()[fidx]);
                    const Deformation::Param & p1(*pt.deformationSet()[fidx1]);
                    for ( int k = 0; k < 4; ++k ) {
                        if ( p(k) > p1(k) ) {
                            return 1;
                        } else if (p(k) < p1(k)) {
                            return -1;
                        }
                    }
                    break;
                }
                case PtNode::IDX_SCALEPRIOR: {
                    const Scaleprior::Param & p(*scalepriorSet()[fidx]);
                    const Scaleprior::Param & p1(*pt.scalepriorSet()[fidx1]);
                    for ( int k = 0; k < 3; ++k ) {
                        if ( p(k) > p1(k) ) {
                            return 1;
                        } else if (p(k) < p1(k)) {
                            return -1;
                        }
                    }
                    break;
                }
                case PtNode::IDX_APP: {
                    const Appearance::Param & p(*appearanceSet()[fidx]);
                    const Appearance::Param & p1(*pt.appearanceSet()[fidx1]);
                    for ( int row=0; row < p.rows(); ++ row )
                        for ( int col=0; col < p.cols(); ++col )
                            for ( int k = 0; k < FeaturePyramid::NbFeatures; ++k ) {
                                if ( p(row, col)(k) > p1(row, col)(k) ) {
                                    return 1;
                                } else if (p(row, col)(k) < p1(row, col)(k)) {
                                    return -1;
                                }
                            }
                    break;
                }
                }
            } else {
                if (fidx != -1 && fidx1 == -1) {
                    return 1;
                }

                if (fidx == -1 && fidx1 != -1) {
                    return -1;
                }
            }
        }

    } // for i

    return 0;
}

Scalar ParseTree::norm() const
{
    Scalar n = 0;

    for ( int i = 0; i < appearanceSet().size(); ++i ) {
        n += FeaturePyramid::Map(*(appearanceSet()[i])).squaredNorm();
    }

    n += std::inner_product(biasSet().begin(), biasSet().end(), biasSet().begin(), 0);

    for ( int i = 0; i < deformationSet().size(); ++i ) {
        n += deformationSet()[i]->squaredNorm();
    }

    for ( int i = 0; i < scalepriorSet().size(); ++i ) {
        n += scalepriorSet()[i]->squaredNorm();
    }

    return std::sqrt(n);
}

Scalar ParseTree::computeOverlapLoss(const Rectangle2i & ref) const
{
    Intersector_<int> inter(ref, 0.5F, true);

    const ParseInfo * p = rootNode()->parseInfo(this);

    Rectangle2i box(p->x(), p->y(), p->width(), p->height());

    Scalar ov = 0;

    inter(box, &ov);

    return 1 - ov;
}

std::vector<const PtNode *> ParseTree::findNode(const Node *n)
{
    RGM_CHECK_NOTNULL(grammar());
    RGM_CHECK_NOTNULL(n);

    int gIdx = grammar()->idxNode(n);
    RGM_CHECK_NOTEQ(gIdx, -1);   

    return findNode(gIdx);
}

std::vector<const PtNode *> ParseTree::findNode(const int idxG)
{
    std::vector<const PtNode *>  n;

    for ( int i = 0; i < nodeSet().size(); ++i ) {
        if ( nodeSet()[i]->idx()[PtNode::IDX_G] == idxG ) {
            n.push_back(nodeSet()[i]);
        }
    }

    return n;
}

std::vector<PtNode *> ParseTree::getNode(const int idxG)
{
    std::vector<PtNode *> n;

    for ( int i = 0; i < nodeSet().size(); ++i ) {
        if ( nodeSet()[i]->idx()[PtNode::IDX_G] == idxG ) {
            n.push_back(getNodeSet()[i]);
        }
    }

    return n;
}


std::vector<const PtNode *> ParseTree::findSingleObjAndNodes() const
{
    const ParseTree &pt(*this);

    std::vector<const PtNode *> sobj;

    for ( int i = 0; i < nodeSet().size(); ++i ) {
        const PtNode * n = nodeSet()[i];
        if ( n->idx()[PtNode::IDX_TYPE] != static_cast<int>(Node::T_NODE) )
            continue;

        // Bottom-up: tnode -> and-node -> or-node -> single obj and-node
        const PtNode * a1 = n->inEdge(0, pt)->fromNode(pt);
        const PtNode * o  = a1->inEdge(0, pt)->fromNode(pt);
        const PtNode * a2 = o->inEdge(0, pt)->fromNode(pt);

        if ( a2->idx()[PtNode::IDX_VALID] <= 0 )
            continue;

        bool found = false;
        for ( int j = 0; j < sobj.size(); ++j ) {
            if ( a2 == sobj[j] ) {
                found = true;
                break;
            }
        }
        if ( !found ) {
            sobj.push_back(a2);
        }
    }

    return sobj;
}

std::vector<PtNode *> ParseTree::getSingleObjAndNodes()
{
    ParseTree &pt(*this);

    std::vector<PtNode *> sobj;

    for ( int i = 0; i < nodeSet().size(); ++i ) {
        PtNode * n = getNodeSet()[i];
        if ( n->idx()[PtNode::IDX_TYPE] != static_cast<int>(Node::T_NODE) )
            continue;

        // Bottom-up: tnode -> and-node -> or-node -> single obj and-node
        PtNode * a1 = n->getInEdge(0, pt)->getFromNode(pt);
        PtNode * o  = a1->getInEdge(0, pt)->getFromNode(pt);
        PtNode * a2 = o->getInEdge(0, pt)->getFromNode(pt);

        if ( a2->idx()[PtNode::IDX_VALID] <= 0 )
            continue;

        bool found = false;
        for ( int j = 0; j < sobj.size(); ++j ) {
            if ( a2 == sobj[j] ) {
                found = true;
                break;
            }
        }
        if ( !found ) {
            sobj.push_back(a2);
        }
    }

    return sobj;
}

void ParseTree::getSingleObjDet(std::vector<Detection> &dets, int ptIdx)
{
    const ParseTree &pt(*this);

    std::vector<const PtNode *> sobj = findSingleObjAndNodes();

    for ( int i = 0; i < sobj.size(); ++i) {
        const ParseInfo * info = sobj[i]->parseInfo(&pt);
        Rectangle_<Scalar> bbox(info->x(), info->y(), info->width(), info->height());
        if ( bbox.width() == 0 || bbox.height() == 0 ) {
            const ParseInfo * info1 = sobj[i]->outEdge(0, pt)->toNode(pt)->parseInfo(&pt);
            bbox.setWidth(info1->width());
            bbox.setHeight(info1->height());
            RGM_CHECK_EQ(info->x(), info1->x());
            RGM_CHECK_EQ(info->y(), info1->y());
        }
        int idxG = sobj[i]->idx()[PtNode::IDX_G];
        //const ParseInfo *info2 = sobj[i]->inEdge(0, pt)->fromNode(pt)->parseInfo(&pt);
        Detection det(idxG, info->l_, info->x_, info->y_, info->score_,
                      bbox, ptIdx, sobj[i]->idx()[PtNode::IDX_MYSELF]);
        if ( det.clipBbox(imgWd(), imgHt()) ) {
            dets.push_back(det);
        }
    }
}

void ParseTree::doBboxPred(std::vector<Detection> & dets, int ptIdx)
{
    if ( grammar() == NULL ) {
        RGM_LOG(warning, "No grammar model is specified");
        return;
    }

    if ( grammar()->bboxPred().size() == 0 )
        return;

    const ParseTree &pt(*this);

    std::vector<const PtNode *> sobj = findSingleObjAndNodes();

    for ( int i = 0; i < sobj.size(); ++i ) {
        const PtNode * n = sobj[i];
        // get detection
        const ParseInfo * info = n->parseInfo(&pt);
        Rectangle_<Scalar> bbox(info->x(), info->y(), info->width(), info->height());
        if ( bbox.width() == 0 || bbox.height() == 0 ) {
            const ParseInfo * info1 = n->outEdge(0, pt)->toNode(pt)->parseInfo(&pt);
            bbox.setWidth(info1->width());
            bbox.setHeight(info1->height());
            RGM_CHECK_EQ(info->x(), info1->x());
            RGM_CHECK_EQ(info->y(), info1->y());
        }
        int idxG = n->idx()[PtNode::IDX_G];
        //const ParseInfo *info2 = sobj[i]->inEdge(0, pt)->fromNode(pt)->parseInfo(&pt);
        Detection det(idxG, info->l_, info->x_, info->y_, info->score_,
                      bbox, ptIdx, sobj[i]->idx()[PtNode::IDX_MYSELF]);
        if ( !det.clipBbox(imgWd(), imgHt()) )
            continue;

        // get prediction model
        std::map<int, Matrix>::const_iterator iter = grammar()->bboxPred().find(idxG);
        RGM_CHECK_NOTEQ(iter, grammar()->bboxPred().end());
        const Matrix & pred(iter->second);

        RGM_CHECK_EQ(pred.rows(), n->idxOutEdges().size() * 2 + 1);

        Scalar wd = bbox.width() - 1;
        Scalar ht = bbox.height() - 1; // bug due to the setting in learning the pred. model, to be fixed
        Scalar rx = bbox.x() + wd / 2.0F;
        Scalar ry = bbox.y() + ht / 2.0F;

        Matrix A(Matrix::Zero(1, pred.rows()));
        int c = 0;

        // get detections of all appearance filters
        for ( int j = 0; j < n->idxOutEdges().size(); ++j ) {
            const PtNode * o = n->outEdge(j, pt)->toNode(pt);
            const PtNode * a = o->outEdge(0, pt)->toNode(pt);
            const PtNode * t = a->outEdge(0, pt)->toNode(pt);

            const ParseInfo * tinfo = t->parseInfo(&pt);

            Scalar tx = tinfo->x() + (tinfo->width() - 1.0F)/2.0F;
            Scalar ty = tinfo->y() + (tinfo->height() - 1.0F)/2.0F;

            A(0, c++) = (tx - rx) / wd;
            A(0, c++) = (ty - ry) / ht;
        }

        A(0, c) = 1;

        // compute the predicted bbox
        Matrix dxy = A * pred;

        int x1 = bbox.x() + dxy(0, 0) * wd;
        int y1 = bbox.y() + dxy(0, 1) * ht;
        int x2 = bbox.right() + dxy(0, 2) * wd;
        int y2 = bbox.bottom() + dxy(0, 3) * ht;

        det.setX(x1);
        det.setY(y1);
        det.setWidth(x2 - x1 + 1);
        det.setHeight(y2 - y1 + 1);

        if ( det.clipBbox(imgWd(), imgHt()) ) {
            dets.push_back(det);
        }
    }
}

template<class Archive>
void ParseTree::serialize(Archive & ar, const unsigned int version)
{
    ar.register_type(static_cast<PtNode *>(NULL));
    ar.register_type(static_cast<PtEdge *>(NULL));
    ar.register_type(static_cast<Appearance::Param *>(NULL));
    ar.register_type(static_cast<Deformation::Param *>(NULL));
    ar.register_type(static_cast<Scaleprior::Param *>(NULL));
    ar.register_type(static_cast<ParseInfo *>(NULL));
    ar.register_type(static_cast<States *>(NULL));

    ar.template register_type<PtNode>();
    ar.template register_type<PtEdge>();
    ar.template register_type<Appearance::Param>();
    ar.template register_type<Deformation::Param>();
    ar.template register_type<Scaleprior::Param>();
    ar.template register_type<ParseInfo>();
    ar.template register_type<States>();

    ar & BOOST_SERIALIZATION_NVP(nodeSet_);
    ar & BOOST_SERIALIZATION_NVP(edgeSet_);
    ar & BOOST_SERIALIZATION_NVP(idxRootNode_);
    ar & BOOST_SERIALIZATION_NVP(appearanceSet_);
    ar & BOOST_SERIALIZATION_NVP(biasSet_);
    ar & BOOST_SERIALIZATION_NVP(deformationSet_);
    ar & BOOST_SERIALIZATION_NVP(scalepriorSet_);
    ar & BOOST_SERIALIZATION_NVP(parseInfoSet_);
    ar & BOOST_SERIALIZATION_NVP(dataId_);
    ar & BOOST_SERIALIZATION_NVP(states_);
    ar & BOOST_SERIALIZATION_NVP(appearanceX_);
    ar & BOOST_SERIALIZATION_NVP(imgWd_);
    ar & BOOST_SERIALIZATION_NVP(imgHt_);
}

INSTANTIATE_BOOST_SERIALIZATION(ParseTree);



// ------ TrainExample ------

TrainExample::TrainExample() :
    marginBound_(0), beliefNorm_(0), maxNonbeliefNorm_(0), nbHist_(0)
{
}

TrainExample::TrainExample(const TrainExample & ex) :
    marginBound_(ex.marginBound()), beliefNorm_(ex.beliefNorm()), maxNonbeliefNorm_(ex.maxNonbeliefNorm()),
    nbHist_(ex.nbHist())
{
    pts_.resize(ex.pts().size());
    for ( int i = 0; i < pts_.size(); ++i ) {
        pts_[i] = ex.pts()[i];
    }
}

TrainExample & TrainExample::operator=(const TrainExample & ex)
{
    if ( this == &ex ) {
        return *this;
    }

    pts_.resize(ex.pts().size());
    for ( int i = 0; i < pts_.size(); ++i ) {
        pts_[i] = ex.pts()[i];
    }

    marginBound_ = ex.marginBound();
    beliefNorm_ = ex.beliefNorm();
    maxNonbeliefNorm_ = ex.maxNonbeliefNorm();
    nbHist_ = ex.nbHist();

    return *this;
}

const std::vector<ParseTree> & TrainExample::pts() const
{
    return pts_;
}

std::vector<ParseTree> & TrainExample::getPts()
{
    return pts_;
}

Scalar   TrainExample::marginBound() const
{
    return marginBound_;
}

Scalar & TrainExample::getMarginBound()
{
    return marginBound_;
}

Scalar   TrainExample::beliefNorm() const
{
    return beliefNorm_;
}

Scalar & TrainExample::getBeliefNorm()
{
    return beliefNorm_;
}

Scalar   TrainExample::maxNonbeliefNorm() const
{
    return maxNonbeliefNorm_;
}

Scalar & TrainExample::getMaxNonbeliefNorm()
{
    return maxNonbeliefNorm_;
}

int   TrainExample::nbHist() const
{
    return nbHist_;
}

int & TrainExample::getNbHist()
{
    return nbHist_;
}

bool TrainExample::isEqual(const ParseTree & pt) const
{    
    assert(pts().size() >= 2);

    const ParseInfo * p = pts()[0].rootNode()->parseInfo(&pts()[0]);
    const ParseInfo * p1 = pt.rootNode()->parseInfo(&pt);

    assert(p);
    assert(p1);

    return (pts()[0].dataId() == pt.dataId()) &&
           (p->c_ == p1->c_) && (p->l_ == p1->l_) &&
            (p->x_ == p1->x_) && (p->y_ == p1->y_ );
}

template<class Archive>
void TrainExample::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(pts_);
    ar & BOOST_SERIALIZATION_NVP(marginBound_);
    ar & BOOST_SERIALIZATION_NVP(beliefNorm_);
    ar & BOOST_SERIALIZATION_NVP(maxNonbeliefNorm_);
    ar & BOOST_SERIALIZATION_NVP(nbHist_);
}

INSTANTIATE_BOOST_SERIALIZATION(TrainExample);



// ------- PtIntersector ------

PtIntersector::PtIntersector(const ParseTree &  reference, Scalar threshold, bool dividedByUnion) :
    reference_(&reference), threshold_(threshold), dividedByUnion_(dividedByUnion)
{
    assert(reference_ != NULL);
}

bool PtIntersector::operator()(const ParseTree &  pt, Scalar * score) const
{
    if (score) {
        *score = 0.0;
    }

    const ParseInfo * ref = reference_->rootNode()->parseInfo( reference_ );
    const ParseInfo * cur = pt.rootNode()->parseInfo( &pt );

    const int left = std::max<int>(ref->left(), cur->left());
    const int right = std::min<int>(ref->right(), cur->right());

    if (right < left) {
        return false;
    }

    const int top = std::max<int>(ref->top(), cur->top());
    const int bottom = std::min<int>(ref->bottom(), cur->bottom());

    if (bottom < top) {
        return false;
    }

    const int intersectionArea = (right - left + 1) * (bottom - top + 1);
    const int rectArea = cur->area();

    if (dividedByUnion_) {
        const int referenceArea = ref->area();
        const int unionArea = referenceArea + rectArea - intersectionArea;

        if (intersectionArea >= unionArea * threshold_) {
            if (score) {
                *score = static_cast<Scalar>(intersectionArea) / unionArea;
            }

            return true;
        }
    } else {
        if (intersectionArea >= rectArea * threshold_) {
            if (score) {
                *score = static_cast<Scalar>(intersectionArea) / rectArea;
            }

            return true;
        }
    }

    return false;
}


} //namespace RGM

