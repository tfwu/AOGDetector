/// This file is adapted from UCLA ImageParser by Brandon Rothrock (rothrock@cs.ucla.edu)

#ifndef RGM_XMLREADER_HPP_
#define RGM_XMLREADER_HPP_

#include <vector>
#include <string>
#include <string.h>

#include "RapidXML/rapidxml.hpp"

namespace RGM
{

/// Read xml file
class XMLData
{
public:
    struct XMLNode {
        XMLNode();
        XMLNode(char* _name, char* _value);

        ~XMLNode();

        XMLNode* Clone();

        void SetName(char* _name);
        void SetValue(char* _value);
        void SetChildren(std::vector<XMLNode*>& _children);
        void AddChild(XMLNode* pNode);

        char* _name;
        char* _value;
        int _childCount;
        XMLNode* _parent;
        XMLNode** _children;
    };

    XMLData();

    void Clear();

    void ReadFromFile(std::string fileName, std::string firstNodeName);
    void ReadFromString(char* str, std::string firstNodeName);
    void WriteToFile(std::string fileName);

    XMLNode* GetNode(std::string path);
    XMLNode* GetNode(std::string path, XMLNode* pNode);
    XMLNode* GetRootNode();

    XMLNode* FindFirst(std::string path, std::string childName, std::string childValue);
    XMLNode* FindFirst(std::string path, std::string childName, std::string field, std::string query);
    XMLNode* FindFirst(std::string basePath, std::string childName, std::string field1, std::string query1, std::string field2, std::string query2);

    std::vector<XMLNode *> GetNodes(const std::string& path);
    std::vector<XMLNode *> GetNodes(const std::string& path, XMLNode* pNode);

    std::string GetNodeName(XMLNode* pNode);
    std::string GetString(std::string path);
    int GetInt(std::string path);
    float GetFloat(std::string path);
    bool GetBoolean(std::string path);

    std::string GetString(XMLNode* pNode);
    int   GetInt(XMLNode* pNode);
    float GetFloat(XMLNode* pNode);
    bool  GetBoolean(XMLNode* pNode);

    std::string GetString(std::string path, XMLNode* pNode);
    int   GetInt(std::string path, XMLNode* pNode);
    float GetFloat(std::string path, XMLNode* pNode);
    bool  GetBoolean(std::string path, XMLNode* pNode);

    void RemoveNode(XMLNode* pNode);
    void InsertNode(XMLNode* pParent, XMLNode* pNode);

private:
    XMLNode* ReadNode(rapidxml::xml_node<>* node);
    void     WriteNode(XMLNode* node, std::ofstream& ofs);

    XMLNode* _pRootNode;
};

} // namespace RGM

#endif // RGM_XMLREADER_HPP_

