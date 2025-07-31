//
// Generated file, do not edit! Created by opp_msgtool 6.1 from clustering.msg.
//

// Disable warnings about unused variables, empty switch stmts, etc:
#ifdef _MSC_VER
#  pragma warning(disable:4101)
#  pragma warning(disable:4065)
#endif

#if defined(__clang__)
#  pragma clang diagnostic ignored "-Wshadow"
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wunused-parameter"
#  pragma clang diagnostic ignored "-Wc++98-compat"
#  pragma clang diagnostic ignored "-Wunreachable-code-break"
#  pragma clang diagnostic ignored "-Wold-style-cast"
#elif defined(__GNUC__)
#  pragma GCC diagnostic ignored "-Wshadow"
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wunused-parameter"
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#  pragma GCC diagnostic ignored "-Wfloat-conversion"
#endif

#include <iostream>
#include <sstream>
#include <memory>
#include <type_traits>
#include "clustering_m.h"

namespace omnetpp {

// Template pack/unpack rules. They are declared *after* a1l type-specific pack functions for multiple reasons.
// They are in the omnetpp namespace, to allow them to be found by argument-dependent lookup via the cCommBuffer argument

// Packing/unpacking an std::vector
template<typename T, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::vector<T,A>& v)
{
    int n = v.size();
    doParsimPacking(buffer, n);
    for (int i = 0; i < n; i++)
        doParsimPacking(buffer, v[i]);
}

template<typename T, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::vector<T,A>& v)
{
    int n;
    doParsimUnpacking(buffer, n);
    v.resize(n);
    for (int i = 0; i < n; i++)
        doParsimUnpacking(buffer, v[i]);
}

// Packing/unpacking an std::list
template<typename T, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::list<T,A>& l)
{
    doParsimPacking(buffer, (int)l.size());
    for (typename std::list<T,A>::const_iterator it = l.begin(); it != l.end(); ++it)
        doParsimPacking(buffer, (T&)*it);
}

template<typename T, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::list<T,A>& l)
{
    int n;
    doParsimUnpacking(buffer, n);
    for (int i = 0; i < n; i++) {
        l.push_back(T());
        doParsimUnpacking(buffer, l.back());
    }
}

// Packing/unpacking an std::set
template<typename T, typename Tr, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::set<T,Tr,A>& s)
{
    doParsimPacking(buffer, (int)s.size());
    for (typename std::set<T,Tr,A>::const_iterator it = s.begin(); it != s.end(); ++it)
        doParsimPacking(buffer, *it);
}

template<typename T, typename Tr, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::set<T,Tr,A>& s)
{
    int n;
    doParsimUnpacking(buffer, n);
    for (int i = 0; i < n; i++) {
        T x;
        doParsimUnpacking(buffer, x);
        s.insert(x);
    }
}

// Packing/unpacking an std::map
template<typename K, typename V, typename Tr, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::map<K,V,Tr,A>& m)
{
    doParsimPacking(buffer, (int)m.size());
    for (typename std::map<K,V,Tr,A>::const_iterator it = m.begin(); it != m.end(); ++it) {
        doParsimPacking(buffer, it->first);
        doParsimPacking(buffer, it->second);
    }
}

template<typename K, typename V, typename Tr, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::map<K,V,Tr,A>& m)
{
    int n;
    doParsimUnpacking(buffer, n);
    for (int i = 0; i < n; i++) {
        K k; V v;
        doParsimUnpacking(buffer, k);
        doParsimUnpacking(buffer, v);
        m[k] = v;
    }
}

// Default pack/unpack function for arrays
template<typename T>
void doParsimArrayPacking(omnetpp::cCommBuffer *b, const T *t, int n)
{
    for (int i = 0; i < n; i++)
        doParsimPacking(b, t[i]);
}

template<typename T>
void doParsimArrayUnpacking(omnetpp::cCommBuffer *b, T *t, int n)
{
    for (int i = 0; i < n; i++)
        doParsimUnpacking(b, t[i]);
}

// Default rule to prevent compiler from choosing base class' doParsimPacking() function
template<typename T>
void doParsimPacking(omnetpp::cCommBuffer *, const T& t)
{
    throw omnetpp::cRuntimeError("Parsim error: No doParsimPacking() function for type %s", omnetpp::opp_typename(typeid(t)));
}

template<typename T>
void doParsimUnpacking(omnetpp::cCommBuffer *, T& t)
{
    throw omnetpp::cRuntimeError("Parsim error: No doParsimUnpacking() function for type %s", omnetpp::opp_typename(typeid(t)));
}

}  // namespace omnetpp

Register_Class(ClusterInfoMsg)

ClusterInfoMsg::ClusterInfoMsg(const char *name, short kind) : ::omnetpp::cMessage(name, kind)
{
}

ClusterInfoMsg::ClusterInfoMsg(const ClusterInfoMsg& other) : ::omnetpp::cMessage(other)
{
    copy(other);
}

ClusterInfoMsg::~ClusterInfoMsg()
{
}

ClusterInfoMsg& ClusterInfoMsg::operator=(const ClusterInfoMsg& other)
{
    if (this == &other) return *this;
    ::omnetpp::cMessage::operator=(other);
    copy(other);
    return *this;
}

void ClusterInfoMsg::copy(const ClusterInfoMsg& other)
{
    this->sender = other.sender;
    this->clusterID = other.clusterID;
    this->CIP = other.CIP;
    this->numNodes = other.numNodes;
    this->numEdges = other.numEdges;
    this->mergedNumNodes = other.mergedNumNodes;
    this->mergedNumEdges = other.mergedNumEdges;
}

void ClusterInfoMsg::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::omnetpp::cMessage::parsimPack(b);
    doParsimPacking(b,this->sender);
    doParsimPacking(b,this->clusterID);
    doParsimPacking(b,this->CIP);
    doParsimPacking(b,this->numNodes);
    doParsimPacking(b,this->numEdges);
    doParsimPacking(b,this->mergedNumNodes);
    doParsimPacking(b,this->mergedNumEdges);
}

void ClusterInfoMsg::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::omnetpp::cMessage::parsimUnpack(b);
    doParsimUnpacking(b,this->sender);
    doParsimUnpacking(b,this->clusterID);
    doParsimUnpacking(b,this->CIP);
    doParsimUnpacking(b,this->numNodes);
    doParsimUnpacking(b,this->numEdges);
    doParsimUnpacking(b,this->mergedNumNodes);
    doParsimUnpacking(b,this->mergedNumEdges);
}

const char * ClusterInfoMsg::getSender() const
{
    return this->sender.c_str();
}

void ClusterInfoMsg::setSender(const char * sender)
{
    this->sender = sender;
}

const char * ClusterInfoMsg::getClusterID() const
{
    return this->clusterID.c_str();
}

void ClusterInfoMsg::setClusterID(const char * clusterID)
{
    this->clusterID = clusterID;
}

const char * ClusterInfoMsg::getCIP() const
{
    return this->CIP.c_str();
}

void ClusterInfoMsg::setCIP(const char * CIP)
{
    this->CIP = CIP;
}

int ClusterInfoMsg::getNumNodes() const
{
    return this->numNodes;
}

void ClusterInfoMsg::setNumNodes(int numNodes)
{
    this->numNodes = numNodes;
}

int ClusterInfoMsg::getNumEdges() const
{
    return this->numEdges;
}

void ClusterInfoMsg::setNumEdges(int numEdges)
{
    this->numEdges = numEdges;
}

int ClusterInfoMsg::getMergedNumNodes() const
{
    return this->mergedNumNodes;
}

void ClusterInfoMsg::setMergedNumNodes(int mergedNumNodes)
{
    this->mergedNumNodes = mergedNumNodes;
}

int ClusterInfoMsg::getMergedNumEdges() const
{
    return this->mergedNumEdges;
}

void ClusterInfoMsg::setMergedNumEdges(int mergedNumEdges)
{
    this->mergedNumEdges = mergedNumEdges;
}

class ClusterInfoMsgDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_sender,
        FIELD_clusterID,
        FIELD_CIP,
        FIELD_numNodes,
        FIELD_numEdges,
        FIELD_mergedNumNodes,
        FIELD_mergedNumEdges,
    };
  public:
    ClusterInfoMsgDescriptor();
    virtual ~ClusterInfoMsgDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(ClusterInfoMsgDescriptor)

ClusterInfoMsgDescriptor::ClusterInfoMsgDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(ClusterInfoMsg)), "omnetpp::cMessage")
{
    propertyNames = nullptr;
}

ClusterInfoMsgDescriptor::~ClusterInfoMsgDescriptor()
{
    delete[] propertyNames;
}

bool ClusterInfoMsgDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<ClusterInfoMsg *>(obj)!=nullptr;
}

const char **ClusterInfoMsgDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *ClusterInfoMsgDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int ClusterInfoMsgDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 7+base->getFieldCount() : 7;
}

unsigned int ClusterInfoMsgDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,    // FIELD_sender
        FD_ISEDITABLE,    // FIELD_clusterID
        FD_ISEDITABLE,    // FIELD_CIP
        FD_ISEDITABLE,    // FIELD_numNodes
        FD_ISEDITABLE,    // FIELD_numEdges
        FD_ISEDITABLE,    // FIELD_mergedNumNodes
        FD_ISEDITABLE,    // FIELD_mergedNumEdges
    };
    return (field >= 0 && field < 7) ? fieldTypeFlags[field] : 0;
}

const char *ClusterInfoMsgDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "sender",
        "clusterID",
        "CIP",
        "numNodes",
        "numEdges",
        "mergedNumNodes",
        "mergedNumEdges",
    };
    return (field >= 0 && field < 7) ? fieldNames[field] : nullptr;
}

int ClusterInfoMsgDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "sender") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "clusterID") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "CIP") == 0) return baseIndex + 2;
    if (strcmp(fieldName, "numNodes") == 0) return baseIndex + 3;
    if (strcmp(fieldName, "numEdges") == 0) return baseIndex + 4;
    if (strcmp(fieldName, "mergedNumNodes") == 0) return baseIndex + 5;
    if (strcmp(fieldName, "mergedNumEdges") == 0) return baseIndex + 6;
    return base ? base->findField(fieldName) : -1;
}

const char *ClusterInfoMsgDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "string",    // FIELD_sender
        "string",    // FIELD_clusterID
        "string",    // FIELD_CIP
        "int",    // FIELD_numNodes
        "int",    // FIELD_numEdges
        "int",    // FIELD_mergedNumNodes
        "int",    // FIELD_mergedNumEdges
    };
    return (field >= 0 && field < 7) ? fieldTypeStrings[field] : nullptr;
}

const char **ClusterInfoMsgDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *ClusterInfoMsgDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int ClusterInfoMsgDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    ClusterInfoMsg *pp = omnetpp::fromAnyPtr<ClusterInfoMsg>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void ClusterInfoMsgDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    ClusterInfoMsg *pp = omnetpp::fromAnyPtr<ClusterInfoMsg>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'ClusterInfoMsg'", field);
    }
}

const char *ClusterInfoMsgDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    ClusterInfoMsg *pp = omnetpp::fromAnyPtr<ClusterInfoMsg>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string ClusterInfoMsgDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    ClusterInfoMsg *pp = omnetpp::fromAnyPtr<ClusterInfoMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: return oppstring2string(pp->getSender());
        case FIELD_clusterID: return oppstring2string(pp->getClusterID());
        case FIELD_CIP: return oppstring2string(pp->getCIP());
        case FIELD_numNodes: return long2string(pp->getNumNodes());
        case FIELD_numEdges: return long2string(pp->getNumEdges());
        case FIELD_mergedNumNodes: return long2string(pp->getMergedNumNodes());
        case FIELD_mergedNumEdges: return long2string(pp->getMergedNumEdges());
        default: return "";
    }
}

void ClusterInfoMsgDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    ClusterInfoMsg *pp = omnetpp::fromAnyPtr<ClusterInfoMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: pp->setSender((value)); break;
        case FIELD_clusterID: pp->setClusterID((value)); break;
        case FIELD_CIP: pp->setCIP((value)); break;
        case FIELD_numNodes: pp->setNumNodes(string2long(value)); break;
        case FIELD_numEdges: pp->setNumEdges(string2long(value)); break;
        case FIELD_mergedNumNodes: pp->setMergedNumNodes(string2long(value)); break;
        case FIELD_mergedNumEdges: pp->setMergedNumEdges(string2long(value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'ClusterInfoMsg'", field);
    }
}

omnetpp::cValue ClusterInfoMsgDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    ClusterInfoMsg *pp = omnetpp::fromAnyPtr<ClusterInfoMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: return pp->getSender();
        case FIELD_clusterID: return pp->getClusterID();
        case FIELD_CIP: return pp->getCIP();
        case FIELD_numNodes: return pp->getNumNodes();
        case FIELD_numEdges: return pp->getNumEdges();
        case FIELD_mergedNumNodes: return pp->getMergedNumNodes();
        case FIELD_mergedNumEdges: return pp->getMergedNumEdges();
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'ClusterInfoMsg' as cValue -- field index out of range?", field);
    }
}

void ClusterInfoMsgDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    ClusterInfoMsg *pp = omnetpp::fromAnyPtr<ClusterInfoMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: pp->setSender(value.stringValue()); break;
        case FIELD_clusterID: pp->setClusterID(value.stringValue()); break;
        case FIELD_CIP: pp->setCIP(value.stringValue()); break;
        case FIELD_numNodes: pp->setNumNodes(omnetpp::checked_int_cast<int>(value.intValue())); break;
        case FIELD_numEdges: pp->setNumEdges(omnetpp::checked_int_cast<int>(value.intValue())); break;
        case FIELD_mergedNumNodes: pp->setMergedNumNodes(omnetpp::checked_int_cast<int>(value.intValue())); break;
        case FIELD_mergedNumEdges: pp->setMergedNumEdges(omnetpp::checked_int_cast<int>(value.intValue())); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'ClusterInfoMsg'", field);
    }
}

const char *ClusterInfoMsgDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    };
}

omnetpp::any_ptr ClusterInfoMsgDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    ClusterInfoMsg *pp = omnetpp::fromAnyPtr<ClusterInfoMsg>(object); (void)pp;
    switch (field) {
        default: return omnetpp::any_ptr(nullptr);
    }
}

void ClusterInfoMsgDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    ClusterInfoMsg *pp = omnetpp::fromAnyPtr<ClusterInfoMsg>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'ClusterInfoMsg'", field);
    }
}

Register_Class(ClusterUpdateMsg)

ClusterUpdateMsg::ClusterUpdateMsg(const char *name, short kind) : ::omnetpp::cMessage(name, kind)
{
}

ClusterUpdateMsg::ClusterUpdateMsg(const ClusterUpdateMsg& other) : ::omnetpp::cMessage(other)
{
    copy(other);
}

ClusterUpdateMsg::~ClusterUpdateMsg()
{
}

ClusterUpdateMsg& ClusterUpdateMsg::operator=(const ClusterUpdateMsg& other)
{
    if (this == &other) return *this;
    ::omnetpp::cMessage::operator=(other);
    copy(other);
    return *this;
}

void ClusterUpdateMsg::copy(const ClusterUpdateMsg& other)
{
    this->newClusterID = other.newClusterID;
    this->newCIP = other.newCIP;
}

void ClusterUpdateMsg::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::omnetpp::cMessage::parsimPack(b);
    doParsimPacking(b,this->newClusterID);
    doParsimPacking(b,this->newCIP);
}

void ClusterUpdateMsg::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::omnetpp::cMessage::parsimUnpack(b);
    doParsimUnpacking(b,this->newClusterID);
    doParsimUnpacking(b,this->newCIP);
}

const char * ClusterUpdateMsg::getNewClusterID() const
{
    return this->newClusterID.c_str();
}

void ClusterUpdateMsg::setNewClusterID(const char * newClusterID)
{
    this->newClusterID = newClusterID;
}

const char * ClusterUpdateMsg::getNewCIP() const
{
    return this->newCIP.c_str();
}

void ClusterUpdateMsg::setNewCIP(const char * newCIP)
{
    this->newCIP = newCIP;
}

class ClusterUpdateMsgDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_newClusterID,
        FIELD_newCIP,
    };
  public:
    ClusterUpdateMsgDescriptor();
    virtual ~ClusterUpdateMsgDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(ClusterUpdateMsgDescriptor)

ClusterUpdateMsgDescriptor::ClusterUpdateMsgDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(ClusterUpdateMsg)), "omnetpp::cMessage")
{
    propertyNames = nullptr;
}

ClusterUpdateMsgDescriptor::~ClusterUpdateMsgDescriptor()
{
    delete[] propertyNames;
}

bool ClusterUpdateMsgDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<ClusterUpdateMsg *>(obj)!=nullptr;
}

const char **ClusterUpdateMsgDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *ClusterUpdateMsgDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int ClusterUpdateMsgDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 2+base->getFieldCount() : 2;
}

unsigned int ClusterUpdateMsgDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,    // FIELD_newClusterID
        FD_ISEDITABLE,    // FIELD_newCIP
    };
    return (field >= 0 && field < 2) ? fieldTypeFlags[field] : 0;
}

const char *ClusterUpdateMsgDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "newClusterID",
        "newCIP",
    };
    return (field >= 0 && field < 2) ? fieldNames[field] : nullptr;
}

int ClusterUpdateMsgDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "newClusterID") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "newCIP") == 0) return baseIndex + 1;
    return base ? base->findField(fieldName) : -1;
}

const char *ClusterUpdateMsgDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "string",    // FIELD_newClusterID
        "string",    // FIELD_newCIP
    };
    return (field >= 0 && field < 2) ? fieldTypeStrings[field] : nullptr;
}

const char **ClusterUpdateMsgDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *ClusterUpdateMsgDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int ClusterUpdateMsgDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    ClusterUpdateMsg *pp = omnetpp::fromAnyPtr<ClusterUpdateMsg>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void ClusterUpdateMsgDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    ClusterUpdateMsg *pp = omnetpp::fromAnyPtr<ClusterUpdateMsg>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'ClusterUpdateMsg'", field);
    }
}

const char *ClusterUpdateMsgDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    ClusterUpdateMsg *pp = omnetpp::fromAnyPtr<ClusterUpdateMsg>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string ClusterUpdateMsgDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    ClusterUpdateMsg *pp = omnetpp::fromAnyPtr<ClusterUpdateMsg>(object); (void)pp;
    switch (field) {
        case FIELD_newClusterID: return oppstring2string(pp->getNewClusterID());
        case FIELD_newCIP: return oppstring2string(pp->getNewCIP());
        default: return "";
    }
}

void ClusterUpdateMsgDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    ClusterUpdateMsg *pp = omnetpp::fromAnyPtr<ClusterUpdateMsg>(object); (void)pp;
    switch (field) {
        case FIELD_newClusterID: pp->setNewClusterID((value)); break;
        case FIELD_newCIP: pp->setNewCIP((value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'ClusterUpdateMsg'", field);
    }
}

omnetpp::cValue ClusterUpdateMsgDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    ClusterUpdateMsg *pp = omnetpp::fromAnyPtr<ClusterUpdateMsg>(object); (void)pp;
    switch (field) {
        case FIELD_newClusterID: return pp->getNewClusterID();
        case FIELD_newCIP: return pp->getNewCIP();
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'ClusterUpdateMsg' as cValue -- field index out of range?", field);
    }
}

void ClusterUpdateMsgDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    ClusterUpdateMsg *pp = omnetpp::fromAnyPtr<ClusterUpdateMsg>(object); (void)pp;
    switch (field) {
        case FIELD_newClusterID: pp->setNewClusterID(value.stringValue()); break;
        case FIELD_newCIP: pp->setNewCIP(value.stringValue()); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'ClusterUpdateMsg'", field);
    }
}

const char *ClusterUpdateMsgDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    };
}

omnetpp::any_ptr ClusterUpdateMsgDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    ClusterUpdateMsg *pp = omnetpp::fromAnyPtr<ClusterUpdateMsg>(object); (void)pp;
    switch (field) {
        default: return omnetpp::any_ptr(nullptr);
    }
}

void ClusterUpdateMsgDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    ClusterUpdateMsg *pp = omnetpp::fromAnyPtr<ClusterUpdateMsg>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'ClusterUpdateMsg'", field);
    }
}

Register_Class(MergeRequestMsg)

MergeRequestMsg::MergeRequestMsg(const char *name, short kind) : ::omnetpp::cMessage(name, kind)
{
}

MergeRequestMsg::MergeRequestMsg(const MergeRequestMsg& other) : ::omnetpp::cMessage(other)
{
    copy(other);
}

MergeRequestMsg::~MergeRequestMsg()
{
}

MergeRequestMsg& MergeRequestMsg::operator=(const MergeRequestMsg& other)
{
    if (this == &other) return *this;
    ::omnetpp::cMessage::operator=(other);
    copy(other);
    return *this;
}

void MergeRequestMsg::copy(const MergeRequestMsg& other)
{
    this->sender = other.sender;
    this->clusterID = other.clusterID;
}

void MergeRequestMsg::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::omnetpp::cMessage::parsimPack(b);
    doParsimPacking(b,this->sender);
    doParsimPacking(b,this->clusterID);
}

void MergeRequestMsg::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::omnetpp::cMessage::parsimUnpack(b);
    doParsimUnpacking(b,this->sender);
    doParsimUnpacking(b,this->clusterID);
}

const char * MergeRequestMsg::getSender() const
{
    return this->sender.c_str();
}

void MergeRequestMsg::setSender(const char * sender)
{
    this->sender = sender;
}

const char * MergeRequestMsg::getClusterID() const
{
    return this->clusterID.c_str();
}

void MergeRequestMsg::setClusterID(const char * clusterID)
{
    this->clusterID = clusterID;
}

class MergeRequestMsgDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_sender,
        FIELD_clusterID,
    };
  public:
    MergeRequestMsgDescriptor();
    virtual ~MergeRequestMsgDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(MergeRequestMsgDescriptor)

MergeRequestMsgDescriptor::MergeRequestMsgDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(MergeRequestMsg)), "omnetpp::cMessage")
{
    propertyNames = nullptr;
}

MergeRequestMsgDescriptor::~MergeRequestMsgDescriptor()
{
    delete[] propertyNames;
}

bool MergeRequestMsgDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<MergeRequestMsg *>(obj)!=nullptr;
}

const char **MergeRequestMsgDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *MergeRequestMsgDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int MergeRequestMsgDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 2+base->getFieldCount() : 2;
}

unsigned int MergeRequestMsgDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,    // FIELD_sender
        FD_ISEDITABLE,    // FIELD_clusterID
    };
    return (field >= 0 && field < 2) ? fieldTypeFlags[field] : 0;
}

const char *MergeRequestMsgDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "sender",
        "clusterID",
    };
    return (field >= 0 && field < 2) ? fieldNames[field] : nullptr;
}

int MergeRequestMsgDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "sender") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "clusterID") == 0) return baseIndex + 1;
    return base ? base->findField(fieldName) : -1;
}

const char *MergeRequestMsgDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "string",    // FIELD_sender
        "string",    // FIELD_clusterID
    };
    return (field >= 0 && field < 2) ? fieldTypeStrings[field] : nullptr;
}

const char **MergeRequestMsgDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *MergeRequestMsgDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int MergeRequestMsgDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    MergeRequestMsg *pp = omnetpp::fromAnyPtr<MergeRequestMsg>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void MergeRequestMsgDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeRequestMsg *pp = omnetpp::fromAnyPtr<MergeRequestMsg>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'MergeRequestMsg'", field);
    }
}

const char *MergeRequestMsgDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    MergeRequestMsg *pp = omnetpp::fromAnyPtr<MergeRequestMsg>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string MergeRequestMsgDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    MergeRequestMsg *pp = omnetpp::fromAnyPtr<MergeRequestMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: return oppstring2string(pp->getSender());
        case FIELD_clusterID: return oppstring2string(pp->getClusterID());
        default: return "";
    }
}

void MergeRequestMsgDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeRequestMsg *pp = omnetpp::fromAnyPtr<MergeRequestMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: pp->setSender((value)); break;
        case FIELD_clusterID: pp->setClusterID((value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'MergeRequestMsg'", field);
    }
}

omnetpp::cValue MergeRequestMsgDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    MergeRequestMsg *pp = omnetpp::fromAnyPtr<MergeRequestMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: return pp->getSender();
        case FIELD_clusterID: return pp->getClusterID();
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'MergeRequestMsg' as cValue -- field index out of range?", field);
    }
}

void MergeRequestMsgDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeRequestMsg *pp = omnetpp::fromAnyPtr<MergeRequestMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: pp->setSender(value.stringValue()); break;
        case FIELD_clusterID: pp->setClusterID(value.stringValue()); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'MergeRequestMsg'", field);
    }
}

const char *MergeRequestMsgDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    };
}

omnetpp::any_ptr MergeRequestMsgDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    MergeRequestMsg *pp = omnetpp::fromAnyPtr<MergeRequestMsg>(object); (void)pp;
    switch (field) {
        default: return omnetpp::any_ptr(nullptr);
    }
}

void MergeRequestMsgDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeRequestMsg *pp = omnetpp::fromAnyPtr<MergeRequestMsg>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'MergeRequestMsg'", field);
    }
}

Register_Class(MergeAckMsg)

MergeAckMsg::MergeAckMsg(const char *name, short kind) : ::omnetpp::cMessage(name, kind)
{
}

MergeAckMsg::MergeAckMsg(const MergeAckMsg& other) : ::omnetpp::cMessage(other)
{
    copy(other);
}

MergeAckMsg::~MergeAckMsg()
{
}

MergeAckMsg& MergeAckMsg::operator=(const MergeAckMsg& other)
{
    if (this == &other) return *this;
    ::omnetpp::cMessage::operator=(other);
    copy(other);
    return *this;
}

void MergeAckMsg::copy(const MergeAckMsg& other)
{
    this->sender = other.sender;
    this->acceptedCluster = other.acceptedCluster;
}

void MergeAckMsg::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::omnetpp::cMessage::parsimPack(b);
    doParsimPacking(b,this->sender);
    doParsimPacking(b,this->acceptedCluster);
}

void MergeAckMsg::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::omnetpp::cMessage::parsimUnpack(b);
    doParsimUnpacking(b,this->sender);
    doParsimUnpacking(b,this->acceptedCluster);
}

const char * MergeAckMsg::getSender() const
{
    return this->sender.c_str();
}

void MergeAckMsg::setSender(const char * sender)
{
    this->sender = sender;
}

const char * MergeAckMsg::getAcceptedCluster() const
{
    return this->acceptedCluster.c_str();
}

void MergeAckMsg::setAcceptedCluster(const char * acceptedCluster)
{
    this->acceptedCluster = acceptedCluster;
}

class MergeAckMsgDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_sender,
        FIELD_acceptedCluster,
    };
  public:
    MergeAckMsgDescriptor();
    virtual ~MergeAckMsgDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(MergeAckMsgDescriptor)

MergeAckMsgDescriptor::MergeAckMsgDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(MergeAckMsg)), "omnetpp::cMessage")
{
    propertyNames = nullptr;
}

MergeAckMsgDescriptor::~MergeAckMsgDescriptor()
{
    delete[] propertyNames;
}

bool MergeAckMsgDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<MergeAckMsg *>(obj)!=nullptr;
}

const char **MergeAckMsgDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *MergeAckMsgDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int MergeAckMsgDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 2+base->getFieldCount() : 2;
}

unsigned int MergeAckMsgDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,    // FIELD_sender
        FD_ISEDITABLE,    // FIELD_acceptedCluster
    };
    return (field >= 0 && field < 2) ? fieldTypeFlags[field] : 0;
}

const char *MergeAckMsgDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "sender",
        "acceptedCluster",
    };
    return (field >= 0 && field < 2) ? fieldNames[field] : nullptr;
}

int MergeAckMsgDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "sender") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "acceptedCluster") == 0) return baseIndex + 1;
    return base ? base->findField(fieldName) : -1;
}

const char *MergeAckMsgDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "string",    // FIELD_sender
        "string",    // FIELD_acceptedCluster
    };
    return (field >= 0 && field < 2) ? fieldTypeStrings[field] : nullptr;
}

const char **MergeAckMsgDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *MergeAckMsgDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int MergeAckMsgDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    MergeAckMsg *pp = omnetpp::fromAnyPtr<MergeAckMsg>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void MergeAckMsgDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeAckMsg *pp = omnetpp::fromAnyPtr<MergeAckMsg>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'MergeAckMsg'", field);
    }
}

const char *MergeAckMsgDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    MergeAckMsg *pp = omnetpp::fromAnyPtr<MergeAckMsg>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string MergeAckMsgDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    MergeAckMsg *pp = omnetpp::fromAnyPtr<MergeAckMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: return oppstring2string(pp->getSender());
        case FIELD_acceptedCluster: return oppstring2string(pp->getAcceptedCluster());
        default: return "";
    }
}

void MergeAckMsgDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeAckMsg *pp = omnetpp::fromAnyPtr<MergeAckMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: pp->setSender((value)); break;
        case FIELD_acceptedCluster: pp->setAcceptedCluster((value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'MergeAckMsg'", field);
    }
}

omnetpp::cValue MergeAckMsgDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    MergeAckMsg *pp = omnetpp::fromAnyPtr<MergeAckMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: return pp->getSender();
        case FIELD_acceptedCluster: return pp->getAcceptedCluster();
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'MergeAckMsg' as cValue -- field index out of range?", field);
    }
}

void MergeAckMsgDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeAckMsg *pp = omnetpp::fromAnyPtr<MergeAckMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: pp->setSender(value.stringValue()); break;
        case FIELD_acceptedCluster: pp->setAcceptedCluster(value.stringValue()); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'MergeAckMsg'", field);
    }
}

const char *MergeAckMsgDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    };
}

omnetpp::any_ptr MergeAckMsgDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    MergeAckMsg *pp = omnetpp::fromAnyPtr<MergeAckMsg>(object); (void)pp;
    switch (field) {
        default: return omnetpp::any_ptr(nullptr);
    }
}

void MergeAckMsgDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeAckMsg *pp = omnetpp::fromAnyPtr<MergeAckMsg>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'MergeAckMsg'", field);
    }
}

Register_Class(MergeNackMsg)

MergeNackMsg::MergeNackMsg(const char *name, short kind) : ::omnetpp::cMessage(name, kind)
{
}

MergeNackMsg::MergeNackMsg(const MergeNackMsg& other) : ::omnetpp::cMessage(other)
{
    copy(other);
}

MergeNackMsg::~MergeNackMsg()
{
}

MergeNackMsg& MergeNackMsg::operator=(const MergeNackMsg& other)
{
    if (this == &other) return *this;
    ::omnetpp::cMessage::operator=(other);
    copy(other);
    return *this;
}

void MergeNackMsg::copy(const MergeNackMsg& other)
{
    this->sender = other.sender;
}

void MergeNackMsg::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::omnetpp::cMessage::parsimPack(b);
    doParsimPacking(b,this->sender);
}

void MergeNackMsg::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::omnetpp::cMessage::parsimUnpack(b);
    doParsimUnpacking(b,this->sender);
}

const char * MergeNackMsg::getSender() const
{
    return this->sender.c_str();
}

void MergeNackMsg::setSender(const char * sender)
{
    this->sender = sender;
}

class MergeNackMsgDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_sender,
    };
  public:
    MergeNackMsgDescriptor();
    virtual ~MergeNackMsgDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(MergeNackMsgDescriptor)

MergeNackMsgDescriptor::MergeNackMsgDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(MergeNackMsg)), "omnetpp::cMessage")
{
    propertyNames = nullptr;
}

MergeNackMsgDescriptor::~MergeNackMsgDescriptor()
{
    delete[] propertyNames;
}

bool MergeNackMsgDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<MergeNackMsg *>(obj)!=nullptr;
}

const char **MergeNackMsgDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *MergeNackMsgDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int MergeNackMsgDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 1+base->getFieldCount() : 1;
}

unsigned int MergeNackMsgDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,    // FIELD_sender
    };
    return (field >= 0 && field < 1) ? fieldTypeFlags[field] : 0;
}

const char *MergeNackMsgDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "sender",
    };
    return (field >= 0 && field < 1) ? fieldNames[field] : nullptr;
}

int MergeNackMsgDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "sender") == 0) return baseIndex + 0;
    return base ? base->findField(fieldName) : -1;
}

const char *MergeNackMsgDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "string",    // FIELD_sender
    };
    return (field >= 0 && field < 1) ? fieldTypeStrings[field] : nullptr;
}

const char **MergeNackMsgDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *MergeNackMsgDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int MergeNackMsgDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    MergeNackMsg *pp = omnetpp::fromAnyPtr<MergeNackMsg>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void MergeNackMsgDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeNackMsg *pp = omnetpp::fromAnyPtr<MergeNackMsg>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'MergeNackMsg'", field);
    }
}

const char *MergeNackMsgDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    MergeNackMsg *pp = omnetpp::fromAnyPtr<MergeNackMsg>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string MergeNackMsgDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    MergeNackMsg *pp = omnetpp::fromAnyPtr<MergeNackMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: return oppstring2string(pp->getSender());
        default: return "";
    }
}

void MergeNackMsgDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeNackMsg *pp = omnetpp::fromAnyPtr<MergeNackMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: pp->setSender((value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'MergeNackMsg'", field);
    }
}

omnetpp::cValue MergeNackMsgDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    MergeNackMsg *pp = omnetpp::fromAnyPtr<MergeNackMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: return pp->getSender();
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'MergeNackMsg' as cValue -- field index out of range?", field);
    }
}

void MergeNackMsgDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeNackMsg *pp = omnetpp::fromAnyPtr<MergeNackMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: pp->setSender(value.stringValue()); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'MergeNackMsg'", field);
    }
}

const char *MergeNackMsgDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    };
}

omnetpp::any_ptr MergeNackMsgDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    MergeNackMsg *pp = omnetpp::fromAnyPtr<MergeNackMsg>(object); (void)pp;
    switch (field) {
        default: return omnetpp::any_ptr(nullptr);
    }
}

void MergeNackMsgDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeNackMsg *pp = omnetpp::fromAnyPtr<MergeNackMsg>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'MergeNackMsg'", field);
    }
}

Register_Class(MergeAcceptMsg)

MergeAcceptMsg::MergeAcceptMsg(const char *name, short kind) : ::omnetpp::cMessage(name, kind)
{
}

MergeAcceptMsg::MergeAcceptMsg(const MergeAcceptMsg& other) : ::omnetpp::cMessage(other)
{
    copy(other);
}

MergeAcceptMsg::~MergeAcceptMsg()
{
}

MergeAcceptMsg& MergeAcceptMsg::operator=(const MergeAcceptMsg& other)
{
    if (this == &other) return *this;
    ::omnetpp::cMessage::operator=(other);
    copy(other);
    return *this;
}

void MergeAcceptMsg::copy(const MergeAcceptMsg& other)
{
    this->sender = other.sender;
}

void MergeAcceptMsg::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::omnetpp::cMessage::parsimPack(b);
    doParsimPacking(b,this->sender);
}

void MergeAcceptMsg::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::omnetpp::cMessage::parsimUnpack(b);
    doParsimUnpacking(b,this->sender);
}

const char * MergeAcceptMsg::getSender() const
{
    return this->sender.c_str();
}

void MergeAcceptMsg::setSender(const char * sender)
{
    this->sender = sender;
}

class MergeAcceptMsgDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_sender,
    };
  public:
    MergeAcceptMsgDescriptor();
    virtual ~MergeAcceptMsgDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(MergeAcceptMsgDescriptor)

MergeAcceptMsgDescriptor::MergeAcceptMsgDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(MergeAcceptMsg)), "omnetpp::cMessage")
{
    propertyNames = nullptr;
}

MergeAcceptMsgDescriptor::~MergeAcceptMsgDescriptor()
{
    delete[] propertyNames;
}

bool MergeAcceptMsgDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<MergeAcceptMsg *>(obj)!=nullptr;
}

const char **MergeAcceptMsgDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *MergeAcceptMsgDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int MergeAcceptMsgDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 1+base->getFieldCount() : 1;
}

unsigned int MergeAcceptMsgDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,    // FIELD_sender
    };
    return (field >= 0 && field < 1) ? fieldTypeFlags[field] : 0;
}

const char *MergeAcceptMsgDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "sender",
    };
    return (field >= 0 && field < 1) ? fieldNames[field] : nullptr;
}

int MergeAcceptMsgDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "sender") == 0) return baseIndex + 0;
    return base ? base->findField(fieldName) : -1;
}

const char *MergeAcceptMsgDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "string",    // FIELD_sender
    };
    return (field >= 0 && field < 1) ? fieldTypeStrings[field] : nullptr;
}

const char **MergeAcceptMsgDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *MergeAcceptMsgDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int MergeAcceptMsgDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    MergeAcceptMsg *pp = omnetpp::fromAnyPtr<MergeAcceptMsg>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void MergeAcceptMsgDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeAcceptMsg *pp = omnetpp::fromAnyPtr<MergeAcceptMsg>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'MergeAcceptMsg'", field);
    }
}

const char *MergeAcceptMsgDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    MergeAcceptMsg *pp = omnetpp::fromAnyPtr<MergeAcceptMsg>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string MergeAcceptMsgDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    MergeAcceptMsg *pp = omnetpp::fromAnyPtr<MergeAcceptMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: return oppstring2string(pp->getSender());
        default: return "";
    }
}

void MergeAcceptMsgDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeAcceptMsg *pp = omnetpp::fromAnyPtr<MergeAcceptMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: pp->setSender((value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'MergeAcceptMsg'", field);
    }
}

omnetpp::cValue MergeAcceptMsgDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    MergeAcceptMsg *pp = omnetpp::fromAnyPtr<MergeAcceptMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: return pp->getSender();
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'MergeAcceptMsg' as cValue -- field index out of range?", field);
    }
}

void MergeAcceptMsgDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeAcceptMsg *pp = omnetpp::fromAnyPtr<MergeAcceptMsg>(object); (void)pp;
    switch (field) {
        case FIELD_sender: pp->setSender(value.stringValue()); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'MergeAcceptMsg'", field);
    }
}

const char *MergeAcceptMsgDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    };
}

omnetpp::any_ptr MergeAcceptMsgDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    MergeAcceptMsg *pp = omnetpp::fromAnyPtr<MergeAcceptMsg>(object); (void)pp;
    switch (field) {
        default: return omnetpp::any_ptr(nullptr);
    }
}

void MergeAcceptMsgDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    MergeAcceptMsg *pp = omnetpp::fromAnyPtr<MergeAcceptMsg>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'MergeAcceptMsg'", field);
    }
}

namespace omnetpp {

}  // namespace omnetpp

