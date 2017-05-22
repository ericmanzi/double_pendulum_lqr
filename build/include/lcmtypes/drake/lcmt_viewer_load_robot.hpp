/** THIS IS AN AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY
 * BY HAND!!
 *
 * Generated by lcm-gen
 **/

#include <lcm/lcm_coretypes.h>

#ifndef __drake_lcmt_viewer_load_robot_hpp__
#define __drake_lcmt_viewer_load_robot_hpp__

#include <vector>
#include "lcmtypes/drake/lcmt_viewer_link_data.hpp"

namespace drake
{

class lcmt_viewer_load_robot
{
    public:
        int32_t    num_links;
        std::vector< drake::lcmt_viewer_link_data > link;

    public:
        inline int encode(void *buf, int offset, int maxlen) const;
        inline int getEncodedSize() const;
        inline int decode(const void *buf, int offset, int maxlen);
        inline static int64_t getHash();
        inline static const char* getTypeName();

        // LCM support functions. Users should not call these
        inline int _encodeNoHash(void *buf, int offset, int maxlen) const;
        inline int _getEncodedSizeNoHash() const;
        inline int _decodeNoHash(const void *buf, int offset, int maxlen);
        inline static int64_t _computeHash(const __lcm_hash_ptr *p);
};

int lcmt_viewer_load_robot::encode(void *buf, int offset, int maxlen) const
{
    int pos = 0, tlen;
    int64_t hash = getHash();

    tlen = __int64_t_encode_array(buf, offset + pos, maxlen - pos, &hash, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->_encodeNoHash(buf, offset + pos, maxlen - pos);
    if (tlen < 0) return tlen; else pos += tlen;

    return pos;
}

int lcmt_viewer_load_robot::decode(const void *buf, int offset, int maxlen)
{
    int pos = 0, thislen;

    int64_t msg_hash;
    thislen = __int64_t_decode_array(buf, offset + pos, maxlen - pos, &msg_hash, 1);
    if (thislen < 0) return thislen; else pos += thislen;
    if (msg_hash != getHash()) return -1;

    thislen = this->_decodeNoHash(buf, offset + pos, maxlen - pos);
    if (thislen < 0) return thislen; else pos += thislen;

    return pos;
}

int lcmt_viewer_load_robot::getEncodedSize() const
{
    return 8 + _getEncodedSizeNoHash();
}

int64_t lcmt_viewer_load_robot::getHash()
{
    static int64_t hash = _computeHash(NULL);
    return hash;
}

const char* lcmt_viewer_load_robot::getTypeName()
{
    return "lcmt_viewer_load_robot";
}

int lcmt_viewer_load_robot::_encodeNoHash(void *buf, int offset, int maxlen) const
{
    int pos = 0, tlen;

    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->num_links, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    for (int a0 = 0; a0 < this->num_links; a0++) {
        tlen = this->link[a0]._encodeNoHash(buf, offset + pos, maxlen - pos);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    return pos;
}

int lcmt_viewer_load_robot::_decodeNoHash(const void *buf, int offset, int maxlen)
{
    int pos = 0, tlen;

    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->num_links, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    this->link.resize(this->num_links);
    for (int a0 = 0; a0 < this->num_links; a0++) {
        tlen = this->link[a0]._decodeNoHash(buf, offset + pos, maxlen - pos);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    return pos;
}

int lcmt_viewer_load_robot::_getEncodedSizeNoHash() const
{
    int enc_size = 0;
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    for (int a0 = 0; a0 < this->num_links; a0++) {
        enc_size += this->link[a0]._getEncodedSizeNoHash();
    }
    return enc_size;
}

int64_t lcmt_viewer_load_robot::_computeHash(const __lcm_hash_ptr *p)
{
    const __lcm_hash_ptr *fp;
    for(fp = p; fp != NULL; fp = fp->parent)
        if(fp->v == lcmt_viewer_load_robot::getHash)
            return 0;
    const __lcm_hash_ptr cp = { p, (void*)lcmt_viewer_load_robot::getHash };

    int64_t hash = 0x739e6927d8bcec39LL +
         drake::lcmt_viewer_link_data::_computeHash(&cp);

    return (hash<<1) + ((hash>>63)&1);
}

}

#endif
