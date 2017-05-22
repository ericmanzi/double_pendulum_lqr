/** THIS IS AN AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY
 * BY HAND!!
 *
 * Generated by lcm-gen
 **/

#include <lcm/lcm_coretypes.h>

#ifndef __drake_lcmt_viewer_link_data_hpp__
#define __drake_lcmt_viewer_link_data_hpp__

#include <string>
#include <vector>
#include "lcmtypes/drake/lcmt_viewer_geometry_data.hpp"

namespace drake
{

class lcmt_viewer_link_data
{
    public:
        std::string name;
        int32_t    robot_num;
        int32_t    num_geom;
        std::vector< drake::lcmt_viewer_geometry_data > geom;

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

int lcmt_viewer_link_data::encode(void *buf, int offset, int maxlen) const
{
    int pos = 0, tlen;
    int64_t hash = getHash();

    tlen = __int64_t_encode_array(buf, offset + pos, maxlen - pos, &hash, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->_encodeNoHash(buf, offset + pos, maxlen - pos);
    if (tlen < 0) return tlen; else pos += tlen;

    return pos;
}

int lcmt_viewer_link_data::decode(const void *buf, int offset, int maxlen)
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

int lcmt_viewer_link_data::getEncodedSize() const
{
    return 8 + _getEncodedSizeNoHash();
}

int64_t lcmt_viewer_link_data::getHash()
{
    static int64_t hash = _computeHash(NULL);
    return hash;
}

const char* lcmt_viewer_link_data::getTypeName()
{
    return "lcmt_viewer_link_data";
}

int lcmt_viewer_link_data::_encodeNoHash(void *buf, int offset, int maxlen) const
{
    int pos = 0, tlen;

    char* name_cstr = (char*) this->name.c_str();
    tlen = __string_encode_array(buf, offset + pos, maxlen - pos, &name_cstr, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->robot_num, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->num_geom, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    for (int a0 = 0; a0 < this->num_geom; a0++) {
        tlen = this->geom[a0]._encodeNoHash(buf, offset + pos, maxlen - pos);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    return pos;
}

int lcmt_viewer_link_data::_decodeNoHash(const void *buf, int offset, int maxlen)
{
    int pos = 0, tlen;

    int32_t __name_len__;
    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &__name_len__, 1);
    if(tlen < 0) return tlen; else pos += tlen;
    if(__name_len__ > maxlen - pos) return -1;
    this->name.assign(((const char*)buf) + offset + pos, __name_len__ - 1);
    pos += __name_len__;

    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->robot_num, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->num_geom, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    this->geom.resize(this->num_geom);
    for (int a0 = 0; a0 < this->num_geom; a0++) {
        tlen = this->geom[a0]._decodeNoHash(buf, offset + pos, maxlen - pos);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    return pos;
}

int lcmt_viewer_link_data::_getEncodedSizeNoHash() const
{
    int enc_size = 0;
    enc_size += this->name.size() + 4 + 1;
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    for (int a0 = 0; a0 < this->num_geom; a0++) {
        enc_size += this->geom[a0]._getEncodedSizeNoHash();
    }
    return enc_size;
}

int64_t lcmt_viewer_link_data::_computeHash(const __lcm_hash_ptr *p)
{
    const __lcm_hash_ptr *fp;
    for(fp = p; fp != NULL; fp = fp->parent)
        if(fp->v == lcmt_viewer_link_data::getHash)
            return 0;
    const __lcm_hash_ptr cp = { p, (void*)lcmt_viewer_link_data::getHash };

    int64_t hash = 0x4b645ec7a5743a2aLL +
         drake::lcmt_viewer_geometry_data::_computeHash(&cp);

    return (hash<<1) + ((hash>>63)&1);
}

}

#endif
