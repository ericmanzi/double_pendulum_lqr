"""LCM type definitions
This file automatically generated by lcm.
DO NOT MODIFY BY HAND!!!!
"""

import cStringIO as StringIO
import struct

import lcmt_viewer_geometry_data

class lcmt_viewer_link_data(object):
    __slots__ = ["name", "robot_num", "num_geom", "geom"]

    def __init__(self):
        self.name = ""
        self.robot_num = 0
        self.num_geom = 0
        self.geom = []

    def encode(self):
        buf = StringIO.StringIO()
        buf.write(lcmt_viewer_link_data._get_packed_fingerprint())
        self._encode_one(buf)
        return buf.getvalue()

    def _encode_one(self, buf):
        __name_encoded = self.name.encode('utf-8')
        buf.write(struct.pack('>I', len(__name_encoded)+1))
        buf.write(__name_encoded)
        buf.write("\0")
        buf.write(struct.pack(">ii", self.robot_num, self.num_geom))
        for i0 in range(self.num_geom):
            assert self.geom[i0]._get_packed_fingerprint() == lcmt_viewer_geometry_data.lcmt_viewer_geometry_data._get_packed_fingerprint()
            self.geom[i0]._encode_one(buf)

    def decode(data):
        if hasattr(data, 'read'):
            buf = data
        else:
            buf = StringIO.StringIO(data)
        if buf.read(8) != lcmt_viewer_link_data._get_packed_fingerprint():
            raise ValueError("Decode error")
        return lcmt_viewer_link_data._decode_one(buf)
    decode = staticmethod(decode)

    def _decode_one(buf):
        self = lcmt_viewer_link_data()
        __name_len = struct.unpack('>I', buf.read(4))[0]
        self.name = buf.read(__name_len)[:-1].decode('utf-8', 'replace')
        self.robot_num, self.num_geom = struct.unpack(">ii", buf.read(8))
        self.geom = []
        for i0 in range(self.num_geom):
            self.geom.append(lcmt_viewer_geometry_data.lcmt_viewer_geometry_data._decode_one(buf))
        return self
    _decode_one = staticmethod(_decode_one)

    _hash = None
    def _get_hash_recursive(parents):
        if lcmt_viewer_link_data in parents: return 0
        newparents = parents + [lcmt_viewer_link_data]
        tmphash = (0x4b645ec7a5743a2a+ lcmt_viewer_geometry_data.lcmt_viewer_geometry_data._get_hash_recursive(newparents)) & 0xffffffffffffffff
        tmphash  = (((tmphash<<1)&0xffffffffffffffff)  + (tmphash>>63)) & 0xffffffffffffffff
        return tmphash
    _get_hash_recursive = staticmethod(_get_hash_recursive)
    _packed_fingerprint = None

    def _get_packed_fingerprint():
        if lcmt_viewer_link_data._packed_fingerprint is None:
            lcmt_viewer_link_data._packed_fingerprint = struct.pack(">Q", lcmt_viewer_link_data._get_hash_recursive([]))
        return lcmt_viewer_link_data._packed_fingerprint
    _get_packed_fingerprint = staticmethod(_get_packed_fingerprint)

