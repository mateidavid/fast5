#ifndef __FAST5_PACK_HPP
#define __FAST5_PACK_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <limits>
#include <stdexcept>
#include <cassert>

//#include <bitset>
//#include "logger.hpp"

namespace fast5_pack
{
    class Huffman_Coder
    {
    public:
        typedef std::vector< std::uint8_t> Code_Type;
        typedef std::map< std::string, std::string > Code_Params_Type;

        Huffman_Coder() = default;
        Huffman_Coder(Huffman_Coder const &) = delete;
        Huffman_Coder & operator = (Huffman_Coder const &) = delete;
        Huffman_Coder(std::istream & is, std::string const & cwm_name)
        {
            load_codeword_map(is, cwm_name);
        }
        Huffman_Coder(std::vector< std::string > const & v, std::string const & cwm_name)
        {
            load_codeword_map(v, cwm_name);
        }

        void load_codeword_map(std::istream & is, std::string const & cwm_name)
        {
            _cwm_name = cwm_name;
            std::string v_s;
            std::string cw_s;
            while (is >> v_s >> cw_s)
            {
                add_codeword(v_s, cw_s);
            }
        }
        void load_codeword_map(std::vector< std::string > const & v, std::string const & cwm_name)
        {
            _cwm_name = cwm_name;
            for (unsigned i = 0; i < v.size() - 1; i += 2)
            {
                add_codeword(v[i], v[i + 1]);
            }
        }

        template < typename Int_Type >
        std::pair< Code_Type, Code_Params_Type >
        encode(std::vector< Int_Type > const & v, bool encode_diff = false)
        {
            Code_Type res;
            Code_Params_Type res_params = id();
            res_params["code_diff"] = encode_diff? "1" : "0";
            uint64_t buff = 0;
            uint8_t buff_len = 0;
            bool reset = true;
            Int_Type last = 0;
            unsigned i = 0;
            long long int val;
            long long int x;
            while (true)
            {
                assert(buff_len <= 64);
                // flush buffer
                while (buff_len >= 8)
                {
                    res.push_back(buff & 0xFF);
                    buff >>= 8;
                    buff_len -= 8;
                }
                assert(buff_len < 8);
                if (reset)
                {
                    assert(buff_len == 0);
                    if (i == v.size()) break;
                    //LOG(debug) << "absolute value val=" << v[i] << std::endl;
                    for (unsigned j = 0; j < sizeof(Int_Type); ++j)
                    {
                        std::uint8_t y = (v[i] >> (8 * j)) & 0xFF;
                        //LOG(debug) << "byte " << j << ": " << std::bitset<8>(y) << std::endl;
                        res.push_back(y);
                    }
                    reset = false;
                    last = v[i];
                    ++i;
                }
                else // not reset
                {
                    if (i < v.size())
                    {
                        val = v[i];
                        x = encode_diff? val - last : val;
                        reset = _cwm.count(x) == 0;
                        //LOG(debug) << "relative value: val=" << v[i] << " last=" << last << " x=" << x << " reset=" << reset << std::endl;
                    }
                    else
                    {
                        reset = true;
                        //LOG(debug) << "end: reset=1" << std::endl;
                    }
                    auto p = (not reset? _cwm[x] : _cwm[break_cw()]);
                    buff |= (p.first << buff_len);
                    buff_len += p.second;
                    if (not reset)
                    {
                        last = v[i];
                        ++i;
                    }
                    else if ((buff_len % 8) > 0) // and reset
                    {
                        buff_len += 8 - (buff_len % 8);
                    }
                        
                }
            }
            return std::make_pair(std::move(res), std::move(res_params));
        }

        template < typename Int_Type >
        std::vector< Int_Type >
        decode(Code_Type const & v, Code_Params_Type const & v_params)
        {
            check_params(v_params);
            bool decode_diff = v_params.at("code_diff") == "1";
            std::vector< Int_Type > res;
            std::uint64_t buff = 0;
            std::uint8_t buff_len = 0;
            bool reset = true;
            Int_Type last = 0;
            unsigned i = 0;
            while (i < v.size() or buff_len > 0)
            {
                assert(buff_len <= 64);
                // fill buffer
                while (i < v.size() and buff_len <= 56)
                {
                    uint64_t y = v[i];
                    buff |= (y << buff_len);
                    buff_len += 8;
                    ++i;
                }
                assert(buff_len <= 64);
                if (reset)
                {
                    assert((buff_len % 8) == 0);
                    assert(buff_len / 8 >= sizeof(Int_Type));
                    //LOG(debug) << "absolute value" << std::endl;
                    Int_Type x = 0;
                    for (unsigned j = 0; j < sizeof(Int_Type); ++j)
                    {
                        std::uint64_t y = (buff & 0xFF);
                        //LOG(debug) << "byte " << j << ": " << std::bitset<8>(y) << std::endl;
                        x |= (y << (8 * j));
                        buff >>= 8;
                        buff_len -= 8;
                    }
                    //LOG(debug) << "got: val=" << x << std::endl;
                    res.push_back(x);
                    last = x;
                    reset = false;
                }
                else // not reset
                {
                    //LOG(debug) << "reading relative value" << std::endl;
                    // TODO: faster decoding
                    // currently, try all codewords one by one
                    auto it = _cwm.begin();
                    while (it != _cwm.end())
                    {
                        if ((buff & ((1llu << it->second.second) - 1)) == it->second.first)
                        {
                            break;
                        }
                        ++it;
                    }
                    if (it == _cwm.end()) throw std::invalid_argument("decoding failure: codeword not found");
                    auto x = it->first;
                    auto p = it->second;
                    assert(buff_len >= p.second);
                    buff >>= p.second;
                    buff_len -= p.second;
                    if (x != break_cw())
                    {
                        //LOG(debug) << "got: x=" << x << " last=" << last << " val=" << x + last << " cw_len=" << (int)p.second << std::endl;
                        if (decode_diff) x += last;
                        if (sizeof(Int_Type) < 8
                            and (x < (long long)std::numeric_limits< Int_Type >::min()
                                 or x > (long long)std::numeric_limits< Int_Type >::max()))
                        {
                            throw std::invalid_argument("decoding failure: overflow");
                        }
                        res.push_back(x);
                        last = x;
                    }
                    else
                    {
                        //LOG(debug) << "got: break cw_len=" << (int)p.second << std::endl;
                        reset = true;
                        if ((buff_len % 8) > 0)
                        {
                            buff >>= (buff_len % 8);
                            buff_len -= (buff_len % 8);
                        }
                    }
                }
            }
            return res;
        }
    private:
        std::map< long long int, std::pair< std::uint64_t, std::uint8_t > > _cwm;
        std::string _cwm_name;
        static long long int break_cw()
        {
            static long long int const _break_cw = std::numeric_limits< long long int >::min();
            return _break_cw;
        }
        Code_Params_Type id() const
        {
            Code_Params_Type res;
            res["packer"] = "huffman_coder";
            res["format_version"] = "2";
            res["codeword_map_name"] = _cwm_name;
            return res;
        }
        void check_params(Code_Params_Type const & params) const
        {
            auto _id = id();
            if (params.at("packer") != _id.at("packer")
                or params.at("format_version") != _id.at("format_version")
                or params.at("codeword_map_name") != _id.at("codeword_map_name"))
            {
                throw std::invalid_argument("decode id mismatch");
            }
        }
        void add_codeword(std::string const & v_s, std::string const & cw_s)
        {
            long long int v;
            if (v_s != ".")
            {
                std::istringstream(v_s) >> v;
            }
            else
            {
                v = break_cw();
            }
            std::uint64_t cw = 0;
            if (cw_s.size() > 57)
            {
                throw std::invalid_argument(std::string("codeword too long: ") + v_s + " " + cw_s);
            }
            std::uint8_t cw_l = cw_s.size();
            for (int i = cw_s.size() - 1; i >= 0; --i)
            {
                cw <<= 1;
                cw |= (cw_s[i] == '1');
            }
            _cwm[v] = std::make_pair(cw, cw_l);
        }
    }; // class Huffman_Coder
} // fast5_pack

#endif
