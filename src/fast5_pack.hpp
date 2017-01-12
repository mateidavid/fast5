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
    class Huffman_Diff_Coder
    {
    public:
        Huffman_Diff_Coder() = default;
        Huffman_Diff_Coder(Huffman_Diff_Coder const &) = delete;
        Huffman_Diff_Coder & operator = (Huffman_Diff_Coder const &) = delete;
        Huffman_Diff_Coder(std::istream & cw_is, std::string const & cw_m_name)
        {
            load_codeword_map(cw_is, cw_m_name);
        }
        Huffman_Diff_Coder(std::vector< std::string > const & cw_v, std::string const & cw_m_name)
        {
            load_codeword_map(cw_v, cw_m_name);
        }

        void load_codeword_map(std::istream & cw_is, std::string const & cw_m_name)
        {
            _cw_m_name = cw_m_name;
            std::string v_s;
            std::string cw_s;
            while (cw_is >> v_s >> cw_s)
            {
                add_codeword(v_s, cw_s);
            }
        }
        void load_codeword_map(std::vector< std::string > const & cw_v, std::string const & cw_m_name)
        {
            _cw_m_name = cw_m_name;
            for (unsigned i = 0; i < cw_v.size() - 1; i += 2)
            {
                add_codeword(cw_v[i], cw_v[i + 1]);
            }
        }

        template < typename Int_Type >
        std::pair< std::vector< uint8_t >, std::map< std::string, std::string > >
        encode(std::vector< Int_Type > const & v)
        {
            std::vector< uint8_t > res;
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
                while (buff_len >= 8)
                {
                    res.push_back(buff & 0xFF);
                    buff >>= 8;
                    buff_len -= 8;
                }
                assert(buff_len < 8);
                if (reset)
                {
                    if (i == v.size()) break;
                    assert(buff_len == 0);
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
                        x = val - last;
                        reset = _cw_m.count(x) == 0;
                        //LOG(debug) << "relative value: val=" << v[i] << " last=" << last << " x=" << x << " reset=" << reset << std::endl;
                    }
                    else
                    {
                        reset = true;
                        //LOG(debug) << "end: reset=1" << std::endl;
                    }
                    auto p = (not reset? _cw_m[x] : _cw_m[break_cw()]);
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
            return std::make_pair(std::move(res), id());
        }

        template < typename Int_Type >
        std::vector< Int_Type > decode(std::vector< uint8_t > const & v,
                                       std::map< std::string, std::string > const & v_id)
        {
            if (v_id != id()) throw std::invalid_argument("decode id mismatch");
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
                    auto it = _cw_m.begin();
                    while (it != _cw_m.end())
                    {
                        if ((buff & ((1llu << it->second.second) - 1)) == it->second.first)
                        {
                            break;
                        }
                        ++it;
                    }
                    if (it == _cw_m.end()) throw std::invalid_argument("decoding failure");
                    auto x = it->first;
                    auto p = it->second;
                    assert(buff_len >= p.second);
                    buff >>= p.second;
                    buff_len -= p.second;
                    if (x != break_cw())
                    {
                        //LOG(debug) << "got: x=" << x << " last=" << last << " val=" << x + last << " cw_len=" << (int)p.second << std::endl;
                        x += last;
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
        std::map< long long int, std::pair< std::uint64_t, std::uint8_t > > _cw_m;
        std::string _cw_m_name;
        static long long int break_cw()
        {
            static long long int const _break_cw = std::numeric_limits< long long int >::min();
            return _break_cw;
        }
        std::map< std::string, std::string > id() const
        {
            std::map< std::string, std::string > res;
            res["packer"] = "huffman_diff_coder";
            res["format_version"] = "1";
            res["codeword_map_name"] = _cw_m_name;
            return res;
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
            _cw_m[v] = std::make_pair(cw, cw_l);
        }
    }; // class Huffman_Diff_Coder
}

#endif
