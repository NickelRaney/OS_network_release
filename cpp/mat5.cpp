#include "mat5.hpp"
#include "third_party/miniz/miniz_tinfl.h"
#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <utility>

namespace mat5 {
namespace {
constexpr std::uint32_t miINT8=1, miUINT8=2, miINT16=3, miUINT16=4,
    miINT32=5, miUINT32=6, miSINGLE=7, miDOUBLE=9, miINT64=12,
    miUINT64=13, miMATRIX=14, miCOMPRESSED=15;
constexpr std::uint32_t mxSTRUCT=2, mxSPARSE=5, mxDOUBLE=6;
[[noreturn]] void fail(const std::string& reason) {
    throw std::runtime_error("MAT-v5: " + reason);
}
std::size_t add(std::size_t a, std::size_t b) {
    if (b > std::numeric_limits<std::size_t>::max()-a) fail("size overflow");
    return a+b;
}
std::size_t mul(std::size_t a, std::size_t b) {
    if (a && b > std::numeric_limits<std::size_t>::max()/a) fail("size overflow");
    return a*b;
}
std::size_t padded(std::size_t n) { return add(n,7)&~std::size_t(7); }
std::uint32_t u32(const std::uint8_t* p) {
    return std::uint32_t(p[0]) | std::uint32_t(p[1])<<8 |
           std::uint32_t(p[2])<<16 | std::uint32_t(p[3])<<24;
}
struct Element { std::uint32_t type; const std::uint8_t* data; std::size_t size; };
struct Cursor {
    const std::uint8_t* data;
    std::size_t size, pos=0;
    Element next() {
        if (size-pos<8) fail("truncated data-element tag");
        const auto* tag=data+pos;
        std::uint32_t word=u32(tag), type=word;
        std::size_t count;
        const std::uint8_t* payload;
        if (word>>16) {
            type=word&0xffff; count=word>>16;
            if (count>4) fail("invalid small data element");
            payload=tag+4; pos+=8;
        } else {
            count=u32(tag+4); payload=tag+8;
            auto total=add(8,type==miCOMPRESSED?count:padded(count));
            if (total>size-pos) fail("data element exceeds enclosing matrix");
            pos+=total;
        }
        return {type,payload,count};
    }
};
std::string nameOf(const Element& e) {
    if (e.type!=miINT8 && e.type!=miUINT8) fail("unsupported variable-name encoding");
    return std::string(reinterpret_cast<const char*>(e.data),e.size);
}
std::size_t typeWidth(std::uint32_t type) {
    switch(type) {
    case miINT8: case miUINT8:return 1;
    case miINT16:case miUINT16:return 2;
    case miINT32:case miUINT32:case miSINGLE:return 4;
    case miDOUBLE:case miINT64:case miUINT64:return 8;
    default:fail("unsupported numeric storage type " + std::to_string(type));
    }
}
template<class T> double number(const std::uint8_t* p) {
    T value; std::memcpy(&value,p,sizeof value); return static_cast<double>(value);
}
double number(const Element& e,std::size_t i) {
    const auto* p=e.data+i*typeWidth(e.type);
    switch(e.type) {
    case miINT8:return number<std::int8_t>(p);
    case miUINT8:return number<std::uint8_t>(p);
    case miINT16:return number<std::int16_t>(p);
    case miUINT16:return number<std::uint16_t>(p);
    case miINT32:return number<std::int32_t>(p);
    case miUINT32:return number<std::uint32_t>(p);
    case miSINGLE:return number<float>(p);
    case miDOUBLE:return number<double>(p);
    case miINT64:return number<std::int64_t>(p);
    case miUINT64:return number<std::uint64_t>(p);
    default:fail("unsupported numeric storage type");
    }
}
void numeric(Node& node,const Element& data,std::size_t count) {
    if (mul(count,typeWidth(data.type))!=data.size) fail("numeric dimensions and byte count disagree");
    if(node.logical) {
        node.bytes.resize(count);
        if(data.type==miUINT8 || data.type==miINT8) {
            // Logical connectivity is by far the largest input. Keep this
            // simple loop available for the compiler's byte-vectorization.
            for(std::size_t i=0;i<count;++i)node.bytes[i]=data.data[i]!=0;
        } else for(std::size_t i=0;i<count;++i) node.bytes[i]=number(data,i)!=0;
    } else {
        node.values.resize(count);
        if(data.type==miDOUBLE && count) std::memcpy(node.values.data(),data.data,data.size);
        else for(std::size_t i=0;i<count;++i) node.values[i]=number(data,i);
    }
}
std::vector<std::uint32_t> indices(const Element& data) {
    if((data.type!=miINT32 && data.type!=miUINT32) || data.size%4)
        fail("sparse indices must contain 32-bit integers");
    std::vector<std::uint32_t> out(data.size/4);
    for(std::size_t i=0;i<out.size();++i) {
        out[i]=u32(data.data+4*i);
        if(out[i]>static_cast<std::uint32_t>(std::numeric_limits<std::int32_t>::max())) fail("negative sparse index");
    }
    return out;
}
std::pair<std::string,Node> matrix(const Element& matrixData,unsigned depth=0) {
    if(matrixData.type!=miMATRIX) fail("expected matrix element");
    if(depth>32) fail("struct nesting exceeds supported depth");
    Cursor c{matrixData.data,matrixData.size};
    auto flags=c.next();
    if(flags.type!=miUINT32 || flags.size!=8) fail("invalid array flags");
    const auto bits=u32(flags.data), cls=bits&0xff;
    if(bits&0x0800) fail("complex arrays are unsupported");
    Node result; result.logical=(bits&0x0200)!=0;
    auto dims=c.next();
    if(dims.type!=miINT32 || dims.size%4 || dims.size<8) fail("invalid dimensions");
    for(std::size_t i=0;i<dims.size;i+=4) {
        auto n=u32(dims.data+i);
        if(n>static_cast<std::uint32_t>(std::numeric_limits<std::int32_t>::max())) fail("negative array dimension");
        result.dims.push_back(n);
    }
    std::string name=nameOf(c.next());
    if(cls==mxSTRUCT) {
        result.kind=Node::Struct;
        if(result.numel()!=1) fail("only scalar structs are supported: " + name);
        auto len=c.next();
        if(len.type!=miINT32 || len.size!=4) fail("invalid struct field-name length");
        auto width=u32(len.data);
        if(!width || width>65536) fail("invalid struct field-name width");
        auto names=c.next();
        if(names.type!=miINT8 || names.size%width) fail("invalid struct field names");
        for(std::size_t i=0;i<names.size;i+=width) {
            const char* p=reinterpret_cast<const char*>(names.data+i);
            auto zero=static_cast<const char*>(std::memchr(p,0,width));
            if(!zero) fail("unterminated struct field name");
            std::string field(p,zero);
            if(field.empty()) fail("empty struct field name");
            auto value=matrix(c.next(),depth+1);
            if(!result.fields.emplace(field,std::move(value.second)).second) fail("duplicate struct field: " + field);
        }
    } else if(cls==mxSPARSE) {
        result.kind=Node::Sparse;
        if(result.dims.size()!=2) fail("sparse matrix must be two-dimensional");
        result.ir=indices(c.next()); result.jc=indices(c.next());
        if(result.jc.size()!=result.dims[1]+1 || result.jc.front()!=0)
            fail("invalid sparse column pointers");
        for(std::size_t i=1;i<result.jc.size();++i)
            if(result.jc[i]<result.jc[i-1]) fail("sparse column pointers are not ordered");
        const std::size_t nnz=result.jc.back();
        if(nnz>result.ir.size()) fail("sparse row index array is too short");
        auto data=c.next();
        if(data.size%typeWidth(data.type) || data.size/typeWidth(data.type)<nnz)
            fail("sparse value array is too short");
        // MATLAB can save unused capacity after the final used entry.
        data.size=mul(nnz,typeWidth(data.type)); numeric(result,data,nnz);
        result.ir.resize(nnz);
        for(auto row:result.ir) if(row>=result.dims[0]) fail("sparse row index outside matrix");
    } else if(cls>=6 && cls<=15) {
        numeric(result,c.next(),result.numel());
    } else fail("unsupported MATLAB array class " + std::to_string(cls) + " in " + name);
    if(c.pos!=c.size) fail("unexpected trailing matrix data in " + name);
    return {std::move(name),std::move(result)};
}
// Prefix parsing intentionally ignores the declared matrix payload length:
// only the flags, dimensions, and variable-name header must be available.
std::string prefixName(const std::uint8_t* data,std::size_t n) {
    if(n<8 || u32(data)!=miMATRIX) fail("compressed element does not contain a matrix");
    Cursor c{data+8,n-8}; c.next(); c.next(); return nameOf(c.next());
}
std::vector<std::uint8_t> inflateSelected(const std::vector<std::uint8_t>& input,
                                        const std::set<std::string>& wanted) {
    tinfl_decompressor decompressor; tinfl_init(&decompressor);
    std::vector<std::uint8_t> output(65536);
    std::size_t inPos=0,outPos=0; bool inspected=false;
    for(;;) {
        auto inCount=input.size()-inPos, outCount=output.size()-outPos;
        auto status=tinfl_decompress(&decompressor,input.data()+inPos,&inCount,
            output.data(),output.data()+outPos,&outCount,
            TINFL_FLAG_PARSE_ZLIB_HEADER|TINFL_FLAG_USING_NON_WRAPPING_OUTPUT_BUF);
        inPos+=inCount; outPos+=outCount;
        if(status<TINFL_STATUS_DONE) fail("invalid or truncated zlib-compressed matrix");
        if(!inspected) {
            auto name=prefixName(output.data(),outPos); inspected=true;
            if(!wanted.empty() && !wanted.count(name)) return {};
        }
        if(status==TINFL_STATUS_DONE) {
            if(inPos!=input.size()) fail("trailing data after zlib stream");
            output.resize(outPos); return output;
        }
        if(status!=TINFL_STATUS_HAS_MORE_OUTPUT) fail("truncated zlib stream");
        // The outer matrix tag declares the exact decompressed size. Allocate
        // once after reading its prefix instead of repeatedly growing a heap.
        auto expected=add(8,padded(u32(output.data()+4)));
        if(expected<=outPos || expected>std::size_t(UINT32_MAX)+8)
            fail("invalid decompressed matrix size");
        output.resize(expected);
    }
}
void readBytes(std::ifstream& in,void* data,std::size_t size) {
    if(size && !in.read(static_cast<char*>(data),static_cast<std::streamsize>(size))) fail("truncated file");
}
void put32(std::ofstream& out,std::uint32_t v) {
    const char b[4]={char(v),char(v>>8),char(v>>16),char(v>>24)}; out.write(b,4);
}
std::uint32_t size32(std::size_t n) {
    // Level-5 tags have unsigned 32-bit lengths, but MATLAB's v6/v7 loader
    // supports at most 2 GiB per variable. Keep the complete outer struct
    // below that compatibility limit, rather than merely fitting the tag.
    if(n>static_cast<std::size_t>(INT32_MAX)) fail("variable exceeds MAT-v5 2 GiB compatibility limit");
    return static_cast<std::uint32_t>(n);
}
void zeros(std::ofstream& out,std::size_t n) {
    const std::array<char,128> zero{};
    while(n) { const auto chunk=std::min(n,zero.size());out.write(zero.data(),chunk);n-=chunk; }
}
void tag(std::ofstream& out,std::uint32_t type,std::size_t n) {put32(out,type);put32(out,size32(n));}
void bytes(std::ofstream& out,std::uint32_t type,const void* data,std::size_t n) {
    tag(out,type,n); if(n)out.write(static_cast<const char*>(data),static_cast<std::streamsize>(n)); zeros(out,padded(n)-n);
}
std::size_t matrixHeaderSize(const std::string& name,std::size_t rank) {
    return add(add(16,add(8,padded(mul(rank,4)))),add(8,padded(name.size())));
}
void matrixHeader(std::ofstream& out,std::uint32_t cls,const std::string& name,
                  const std::vector<std::size_t>& dims) {
    tag(out,miUINT32,8);put32(out,cls);put32(out,0);
    tag(out,miINT32,mul(dims.size(),4));
    for(auto d:dims) { if(d>INT32_MAX)fail("dimension exceeds MAT-v5 INT32_MAX");put32(out,static_cast<std::uint32_t>(d)); }
    zeros(out,padded(mul(dims.size(),4))-mul(dims.size(),4));
    bytes(out,miINT8,name.data(),name.size());
}
}

const Node& Node::field(const std::string& name) const {
    if(kind!=Struct)fail("field access on a non-struct");
    auto found=fields.find(name);if(found==fields.end())fail("missing field " + name);return found->second;
}
Node& Node::field(const std::string& name) {
    return const_cast<Node&>(static_cast<const Node&>(*this).field(name));
}
std::size_t Node::numel() const {
    if(dims.empty())return 0;
    std::size_t count=1;for(auto d:dims)count=mul(count,d);return count;
}
double Node::doubleAt(std::size_t index) const {
    if(kind==Struct)fail("numeric access on a struct");
    if(logical) {if(index>=bytes.size())fail("numeric index out of bounds");return bytes[index];}
    if(index>=values.size())fail("numeric index out of bounds");return values[index];
}
double Node::scalar() const {
    if(numel()!=1 || kind!=Numeric)fail("expected a numeric scalar");return doubleAt(0);
}
Variables read(const std::filesystem::path& path,const std::set<std::string>& wanted) {
    std::ifstream in(path,std::ios::binary); if(!in)fail("cannot open " + path.string());
    in.seekg(0,std::ios::end);auto end=in.tellg();in.seekg(0);
    if(end<128)fail("file is too short for a MAT-v5 header");
    std::array<std::uint8_t,128> header{};readBytes(in,header.data(),header.size());
    if(header[126]!='I' || header[127]!='M')fail("only little-endian MAT-v5 files are supported (not v7.3/HDF5)");
    if(header[124]!=0 || header[125]!=1)fail("unsupported MAT-file version");
    Variables result;
    while(in.tellg()<end) {
        std::array<std::uint8_t,8> rawTag{};readBytes(in,rawTag.data(),8);
        const auto type=u32(rawTag.data()), size=u32(rawTag.data()+4);
        if(type!=miMATRIX && type!=miCOMPRESSED)fail("unsupported top-level data element " + std::to_string(type));
        const auto remaining=static_cast<std::uint64_t>(end-in.tellg());
        if(size>remaining)fail("top-level element exceeds file length");
        std::vector<std::uint8_t> raw;
        if(type==miCOMPRESSED) {
            raw.resize(size);readBytes(in,raw.data(),raw.size());
            raw=inflateSelected(raw,wanted);if(raw.empty())continue;
        } else {
            const auto beginning=in.tellg();
            const std::size_t prefix=std::min<std::size_t>(size,65528);
            raw.resize(prefix+8);std::copy(rawTag.begin(),rawTag.end(),raw.begin());readBytes(in,raw.data()+8,prefix);
            const auto name=prefixName(raw.data(),raw.size());
            if(!wanted.empty() && !wanted.count(name)) {in.seekg(beginning+static_cast<std::streamoff>(padded(size)));continue;}
            raw.resize(add(8,padded(size)));readBytes(in,raw.data()+8+prefix,padded(size)-prefix);
        }
        Cursor c{raw.data(),raw.size()};auto value=matrix(c.next());
        if(c.pos!=c.size)fail("unexpected trailing top-level matrix data");
        if(!result.emplace(value.first,std::move(value.second)).second)fail("duplicate top-level variable " + value.first);
        if(!wanted.empty() && result.size()==wanted.size())break;
    }
    for(const auto& name:wanted)if(!result.count(name))fail("requested variable was not found: " + name);
    return result;
}
void writeStruct(const std::filesystem::path& path,const std::string& variableName,
                 const std::vector<FieldView>& fields) {
    std::size_t width=1;
    std::set<std::string> seen;
    for(const auto& field:fields) {
        if(field.name.empty() || field.name.find('\0')!=std::string::npos || !seen.insert(field.name).second)
            fail("invalid or duplicate output field name");
        width=std::max(width,field.name.size()+1);
    }
    std::size_t payload=matrixHeaderSize(variableName,2);
    payload=add(payload,16); // Ordinary 4-byte field-name-length element plus pad.
    payload=add(payload,add(8,padded(mul(width,fields.size()))));
    std::vector<std::size_t> childSizes;
    for(const auto& field:fields) {
        if(field.dims.size()<2)fail("output array needs at least two dimensions");
        std::size_t count=1;for(auto d:field.dims) {if(d>INT32_MAX)fail("output dimension too large");count=mul(count,d);}
        if(count!=field.count || (count && !field.data))fail("output field dimensions and data length disagree");
        auto child=add(matrixHeaderSize("",field.dims.size()),add(8,mul(count,8)));
        size32(child);childSizes.push_back(child);payload=add(payload,add(8,padded(child)));
    }
    size32(payload);
    std::ofstream out(path,std::ios::binary|std::ios::trunc);if(!out)fail("cannot create " + path.string());
    std::array<char,128> header{};
    const char* text="MATLAB 5.0 MAT-file, Platform: standalone C++17, Created by V1 standalone simulator";
    std::fill(header.begin(),header.begin()+116,' ');std::memcpy(header.data(),text,std::strlen(text));
    header[124]=0;header[125]=1;header[126]='I';header[127]='M';out.write(header.data(),header.size());
    tag(out,miMATRIX,payload);matrixHeader(out,mxSTRUCT,variableName,{1,1});
    tag(out,miINT32,4);put32(out,size32(width));zeros(out,4);
    tag(out,miINT8,mul(width,fields.size()));
    for(const auto& field:fields) {out.write(field.name.data(),field.name.size());zeros(out,width-field.name.size());}
    zeros(out,padded(mul(width,fields.size()))-mul(width,fields.size()));
    for(std::size_t i=0;i<fields.size();++i) {
        const auto& field=fields[i];tag(out,miMATRIX,childSizes[i]);
        matrixHeader(out,mxDOUBLE,"",field.dims);bytes(out,miDOUBLE,field.data,mul(field.count,8));
    }
    out.flush();if(!out)fail("failed writing " + path.string());
}
} // namespace mat5
