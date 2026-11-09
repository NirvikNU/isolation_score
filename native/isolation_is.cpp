// Full-population IS. Stream noise reference tiles; retain every feature.
// Double-precision cuBLAS evaluates full 144-D dot products. CUDA kernels
// reduce distance blocks immediately; no full pairwise matrix is stored.
#include <cuda.h>
#include <nvrtc.h>
#include <cublas_v2.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
void ck(CUresult x){if(x!=CUDA_SUCCESS){const char* s;cuGetErrorString(x,&s);throw std::runtime_error(s);}}
void cb(cublasStatus_t x){if(x!=CUBLAS_STATUS_SUCCESS)throw std::runtime_error("cuBLAS error "+std::to_string(x));}
void cr(nvrtcResult x){if(x!=NVRTC_SUCCESS)throw std::runtime_error(nvrtcGetErrorString(x));}
std::vector<double> read(const std::string& path){
    std::ifstream f(path,std::ios::binary|std::ios::ate);if(!f)throw std::runtime_error("Cannot read "+path);
    auto n=f.tellg();if(n<0 || uint64_t(n)%8)throw std::runtime_error("Invalid binary size");
    std::vector<double> v(size_t(n)/8);f.seekg(0);f.read((char*)v.data(),n);if(!f)throw std::runtime_error("Read failed");return v;
}
template<class T>CUdeviceptr upload(const std::vector<T>& x){CUdeviceptr p;ck(cuMemAlloc(&p,x.size()*sizeof(T)));ck(cuMemcpyHtoD(p,x.data(),x.size()*sizeof(T)));return p;}
CUdeviceptr alloc(size_t n){CUdeviceptr p;ck(cuMemAlloc(&p,n));return p;}
void launch(CUfunction f,unsigned blocks,void** args){ck(cuLaunchKernel(f,blocks,1,1,256,1,1,0,nullptr,args,nullptr));}
struct Matrix{CUdeviceptr x,norm;int n;CUdeviceptr xf=0;};
struct Timing{double seconds=0,gemm=0,reduction=0;};
struct Engine{
    CUcontext context;CUmodule module;CUfunction fnNorm,fnScale,fnAffinity,fnAffinityFloat;bool mixed;
    cublasHandle_t handle;CUdeviceptr buffer;size_t capacity=0;std::string name;
    Engine(const std::filesystem::path& dir,bool useMixed):mixed(useMixed){
        ck(cuInit(0));CUdevice dev;ck(cuDeviceGet(&dev,0));ck(cuCtxCreate(&context,0,dev));
        char device[256];ck(cuDeviceGetName(device,256,dev));name=device;
        std::ifstream in(dir/"isolation_kernels.cu");if(!in)throw std::runtime_error("Missing kernel source");
        std::string source((std::istreambuf_iterator<char>(in)),std::istreambuf_iterator<char>());
        nvrtcProgram program;cr(nvrtcCreateProgram(&program,source.c_str(),"is_kernels.cu",0,nullptr,nullptr));
        int major,minor;ck(cuDeviceGetAttribute(&major,CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR,dev));
        ck(cuDeviceGetAttribute(&minor,CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR,dev));
        if(major<6)throw std::runtime_error("CUDA compute capability 6.0 or newer required");
        std::string arch="--gpu-architecture=compute_"+std::to_string(std::min(major,8))+(major>=8?"0":std::to_string(minor));
        const char* opts[]={arch.c_str(),"--std=c++14","--fmad=true"};
        nvrtcResult code=nvrtcCompileProgram(program,3,opts);size_t len;
        cr(nvrtcGetProgramLogSize(program,&len));std::string log(len,'\0');cr(nvrtcGetProgramLog(program,log.data()));
        if(code!=NVRTC_SUCCESS)throw std::runtime_error(log);
        cr(nvrtcGetPTXSize(program,&len));std::string ptx(len,'\0');cr(nvrtcGetPTX(program,ptx.data()));cr(nvrtcDestroyProgram(&program));
        ck(cuModuleLoadData(&module,ptx.c_str()));
        ck(cuModuleGetFunction(&fnNorm,module,"norms"));
        ck(cuModuleGetFunction(&fnScale,module,"scale_reduce"));ck(cuModuleGetFunction(&fnAffinity,module,"affinity_reduce"));
        ck(cuModuleGetFunction(&fnAffinityFloat,module,"affinity_reduce_float"));
        cb(cublasCreate(&handle));cb(cublasSetMathMode(handle,CUBLAS_PEDANTIC_MATH));
    }
    ~Engine(){cublasDestroy(handle);cuCtxDestroy(context);}
    void reserve(size_t count){if(count<=capacity)return;if(capacity)ck(cuMemFree(buffer));buffer=alloc(count*8);capacity=count;}
    Matrix matrix(const std::vector<double>& x){
        if(x.size()%144 || x.size()/144>INT_MAX)throw std::runtime_error("Invalid feature shape");
        Matrix m{upload(x),alloc(x.size()/144*8),int(x.size()/144)};
        void* args[]={&m.x,&m.n,&m.norm};launch(fnNorm,(m.n+255)/256,args);ck(cuCtxSynchronize());
        if(mixed){std::vector<float> values(x.size());std::transform(x.begin(),x.end(),values.begin(),[](double v){return float(v);});m.xf=upload(values);}
        return m;
    }
    void dot(CUdeviceptr q,CUdeviceptr r,int Q,int R){
        reserve(size_t(Q)*R);const double one=1,zero=0;
        // Event-major input is column-major [144, event]. Output is [R,Q].
        cb(cublasDgemm(handle,CUBLAS_OP_T,CUBLAS_OP_N,R,Q,144,&one,
                      (const double*)r,144,(const double*)q,144,&zero,(double*)buffer,R));
    }
    void dotFloat(CUdeviceptr q,CUdeviceptr r,int Q,int R){
        reserve(size_t(Q)*R);const float one=1,zero=0;
        cb(cublasSgemm(handle,CUBLAS_OP_T,CUBLAS_OP_N,R,Q,144,&one,
                      (const float*)r,144,(const float*)q,144,&zero,(float*)buffer,R));
    }
    double scale(const Matrix& ap,int block,Timing& timing){
        auto start=Clock::now();double lastReport=-1;CUdeviceptr sums=alloc(size_t(ap.n)*8);ck(cuMemsetD8(sums,0,size_t(ap.n)*8));
        for(int qs=0;qs<ap.n;qs+=block){
            int Q=std::min(block,ap.n-qs);CUdeviceptr q=ap.x+size_t(qs)*144*8,nq=ap.norm+size_t(qs)*8;
            for(int rs=qs;rs<ap.n;rs+=block){
                int R=std::min(block,ap.n-rs);CUdeviceptr r=ap.x+size_t(rs)*144*8,nr=ap.norm+size_t(rs)*8;
                auto t=Clock::now();dot(q,r,Q,R);ck(cuCtxSynchronize());timing.gemm+=elapsed(t);
                double factor=qs==rs?1.0:2.0;t=Clock::now();
                void* args[]={&buffer,&q,&r,&nq,&nr,&Q,&R,&qs,&rs,&factor,&sums};
                launch(fnScale,Q,args);ck(cuCtxSynchronize());timing.reduction+=elapsed(t);
            }
            double done=1-std::pow(double(ap.n-qs-Q)/ap.n,2),now=elapsed(start);
            if(lastReport<0 || now-lastReport>=10 || qs+Q==ap.n){
                std::cout<<std::fixed<<std::setprecision(1)<<"  IS scale "<<done*100<<"%; stage remaining ~"<<now*(1-done)/done/60<<" min"<<std::defaultfloat<<std::setprecision(8)<<std::endl;lastReport=now;
            }
        }
        std::vector<double> host(ap.n);ck(cuMemcpyDtoH(host.data(),sums,host.size()*8));ck(cuMemFree(sums));
        long double total=0;for(double v:host)total+=v;
        timing.seconds=elapsed(start);return double(0.1L*total/(static_cast<long double>(ap.n)*ap.n));
    }
    void release(Matrix& m){for(auto p:{m.x,m.norm,m.xf})if(p)ck(cuMemFree(p));}
};
// Reference chunks are uploaded once, then compared with every AP. Only AP
// features and per-AP running sums remain resident. No noise-noise work.
int main(int argc,char** argv){try{
    if(argc!=5)throw std::runtime_error("Usage: isolation_is CACHE OUTPUT mixed|double scale_factor");
    std::string data=argv[1],out=argv[2];bool mixed=std::string(argv[3])=="mixed";
    double factor=std::stod(argv[4]);if(!(factor>0))throw std::runtime_error("Invalid scale factor");
    auto start=Clock::now();std::cout<<std::setprecision(8);
    auto ax=read(data+"/ap.bin");
    Engine e(std::filesystem::absolute(argv[0]).parent_path(),mixed);
    size_t freeBytes,totalBytes;ck(cuMemGetInfo(&freeBytes,&totalBytes));
    // Leave 512 MiB free beyond the AP cache and bounded GEMM/reference tiles.
    if(ax.size()*(mixed?12:8)+size_t(768)*1024*1024>freeBytes)
        throw std::runtime_error("Insufficient GPU memory for AP features; use MATLAB for small populations or a larger GPU.");
    Matrix ap=e.matrix(ax);std::vector<double>().swap(ax);
    std::ifstream nf(data+"/noise.bin",std::ios::binary|std::ios::ate);
    if(!nf || uint64_t(nf.tellg())%(144*8))throw std::runtime_error("Invalid noise file");
    uint64_t nn=uint64_t(nf.tellg())/(144*8);nf.seekg(0);
    if(ap.n<2 || nn<1)throw std::runtime_error("At least two APs and one noise event required");
    std::cout<<"  IS on "<<e.name<<": "<<ap.n<<" AP / "<<nn<<" noise, all 144 dimensions"<<std::endl;
    Timing t;double h=e.scale(ap,4096,t)*factor/0.1;
    if(!(h>0 && std::isfinite(h)))throw std::runtime_error("Zero or invalid AP spread; IS undefined");
    CUdeviceptr minima=upload(std::vector<double>(ap.n,std::numeric_limits<double>::infinity()));
    CUdeviceptr sa=upload(std::vector<double>(ap.n,0)),sb=upload(std::vector<double>(ap.n,0));
    std::vector<int> ids(ap.n);std::iota(ids.begin(),ids.end(),0);CUdeviceptr qi=upload(ids);
    const int refs=32768,queries=1024;auto scoring=Clock::now();double lastReport=-1;
    for(int isAP=0;isAP<2;++isAP){
        uint64_t count=isAP?uint64_t(ap.n):nn;
        for(uint64_t first=0;first<count;first+=refs){
            int R=int(std::min(uint64_t(refs),count-first));Matrix chunk{};
            if(!isAP){
                std::vector<double> host(size_t(R)*144);nf.read((char*)host.data(),host.size()*8);
                if(!nf)throw std::runtime_error("Noise read failed");chunk=e.matrix(host);
            }else chunk={ap.x+first*144*8,ap.norm+first*8,R,mixed?ap.xf+first*144*4:0};
            int rstart=isAP?int(first):0;
            for(int qs=0;qs<ap.n;qs+=queries){
                int Q=std::min(queries,ap.n-qs);
                CUdeviceptr q=ap.x+size_t(qs)*144*8,nq=ap.norm+size_t(qs)*8,id=qi+size_t(qs)*4;
                CUdeviceptr m=minima+size_t(qs)*8,a=sa+size_t(qs)*8,b=sb+size_t(qs)*8;
                if(mixed)e.dotFloat(ap.xf+size_t(qs)*144*4,chunk.xf,Q,R);else e.dot(q,chunk.x,Q,R);
                void* args[]={&e.buffer,&q,&chunk.x,&nq,&chunk.norm,&id,&Q,&R,&rstart,&isAP,&h,&m,&a,&b};
                launch(mixed?e.fnAffinityFloat:e.fnAffinity,Q,args);ck(cuCtxSynchronize());
            }
            if(!isAP)e.release(chunk);
            double done=double(isAP?nn:0)+first+R,now=elapsed(scoring);
            if(lastReport<0 || now-lastReport>=10 || done==double(nn)+ap.n){
                std::cout<<std::fixed<<std::setprecision(1)<<"  IS affinities "<<100*done/(double(nn)+ap.n)<<"%; stage remaining ~"
                    <<now*((double(nn)+ap.n)/done-1)/60<<" min"<<std::defaultfloat<<std::setprecision(8)<<std::endl;lastReport=now;
            }
        }
    }
    std::vector<double> a(ap.n),b(ap.n),p(ap.n);
    ck(cuMemcpyDtoH(a.data(),sa,ap.n*8));ck(cuMemcpyDtoH(b.data(),sb,ap.n*8));
    for(int i=0;i<ap.n;++i){p[i]=a[i]/(a[i]+b[i]);if(!std::isfinite(p[i]))throw std::runtime_error("Nonfinite affinity");}
    std::ofstream pf(out+"/affinities.bin",std::ios::binary);pf.write((char*)p.data(),p.size()*8);
    std::ofstream jf(out+"/is_result.json");
    jf<<std::setprecision(17)<<"{\"IS\":"<<std::accumulate(p.begin(),p.end(),0.0)/p.size()
      <<",\"CharacteristicDistance\":"<<h<<",\"APCount\":"<<ap.n<<",\"NoiseCount\":"<<nn
      <<",\"Seconds\":"<<elapsed(start)<<",\"ScaleSeconds\":"<<t.seconds<<"}";
    if(!pf || !jf)throw std::runtime_error("Cannot write IS results");
    std::cout<<"  IS complete in "<<elapsed(start)/60<<" min"<<std::defaultfloat<<std::setprecision(8)<<std::endl;
}catch(const std::exception& ex){std::cerr<<"ERROR "<<ex.what()<<std::endl;return 1;}return 0;}
