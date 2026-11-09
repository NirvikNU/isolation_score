// Full-population FN/FP classifier; K=31 with AP-first distance ties.
// Full 144-coordinate Euclidean distances determine labels. The first P
// orthonormal PCA coordinates are used ONLY as conservative search bounds.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>
#include <immintrin.h>
using Clock=std::chrono::steady_clock;
double seconds(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
#ifndef BOUND_DIMS
#define BOUND_DIMS 16
#endif
constexpr int D=144,P=BOUND_DIMS,H=16,T=16;
constexpr double INF=std::numeric_limits<double>::infinity();
struct Stats{
    uint64_t nodes=0, distances=0, initialCertified=0, refined=0, screened=0;
    Stats& operator+=(const Stats& b){screened+=b.screened;nodes+=b.nodes;distances+=b.distances;initialCertified+=b.initialCertified;refined+=b.refined;return *this;}
};
struct Node{
    uint32_t begin=0,end=0,left=UINT32_MAX,right=UINT32_MAX;
    std::array<double,T> lo{},hi{};
    double maxnorm=0;
};
struct Data{
    std::vector<double> x,z,norm; std::vector<float> zf;
    std::vector<uint32_t> order,leaf,originalToRank,originalID;
    std::vector<Node> nodes;
    size_t n=0;
    int bucket=64;
    const double* point(uint32_t i)const{return x.data()+size_t(i)*D;}
    const double* projected(uint32_t i)const{return z.data()+size_t(i)*P;}
    void prepare(const std::vector<double>& q,int threads){
        n=x.size()/D;
        if(x.size()%D || n<1 || n>=UINT32_MAX)throw std::runtime_error("Invalid waveform population size.");
        z.resize(n*P);zf.resize(n*T);norm.resize(n);leaf.resize(n);order.resize(n);
        std::iota(order.begin(),order.end(),0);
        #pragma omp parallel for num_threads(threads) schedule(static)
        for(int64_t ii=0;ii<(int64_t)n;++ii){
            const double* a=point((uint32_t)ii);double* p=z.data()+size_t(ii)*P;
            double sq=0;
            for(int d=0;d<D;++d){sq+=a[d]*a[d];for(int j=0;j<P;++j)p[j]+=a[d]*q[d*P+j];}
            norm[ii]=sq;for(int j=0;j<T;++j)zf[size_t(ii)*T+j]=float(p[j]);
        }
        nodes.reserve(2*n/size_t(bucket/2)+8); build(0,(uint32_t)n);
#ifdef CONTIGUOUS_LEAVES
        // Adjacent tree leaves become adjacent in memory. Geometry, event
        // identity and distances are unchanged; only storage order changes.
        std::vector<double> xx(n*D),zz(n*P),norms(n);std::vector<float> ff(n*T);
        std::vector<uint32_t> leaves(n);originalToRank.resize(n);originalID=order;
        #pragma omp parallel for num_threads(threads) schedule(static)
        for(int64_t k=0;k<(int64_t)n;++k){
            uint32_t old=order[k];originalToRank[old]=(uint32_t)k;
            std::copy_n(point(old),D,xx.data()+size_t(k)*D);
            std::copy_n(projected(old),P,zz.data()+size_t(k)*P);
            std::copy_n(zf.data()+size_t(old)*T,T,ff.data()+size_t(k)*T);norms[k]=norm[old];leaves[k]=leaf[old];
        }
        x.swap(xx);z.swap(zz);zf.swap(ff);norm.swap(norms);leaf.swap(leaves);
        std::iota(order.begin(),order.end(),0);
#else
        originalToRank.resize(n);originalID.resize(n);
        std::iota(originalToRank.begin(),originalToRank.end(),0);originalID=originalToRank;
#endif
    }
    uint32_t build(uint32_t begin,uint32_t end){
        uint32_t id=(uint32_t)nodes.size();nodes.emplace_back();
        nodes[id].begin=begin;nodes[id].end=end;
        if(end-begin<=(uint32_t)bucket){
            nodes[id].lo.fill(INF);nodes[id].hi.fill(-INF);
            for(uint32_t k=begin;k<end;++k){
                uint32_t i=order[k];leaf[i]=id;nodes[id].maxnorm=std::max(nodes[id].maxnorm,norm[i]);
                for(int d=0;d<T;++d){nodes[id].lo[d]=std::min(nodes[id].lo[d],projected(i)[d]);nodes[id].hi[d]=std::max(nodes[id].hi[d],projected(i)[d]);}
            }
            return id;
        }
        // Split selection affects speed only, never the validity of bounds.
        std::array<double,T> lo,hi;lo.fill(INF);hi.fill(-INF);
        uint32_t count=std::min(uint32_t(128),end-begin);
        for(uint32_t j=0;j<count;++j){
            uint32_t i=order[begin+uint64_t(j)*(end-begin)/count];
            for(int d=0;d<T;++d){lo[d]=std::min(lo[d],projected(i)[d]);hi[d]=std::max(hi[d],projected(i)[d]);}
        }
        int dim=0;for(int d=1;d<T;++d)if(hi[d]-lo[d]>hi[dim]-lo[dim])dim=d;
        uint32_t mid=begin+(end-begin)/2;
        std::nth_element(order.begin()+begin,order.begin()+mid,order.begin()+end,[&](uint32_t a,uint32_t b){
            double aa=projected(a)[dim],bb=projected(b)[dim];return aa==bb?a<b:aa<bb;
        });
        uint32_t l=build(begin,mid),r=build(mid,end);
        nodes[id].left=l;nodes[id].right=r;
        nodes[id].maxnorm=std::max(nodes[l].maxnorm,nodes[r].maxnorm);
        for(int d=0;d<T;++d){nodes[id].lo[d]=std::min(nodes[l].lo[d],nodes[r].lo[d]);nodes[id].hi[d]=std::max(nodes[l].hi[d],nodes[r].hi[d]);}
        return id;
    }
};
double distance(const double* a,const double* b,double limit,Stats& s){
    ++s.distances;
    double sum=0;
    for(int j=0;j<D;j+=8){
        __m256d x=_mm256_sub_pd(_mm256_loadu_pd(a+j),_mm256_loadu_pd(b+j));
        __m256d y=_mm256_sub_pd(_mm256_loadu_pd(a+j+4),_mm256_loadu_pd(b+j+4));
        __m256d v=_mm256_add_pd(_mm256_mul_pd(x,x),_mm256_mul_pd(y,y));
        double terms[4];_mm256_storeu_pd(terms,v);
        sum+=(terms[0]+terms[1])+(terms[2]+terms[3]);
        if(sum>limit)return INF; // Nonnegative remaining terms cannot restore eligibility.
    }
    return sum;
}
// A projected Euclidean distance is a lower bound on the full metric.
// Reject only with a conservative rounding allowance. Retained candidates
// are evaluated by the unchanged original-coordinate distance function.
double screened_distance(const Data& ref,uint32_t i,const double* q,const double* z,const float* f,
                         double qnorm,double limit,Stats& s){
    ++s.screened;
    // Fast FP32 screening ONLY. It cannot select neighbours or assign labels.
    // 8e-6*(norm(q)+norm(r)+1) conservatively covers FP32 casting, subtraction,
    // squaring and summation errors for 16 coordinates (over 100 unit roundoffs).
    // Any inconclusive or nonfinite bound goes through the FP64 path below.
    const float* ff=ref.zf.data()+size_t(i)*T;
    __m256 fast=_mm256_setzero_ps();
    for(int j=0;j<T;j+=8){
        __m256 delta=_mm256_sub_ps(_mm256_loadu_ps(ff+j),_mm256_loadu_ps(f+j));
        fast=_mm256_add_ps(fast,_mm256_mul_ps(delta,delta));
    }
    __m128 half=_mm_add_ps(_mm256_castps256_ps128(fast),_mm256_extractf128_ps(fast,1));
    half=_mm_hadd_ps(half,half);half=_mm_hadd_ps(half,half);
    double rough=_mm_cvtss_f32(half),normsum=qnorm+ref.norm[i]+1.0;
    if(std::isfinite(rough) && rough>limit+8e-6*normsum)return INF;
    const double* p=ref.projected(i);__m256d sum=_mm256_setzero_pd();
    for(int j=0;j<P;j+=4){
        __m256d delta=_mm256_sub_pd(_mm256_loadu_pd(p+j),_mm256_loadu_pd(z+j));
        sum=_mm256_add_pd(sum,_mm256_mul_pd(delta,delta));
    }
    double terms[4];_mm256_storeu_pd(terms,sum);
    if((terms[0]+terms[1])+(terms[2]+terms[3])>limit+1e-10*normsum)return INF;
    return distance(q,ref.point(i),limit,s);
}
double lower_bound(const Node& nd,const double* z){
    __m256d sum=_mm256_setzero_pd(),zero=sum;
    for(int j=0;j<T;j+=4){
        __m256d point=_mm256_loadu_pd(z+j);
        __m256d left=_mm256_sub_pd(_mm256_loadu_pd(nd.lo.data()+j),point);
        __m256d right=_mm256_sub_pd(point,_mm256_loadu_pd(nd.hi.data()+j));
        __m256d gap=_mm256_max_pd(zero,_mm256_max_pd(left,right));
        sum=_mm256_add_pd(sum,_mm256_mul_pd(gap,gap));
    }
    double terms[4];_mm256_storeu_pd(terms,sum);
    return (terms[0]+terms[1])+(terms[2]+terms[3]);
}
bool prunable(const Node& nd,const double* z,double qnorm,double radius){
    // Guard floating-point projection/bound errors conservatively. Full
    // distances and boundary ties are always evaluated on the original data.
    return lower_bound(nd,z)>radius+1e-10*(qnorm+nd.maxnorm+1.0);
}
struct Heap{
    std::array<double,H> data{};int n=0;double limit;
    explicit Heap(double bound):limit(bound){}
    double radius()const{return n==H?data[0]:limit;}
    void add(double d){
        if(d>radius())return;
        if(n<H){data[n++]=d;std::push_heap(data.begin(),data.begin()+n);}
        else if(d<data[0]){std::pop_heap(data.begin(),data.begin()+n);data[n-1]=d;std::push_heap(data.begin(),data.begin()+n);}
    }
};
void nearest(const Data& ref,uint32_t id,const double* q,const double* z,const float* f,double qnorm,Heap& heap,Stats& s,uint32_t excluded=UINT32_MAX,double knownBound=-1){
    ++s.nodes;const Node& nd=ref.nodes[id];
    if((knownBound<0?lower_bound(nd,z):knownBound)>heap.radius()+1e-10*(qnorm+nd.maxnorm+1.0))return;
    if(nd.left==UINT32_MAX){
        for(uint32_t k=nd.begin;k<nd.end;++k){
            if(ref.order[k]==excluded)continue;
            double d=screened_distance(ref,ref.order[k],q,z,f,qnorm,heap.radius(),s);
            if(std::isfinite(d))heap.add(d);
        }
        return;
    }
    uint32_t first=nd.left,second=nd.right;
    double lo=lower_bound(ref.nodes[first],z),hi=lower_bound(ref.nodes[second],z);
    if(hi<lo){std::swap(first,second);std::swap(lo,hi);}
    nearest(ref,first,q,z,f,qnorm,heap,s,excluded,lo);
    nearest(ref,second,q,z,f,qnorm,heap,s,excluded,hi);
}
void count_inside(const Data& own,uint32_t node,uint32_t excluded,const double* q,const double* z,const float* f,double qnorm,
                  double radius,bool inclusive,int& count,Stats& s){
    if(count>=H)return;
    ++s.nodes;const Node& nd=own.nodes[node];
    if(prunable(nd,z,qnorm,radius))return;
    if(nd.left==UINT32_MAX){
        for(uint32_t k=nd.begin;k<nd.end && count<H;++k){
            uint32_t i=own.order[k];if(i==excluded)continue;
            double d=screened_distance(own,i,q,z,f,qnorm,radius,s);
            if(inclusive?d<=radius:d<radius)++count;
        }
        return;
    }
    uint32_t first=nd.left,second=nd.right;
    if(lower_bound(own.nodes[second],z)<lower_bound(own.nodes[first],z))std::swap(first,second);
    count_inside(own,first,excluded,q,z,f,qnorm,radius,inclusive,count,s);
    count_inside(own,second,excluded,q,z,f,qnorm,radius,inclusive,count,s);
}
bool classify(const Data& own,const Data& opposite,uint32_t i,bool ownAP,Stats& s){
    // With fewer than 16 available votes a class cannot win K=31. Handle
    // small noise populations without changing K or the detection threshold.
    if(own.n<=H)return !ownAP;
    if(opposite.n<H)return ownAP;
    const double* q=own.point(i);const double* z=own.projected(i);const float* f=own.zf.data()+size_t(i)*T;
    Heap candidates(INF);const Node& leaf=own.nodes[own.leaf[i]];
    for(uint32_t k=leaf.begin;k<leaf.end;++k){
        uint32_t j=own.order[k];if(j==i)continue;
        candidates.add(distance(q,own.point(j),candidates.radius(),s));
    }
    // Small test leaves can contain fewer than 17 points. Fill with distinct
    // same-class events from other leaves; this only loosens the upper bound.
    for(uint32_t j=0;j<own.n && candidates.n<H;++j)
        if(j!=i && own.leaf[j]!=own.leaf[i])candidates.add(distance(q,own.point(j),INF,s));
    if(candidates.n!=H)throw std::runtime_error("Class lacks 16 other events.");
#ifdef EXACT_OWN_AP
    if(ownAP){Heap exact(candidates.radius());nearest(own,0,q,z,f,own.norm[i],exact,s,i);candidates=exact;}
#endif
    // Any 16 same-class points give an upper bound on the true 16th distance.
    Heap other(candidates.radius());
    nearest(opposite,0,q,z,f,own.norm[i],other,s);
    if(other.n<H){++s.initialCertified;return ownAP;}
#ifdef EXACT_OWN_AP
    if(ownAP){++s.refined;return candidates.radius()<=other.radius();}
#endif
    // The opposing 16th distance is now exact: everything outside the initial
    // upper bound is farther away. Count same-class events inside that radius.
    int count=0;
    count_inside(own,0,i,q,z,f,own.norm[i],other.radius(),ownAP,count,s);
    ++s.refined;
    return count>=H?ownAP:!ownAP;
}
struct Neighbor{double d;uint64_t index;};
bool less_neighbor(const Neighbor& a,const Neighbor& b){return a.d==b.d?a.index<b.index:a.d<b.d;}
bool brute(const Data& ap,const Data& noise,bool queryAP,uint32_t i){
    const double* q=queryAP?ap.point(i):noise.point(i);
    std::vector<Neighbor> heap;heap.reserve(31);Stats ignored;
    for(int c=0;c<2;++c){
        const Data& ref=c==0?ap:noise;
        for(uint32_t j=0;j<ref.n;++j){
            if((c==0)==queryAP && j==i)continue;
            double d=distance(q,ref.point(j),INF,ignored);
            Neighbor candidate{d,uint64_t(ref.originalID[j])+(c==0?0:ap.n)};
            if(heap.size()<31){heap.push_back(candidate);std::push_heap(heap.begin(),heap.end(),less_neighbor);}
            else if(less_neighbor(candidate,heap[0])){std::pop_heap(heap.begin(),heap.end(),less_neighbor);heap.back()=candidate;std::push_heap(heap.begin(),heap.end(),less_neighbor);}
        }
    }
    int votes=0;for(auto v:heap)votes+=v.index<ap.n;
    return votes>=16;
}
std::vector<double> read_doubles(const std::string& path){
    std::ifstream f(path,std::ios::binary|std::ios::ate);if(!f)throw std::runtime_error("Cannot read "+path);
    auto bytes=f.tellg();if(bytes<0 || uint64_t(bytes)%8)throw std::runtime_error("Invalid binary size");
    std::vector<double> values(size_t(bytes)/8);f.seekg(0);f.read((char*)values.data(),bytes);if(!f)throw std::runtime_error("Read failed "+path);return values;
}
void selftest(){
    auto started=Clock::now();std::mt19937_64 rng(8831);std::normal_distribution<double> normal;
    std::vector<double> projection(D*P,0);for(int i=0;i<P;++i)projection[i*P+i]=1;
    size_t checked=0;
    for(int trial=0;trial<10;++trial){
        Data a,b;a.bucket=b.bucket=32;
        a.x.resize((40+trial*3)*D);b.x.resize((71+trial*5)*D);
        for(double& v:a.x)v=normal(rng);for(double& v:b.x)v=normal(rng)+(trial%3)*.25;
        if(trial==3 || trial==4){
            for(size_t i=0;i<a.x.size();++i)a.x[i]=std::round(a.x[i]);
            for(size_t i=0;i<b.x.size();++i)b.x[i]=std::round(b.x[i]);
            for(int i=0;i<35;++i)std::copy_n(a.x.data(),D,b.x.data()+i*D);
        }
        if(trial>=5){
            // Projection contains no separating information: the classifier
            // must still use the remaining 128 coordinates correctly.
            for(size_t i=0;i<a.x.size()/D;++i)std::fill_n(a.x.data()+i*D,P,0);
            for(size_t i=0;i<b.x.size()/D;++i)std::fill_n(b.x.data()+i*D,P,0);
        }
        if(trial==8){std::fill(a.x.begin(),a.x.end(),1);std::fill(b.x.begin(),b.x.end(),1);}
        if(trial==9){std::fill(a.x.begin(),a.x.end(),1);std::fill(b.x.begin(),b.x.end(),2);}
        a.prepare(projection,2);b.prepare(projection,2);
        for(int c=0;c<2;++c){
            Data& own=c==0?a:b;Data& other=c==0?b:a;
            for(uint32_t i=0;i<own.n;++i){Stats s;bool got=classify(own,other,i,c==0,s),want=brute(a,b,c==0,i);if(got!=want)throw std::runtime_error("Self-test classification mismatch");++checked;}
        }
    }
    std::cout<<"SELFTEST_PASSED "<<checked<<" labels; "<<seconds(started)<<" seconds\n";
}
struct Run{
    uint64_t events=0,errors=0;Stats stats;double elapsed=0;
};
Run run_queries(const Data& ap,const Data& noise,bool isAP,const std::vector<uint32_t>& ids,int threads,
                std::ofstream* output=nullptr,std::vector<unsigned char>* originalLabels=nullptr){
    const Data& own=isAP?ap:noise;const Data& other=isAP?noise:ap;
    std::vector<Stats> stats(threads);std::vector<uint64_t> errors(threads,0);
    std::vector<unsigned char> labels(ids.size());
    auto start=Clock::now();
    #pragma omp parallel num_threads(threads)
    {
        int t=omp_get_thread_num();
        #pragma omp for schedule(dynamic,16)
        for(int64_t k=0;k<(int64_t)ids.size();++k){
            bool label=classify(own,other,own.originalToRank[ids[k]],isAP,stats[t]);
            labels[k]=(unsigned char)label;if(originalLabels)(*originalLabels)[ids[k]]=(unsigned char)label;errors[t]+=label!=isAP;
        }
    }
    Run r;r.events=ids.size();r.elapsed=seconds(start);
    for(int t=0;t<threads;++t){r.errors+=errors[t];r.stats+=stats[t];}
    if(output)output->write((const char*)labels.data(),labels.size());
    return r;
}
int main(int argc,char** argv){
    try{
        std::cout<<std::setprecision(12);
        if(argc==2 && std::string(argv[1])=="--selftest"){selftest();return 0;}
        if(argc<3)throw std::runtime_error("Usage: isolation_knn CACHE OUTPUT [threads]");
        std::string dir=argv[1],outdir=argv[2];int threads=argc>3?std::stoi(argv[3]):12;
        if(threads<1 || threads>256)throw std::runtime_error("Invalid thread count");
        auto start=Clock::now();Data ap,noise;
        auto projection=read_doubles(dir+"/projection"+std::to_string(P)+".bin");
        if(projection.size()!=D*P)throw std::runtime_error("Projection shape mismatch");
        ap.x=read_doubles(dir+"/ap.bin");noise.x=read_doubles(dir+"/noise.bin");
        if(ap.x.size()/D+noise.x.size()/D<=31)throw std::runtime_error("K=31 requires more than 31 total events");
        std::cout<<"LOADED seconds="<<seconds(start)<<" AP="<<ap.x.size()/D<<" noise="<<noise.x.size()/D<<std::endl;
        auto build=Clock::now();ap.prepare(projection,threads);std::cout<<"AP_TREE_READY seconds="<<seconds(build)<<std::endl;
        noise.prepare(projection,threads);std::cout<<"BOTH_TREES_READY seconds="<<seconds(build)<<std::endl;
        Run total[2];double buildSeconds=seconds(build);
        for(int c=0;c<2;++c){
            const Data& own=c==0?ap:noise;
            std::vector<unsigned char> labels(own.n);
            std::ofstream output(outdir+(c==0?"/ap_labels.bin":"/noise_labels.bin"),std::ios::binary);
            if(!output)throw std::runtime_error("Cannot create label checkpoint");
            auto scoring=Clock::now(); double lastReport=-1;
            for(uint32_t begin=0;begin<own.n;begin+=20000){
                uint32_t end=std::min<uint32_t>((uint32_t)own.n,begin+20000);
                std::vector<uint32_t> ids(end-begin);for(uint32_t j=begin;j<end;++j)ids[j-begin]=own.originalID[j];
                Run r=run_queries(ap,noise,c==0,ids,threads,nullptr,&labels);
                total[c].events+=r.events;total[c].errors+=r.errors;total[c].stats+=r.stats;total[c].elapsed=seconds(scoring);
                double now=seconds(scoring);
                if(lastReport<0 || now-lastReport>=10 || end==own.n){
                    std::cout<<std::fixed<<std::setprecision(1)<<"  FN/FP "<<(c==0?"AP":"noise")<<" "<<end<<"/"<<own.n
                        <<"; stage remaining ~"<<now*(own.n-end)/end/60<<" min"<<std::defaultfloat<<std::setprecision(12)<<std::endl;
                    lastReport=now;
                }
            }
            output.write((const char*)labels.data(),labels.size());
        }
        double fn=double(total[1].errors)/(total[1].errors+ap.n),fp=double(total[0].errors)/ap.n;
        std::ofstream result(outdir+"/full_scores.json");
        result<<std::setprecision(17)<<"{\n\"APCount\":"<<ap.n<<",\"NoiseCount\":"<<noise.n
              <<",\"FNCount\":"<<total[1].errors<<",\"FPCount\":"<<total[0].errors
              <<",\"FN\":"<<fn<<",\"FP\":"<<fp<<",\"Threads\":"<<threads
              <<",\"BuildSeconds\":"<<buildSeconds<<",\"APSeconds\":"<<total[0].elapsed
              <<",\"NoiseSeconds\":"<<total[1].elapsed<<",\"TotalSeconds\":"<<seconds(start)
              <<",\"DistanceEvaluations\":"<<total[0].stats.distances+total[1].stats.distances
              <<",\"InitiallyCertified\":"<<total[0].stats.initialCertified+total[1].stats.initialCertified<<"}\n";
        std::cout<<"FULL_COMPLETE FN="<<fn<<" FP="<<fp<<" seconds="<<seconds(start)<<std::endl;
    }catch(const std::exception& e){std::cerr<<"ERROR "<<e.what()<<std::endl;return 1;}
    return 0;
}
