// Full-dimensional IS: FP64 scale; FP64 or mixed-precision affinities.
#define D 144
#define THREADS 256
extern "C" __global__ void norms(const double* x,int n,double* output){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    double s=0;for(int d=0;d<D;++d){double v=x[(long long)i*D+d];s+=v*v;}output[i]=s;
}
__device__ double distance(double dot,const double* q,const double* r,double qn,double rn){
    double d2=qn+rn-2.0*dot;
    // Recompute nearly coincident vectors directly to avoid cancellation in
    // norm(q)^2+norm(r)^2-2*dot(q,r). This also handles duplicate waveforms.
    if(d2<=512.0*2.2204460492503131e-16*(qn+rn)){
        d2=0;for(int d=0;d<D;++d){double delta=q[d]-r[d];d2+=delta*delta;}
    }
    return sqrt(fmax(0.0,d2));
}
extern "C" __global__ void scale_reduce(const double* dots,const double* q,const double* r,
    const double* nq,const double* nr,int Q,int R,int qstart,int rstart,double factor,double* sums){
    int qi=blockIdx.x,t=threadIdx.x;if(qi>=Q)return;
    double s=0;
    for(int ri=t;ri<R;ri+=THREADS){
        if(qstart+qi==rstart+ri)continue;
        s+=distance(dots[(long long)qi*R+ri],q+(long long)qi*D,r+(long long)ri*D,nq[qi],nr[ri]);
    }
    __shared__ double partial[THREADS];partial[t]=s;__syncthreads();
    for(int k=THREADS/2;k;k/=2){if(t<k)partial[t]+=partial[t+k];__syncthreads();}
    if(t==0)sums[qstart+qi]+=factor*partial[0];
}
extern "C" __global__ void affinity_reduce(double* dots,const double* q,const double* r,
    const double* nq,const double* nr,const int* ids,int Q,int R,int rstart,int isAP,
    double h,double* minima,double* sumAP,double* sumNoise){
    int qi=blockIdx.x,t=threadIdx.x;if(qi>=Q)return;
    double smallest=__longlong_as_double(0x7ff0000000000000LL);
    for(int ri=t;ri<R;ri+=THREADS){
        double d=(isAP && ids[qi]==rstart+ri)?__longlong_as_double(0x7ff0000000000000LL):
            distance(dots[(long long)qi*R+ri],q+(long long)qi*D,r+(long long)ri*D,nq[qi],nr[ri]);
        dots[(long long)qi*R+ri]=d;smallest=fmin(smallest,d);
    }
    __shared__ double part[THREADS];part[t]=smallest;__syncthreads();
    for(int k=THREADS/2;k;k/=2){if(t<k)part[t]=fmin(part[t],part[t+k]);__syncthreads();}
    double tileMin=part[0];__syncthreads();
    if(!isfinite(tileMin))return; // A tile containing only the excluded self.
    double s=0;for(int ri=t;ri<R;ri+=THREADS)s+=exp(-(dots[(long long)qi*R+ri]-tileMin)/h);
    part[t]=s;__syncthreads();
    for(int k=THREADS/2;k;k/=2){if(t<k)part[t]+=part[t+k];__syncthreads();}
    if(t==0){
        double old=minima[qi],next=fmin(old,tileMin);
        double shrink=exp((next-old)/h),addition=part[0]*exp((next-tileMin)/h);
        sumAP[qi]=sumAP[qi]*shrink+(isAP?addition:0.0);
        sumNoise[qi]=sumNoise[qi]*shrink+(isAP?0.0:addition);minima[qi]=next;
    }
}

__device__ double distance_float(double dot,const double* q,const double* r,double qn,double rn){
    double d2=qn+rn-2.0*dot;
    // Recompute nearly coincident vectors directly to avoid cancellation in
    // norm(q)^2+norm(r)^2-2*dot(q,r). This also handles duplicate waveforms.
    if(!isfinite(d2) || qn+rn<1e-20 || d2<=2e-5*(qn+rn)){
        d2=0;for(int d=0;d<D;++d){double delta=q[d]-r[d];d2+=delta*delta;}
    }
    return sqrt(fmax(0.0,d2));
}

extern "C" __global__ void affinity_reduce_float(float* dots,const double* q,const double* r,
    const double* nq,const double* nr,const int* ids,int Q,int R,int rstart,int isAP,
    double h,double* minima,double* sumAP,double* sumNoise){
    int qi=blockIdx.x,t=threadIdx.x;if(qi>=Q)return;
    double smallest=__longlong_as_double(0x7ff0000000000000LL);
    for(int ri=t;ri<R;ri+=THREADS){
        double d=(isAP && ids[qi]==rstart+ri)?__longlong_as_double(0x7ff0000000000000LL):
            distance_float(dots[(long long)qi*R+ri],q+(long long)qi*D,r+(long long)ri*D,nq[qi],nr[ri])/h;
        dots[(long long)qi*R+ri]=d;smallest=fmin(smallest,d);
    }
    __shared__ double part[THREADS];part[t]=smallest;__syncthreads();
    for(int k=THREADS/2;k;k/=2){if(t<k)part[t]=fmin(part[t],part[t+k]);__syncthreads();}
    double tileMin=part[0];__syncthreads();
    if(!isfinite(tileMin))return; // A tile containing only the excluded self.
    double s=0;for(int ri=t;ri<R;ri+=THREADS)s+=double(expf(float(-(double(dots[(long long)qi*R+ri])-tileMin))));
    part[t]=s;__syncthreads();
    for(int k=THREADS/2;k;k/=2){if(t<k)part[t]+=part[t+k];__syncthreads();}
    if(t==0){
        double old=minima[qi],next=fmin(old,tileMin);
        double shrink=exp(next-old),addition=part[0]*exp(next-tileMin);
        sumAP[qi]=sumAP[qi]*shrink+(isAP?addition:0.0);
        sumNoise[qi]=sumNoise[qi]*shrink+(isAP?0.0:addition);minima[qi]=next;
    }
}
