#ifndef DRWINS_MODEL_MATH_H
#define DRWINS_MODEL_MATH_H
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <string>
#include <limits>
namespace drwins_model {
using real = long double;
constexpr int max_local_derivatives = 65;
struct AD {
  real v; int n; std::array<real,max_local_derivatives> g;
  AD(real value=0,int dim=0):v(value),n(dim) { if(dim<0||dim>max_local_derivatives)throw std::runtime_error("Too many endpoints");for(int j=0;j<n;++j)g[j]=0; }
  static AD variable(real value,int dim,int pos){AD x(value,dim);x.g[pos]=1;return x;}
};
inline AD operator+(const AD&a,const AD&b){AD z(a.v+b.v,a.n);for(int j=0;j<z.n;++j)z.g[j]=a.g[j]+b.g[j];return z;}
inline AD operator-(const AD&a,const AD&b){AD z(a.v-b.v,a.n);for(int j=0;j<z.n;++j)z.g[j]=a.g[j]-b.g[j];return z;}
inline AD operator-(const AD&a){AD z(-a.v,a.n);for(int j=0;j<z.n;++j)z.g[j]=-a.g[j];return z;}
inline AD operator*(const AD&a,const AD&b){AD z(a.v*b.v,a.n);for(int j=0;j<z.n;++j)z.g[j]=a.g[j]*b.v+a.v*b.g[j];return z;}
inline AD operator/(const AD&a,const AD&b){AD z(a.v/b.v,a.n);for(int j=0;j<z.n;++j)z.g[j]=(a.g[j]-z.v*b.g[j])/b.v;return z;}
inline AD operator+(const AD&a,real b){AD z=a;z.v+=b;return z;}inline AD operator+(real a,const AD&b){return b+a;}
inline AD operator-(const AD&a,real b){AD z=a;z.v-=b;return z;}inline AD operator-(real a,const AD&b){return AD(a,b.n)-b;}
inline AD operator*(const AD&a,real b){AD z(a.v*b,a.n);for(int j=0;j<z.n;++j)z.g[j]=a.g[j]*b;return z;}inline AD operator*(real a,const AD&b){return b*a;}
inline AD operator/(const AD&a,real b){return a*(1/b);}inline AD operator/(real a,const AD&b){return AD(a,b.n)/b;}
inline AD exp(const AD&a){real e=std::exp(a.v);AD z(e,a.n);for(int j=0;j<z.n;++j)z.g[j]=e*a.g[j];return z;}
inline AD log(const AD&a){if(!(a.v>0))throw std::runtime_error("Nonpositive probability");AD z(std::log(a.v),a.n);for(int j=0;j<z.n;++j)z.g[j]=a.g[j]/a.v;return z;}
inline AD log1p(const AD&a){if(!(a.v>-1))throw std::runtime_error("Invalid probability");AD z(std::log1p(a.v),a.n);for(int j=0;j<z.n;++j)z.g[j]=a.g[j]/(1+a.v);return z;}
inline AD expm1(const AD&a){AD z(std::expm1(a.v),a.n);real e=std::exp(a.v);for(int j=0;j<z.n;++j)z.g[j]=e*a.g[j];return z;}
inline AD pow(const AD&a,const AD&b){if(a.v==0)return AD(0,a.n);return exp(b*log(a));}
inline AD log_survival(const std::vector<AD>&a,const AD&theta,int family){
  int nd=theta.n;std::vector<int>active;for(std::size_t j=0;j<a.size();++j)if(a[j].v>0)active.push_back(j);
  if(active.empty())return AD(0,nd);
  if(active.size()==1)return -a[active[0]];
  if(family==0){AD s(0,nd);for(int j:active)s=s+a[j];return -s;}
  if(family==1){ // Gumbel
    real high=0;for(int j:active)high=std::max(high,a[j].v);
    AD sum(0,nd);for(int j:active)sum=sum+pow(a[j]/high,theta);
    return -high*exp(log(sum)/theta);
  }
  if(family==2){ // Clayton, evaluated by log-sum-exp.
    real high=0;for(int j:active)high=std::max(high,theta.v*a[j].v);
    AD sum(0,nd);for(int j:active)sum=sum+exp(theta*a[j]-high);
    sum=sum-(static_cast<real>(active.size())-1)*std::exp(-high);
    return -(log(sum)+high)/theta;
  }
  if(family==3){ // Frank, positive dependence in arbitrary dimension.
    AD den=-expm1(-theta),lprod(0,nd);
    for(int j:active)lprod=lprod+log(-expm1(-theta*exp(-a[j])))-log(den);
    AD result=-log1p(-den*exp(lprod))/theta;
    return log(result);
  }
  if(family==4){ // Joe.
    AD lprod(0,nd);
    for(int j:active){AD u=exp(-a[j]);lprod=lprod+log(-expm1(theta*log1p(-u)));}
    AD result=-expm1(log1p(-exp(lprod))/theta);
    return log(result);
  }
  throw std::runtime_error("Unknown copula family");
}
inline real log_survival_value(const std::vector<real>&a,real theta,int family){
  std::vector<AD> x;for(real v:a)x.emplace_back(v,0);return log_survival(x,AD(theta,0),family).v;
}
}
namespace drw {
struct RectResult { double value=0; std::vector<double> gradient; };
inline int copula_family_id(const std::string& family) {
  if(family=="independence"||family=="independent")return 0;
  if(family=="gumbel")return 1;
  if(family=="clayton")return 2;
  if(family=="frank")return 3;
  if(family=="joe")return 4;
  throw std::runtime_error("Unknown copula family: "+family);
}
struct CopulaResult {drwins_model::real value=0;std::vector<drwins_model::real> du;drwins_model::real dt=0;};
inline CopulaResult copula_eval(const std::vector<double>&u,int family,double theta,bool derivatives) {
  using namespace drwins_model;
  const int d=u.size();CopulaResult out;if(derivatives)out.du.assign(d,0);
  if(!d){out.value=1;return out;}
  if((family==1||family==4)&&theta<=1)family=0;
  if((family==2||family==3)&&theta<=0)family=0;
  int zeros=0,zero=-1;for(int j=0;j<d;++j){if(u[j]<0||u[j]>1||!std::isfinite(u[j]))throw std::runtime_error("Invalid marginal survival probability");if(u[j]==0){++zeros;zero=j;}}
  if(family==0){out.value=1;for(double v:u)out.value*=v;if(derivatives)for(int j=0;j<d;++j){real v=1;for(int k=0;k<d;++k)if(k!=j)v*=u[k];out.du[j]=v;}return out;}
  if(zeros){if(derivatives&&zeros==1){real v=1;if(family==3){real den=-std::expm1(-static_cast<real>(theta));for(int k=0;k<d;++k)if(k!=zero)v*=(-std::expm1(-static_cast<real>(theta)*u[k]))/den;}else if(family==4){for(int k=0;k<d;++k)if(k!=zero)v*=-std::expm1(theta*std::log1p(-static_cast<real>(u[k])));}out.du[zero]=v;}return out;}
  std::vector<AD>a;for(double v:u)a.emplace_back(-std::log(static_cast<real>(v)),derivatives?1:0);
  AD th(theta,derivatives?1:0);if(derivatives)th.g[0]=1;
  AD lc=log_survival(a,th,family);out.value=std::exp(lc.v);if(!derivatives)return out;out.dt=out.value*lc.g[0];
  if(d==1){out.du[0]=1;out.dt=0;return out;}
  if(family==1){real R=-lc.v;for(int j=0;j<d;++j)out.du[j]=R>0?out.value/u[j]*std::pow(a[j].v/R,theta-1):1;}
  else if(family==2){for(int j=0;j<d;++j)out.du[j]=std::exp((1+theta)*(lc.v-std::log(static_cast<real>(u[j]))));}
  else if(family==3){real num=std::expm1(theta*out.value);for(int j=0;j<d;++j)out.du[j]=num/std::expm1(theta*static_cast<real>(u[j]));}
  else if(family==4){real C=out.value;if(C<1){real scale=(-std::expm1(theta*std::log1p(-C)))/std::pow(1-C,theta-1);for(int j=0;j<d;++j){real v=u[j];out.du[j]=v==1?0:scale*std::pow(1-v,theta-1)/(-std::expm1(theta*std::log1p(-v)));}}else{std::fill(out.du.begin(),out.du.end(),1);}}
  return out;
}
inline double copula_survival(const std::vector<double>&u,const std::string&family,double theta){return static_cast<double>(copula_eval(u,copula_family_id(family),theta,false).value);}
inline RectResult rectangle_survival(const std::vector<double>&lo,const std::vector<double>&hi,const std::string&family,double theta,bool derivatives=true){
  using drwins_model::real;int d=lo.size();if(hi.size()!=lo.size()||d>20)throw std::runtime_error("Unsupported survival rectangle dimension");
  for(int j=0;j<d;++j)if(lo[j]<hi[j]||lo[j]>1||hi[j]<0)throw std::runtime_error("Invalid survival rectangle bounds");
  int f=copula_family_id(family);RectResult out;if(derivatives)out.gradient.assign(2*d+1,0);std::vector<real>g(derivatives?2*d+1:0,0);real v=0;
  if(f==0 || ((f==1||f==4)&&theta<=1) || ((f==2||f==3)&&theta<=0)){
    v=1;for(int j=0;j<d;++j)v*=static_cast<real>(lo[j])-hi[j];
    if(derivatives)for(int j=0;j<d;++j){real a=1;for(int k=0;k<d;++k)if(k!=j)a*=static_cast<real>(lo[k])-hi[k];g[j]=a;g[d+j]=-a;}
  } else {
    const std::size_t count=std::size_t(1)<<d;
    for(std::size_t mask=0;mask<count;++mask){std::vector<double>u(d);int parity=0;for(int j=0;j<d;++j){bool use=mask&(std::size_t(1)<<j);u[j]=use?hi[j]:lo[j];parity+=use;}real sign=parity%2?-1:1;auto z=copula_eval(u,f,theta,derivatives);v+=sign*z.value;if(derivatives){for(int j=0;j<d;++j)g[(mask&(std::size_t(1)<<j))?d+j:j]+=sign*z.du[j];g[2*d]+=sign*z.dt;}}
  }
  if(v<0&&v>=-1e-12L)v=0;
  if(v<0||v>1+1e-10L||!std::isfinite(v))throw std::runtime_error("Invalid fitted outcome rectangle probability");
  out.value=static_cast<double>(std::min(v,static_cast<real>(1)));if(derivatives)for(int j=0;j<2*d+1;++j)out.gradient[j]=static_cast<double>(g[j]);return out;
}
}
#endif
