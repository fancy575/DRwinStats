#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include "model_math.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace drwk {
using Vec = std::vector<double>;
using IVec = std::vector<int>;
struct Error : std::runtime_error { using std::runtime_error::runtime_error; };
SEXP field(SEXP x, const char* name) {
  SEXP nm=Rf_getAttrib(x,R_NamesSymbol);
  for(R_xlen_t j=0;j<XLENGTH(x);++j)
    if(std::string(CHAR(STRING_ELT(nm,j)))==name) return VECTOR_ELT(x,j);
  throw Error(std::string("Missing internal field: ")+name);
}
double number(SEXP x) { return Rf_asReal(x); }
bool boolean(SEXP x) { return Rf_asLogical(x)==TRUE; }
std::string string(SEXP x) { return CHAR(Rf_asChar(x)); }
Vec vec(SEXP x) {
  Vec out(XLENGTH(x));
  for(R_xlen_t i=0;i<XLENGTH(x);++i) out[i]=TYPEOF(x)==INTSXP?INTEGER(x)[i]:REAL(x)[i];
  return out;
}
IVec ivec(SEXP x) { Vec v=vec(x); return IVec(v.begin(),v.end()); }
struct Matrix {
  int n=0,p=0; Vec z;
  explicit Matrix(SEXP x) {
    SEXP dim=Rf_getAttrib(x,R_DimSymbol);
    if(XLENGTH(dim)!=2) throw Error("An internal design matrix is not a matrix.");
    n=INTEGER(dim)[0];p=INTEGER(dim)[1];z=vec(x);
  }
  double at(int i,int j) const { return z[i+n*j]; }
};
struct Margin {
  Vec beta,times,surv;
  int beta_offset=0,surv_offset=0;
};
struct Model {
  std::string family; double theta=1;
  int theta_offset=-1; std::vector<Margin> margins;
};
Margin read_margin(SEXP beta,SEXP times,SEXP surv,int& offset) {
  Margin x;x.beta=vec(beta);x.times=vec(times);x.surv=vec(surv);
  if(x.times.size()!=x.surv.size()) throw Error("Baseline times and survival values have different lengths.");
  x.beta_offset=offset;offset+=x.beta.size();x.surv_offset=offset;offset+=x.surv.size();
  for(std::size_t j=0;j<x.times.size();++j) {
    if(!std::isfinite(x.times[j]) || (j && x.times[j]<=x.times[j-1])) throw Error("Baseline times must be finite and strictly increasing.");
    if(!std::isfinite(x.surv[j]) || x.surv[j]<0 || x.surv[j]>1 || (j&&x.surv[j]>x.surv[j-1]+1e-10)) throw Error("Invalid fitted baseline survival function.");
  }
  return x;
}
Model read_outcome(SEXP x,int d,int& offset) {
  Model z;z.family=string(field(x,"family"));z.theta=number(field(x,"theta"));
  SEXP b=field(x,"beta"),t=field(x,"base_times"),s=field(x,"base_survival");
  if(XLENGTH(b)!=d||XLENGTH(t)!=d||XLENGTH(s)!=d) throw Error("Outcome model dimension does not match the priorities.");
  for(int j=0;j<d;++j) z.margins.push_back(read_margin(VECTOR_ELT(b,j),VECTOR_ELT(t,j),VECTOR_ELT(s,j),offset));
  z.theta_offset=offset++;
  return z;
}
Margin read_censor(SEXP x,int& offset) {
  return read_margin(field(x,"beta"),field(x,"base_times"),field(x,"base_survival"),offset);
}
struct Prediction {
  Vec s,dbase,dlp; const Margin* model=nullptr; const Matrix* X=nullptr; int row=0;
};
Prediction predict(const Margin& model,const Matrix& X,int row) {
  if(static_cast<int>(model.beta.size())!=X.p) throw Error("Model coefficients and design matrix have different dimensions.");
  Prediction out;out.model=&model;out.X=&X;out.row=row;
  double lp=0;for(int j=0;j<X.p;++j)lp+=model.beta[j]*X.at(row,j);
  if(!std::isfinite(lp)||std::fabs(lp)>700) throw Error("A fitted linear predictor is outside the numerical range.");
  const double e=std::exp(lp);
  for(double s:model.surv) {
    if(s<=0){out.s.push_back(0);out.dbase.push_back(e==1?1:0);out.dlp.push_back(0);}
    else {double v=std::exp(e*std::log(s));out.s.push_back(v);out.dbase.push_back(v*e/s);out.dlp.push_back(v*e*std::log(s));}
  }
  return out;
}
double sv(const Prediction& x,int k){return k==-1?1:(k==-2?0:x.s.at(k));}
int before(const Vec& t,double u){return static_cast<int>(std::lower_bound(t.begin(),t.end(),u)-t.begin())-1;}
int through(const Vec& t,double u){return static_cast<int>(std::upper_bound(t.begin(),t.end(),u)-t.begin())-1;}
struct State {Vec y;IVec lo,hi;int index=-1;};
void enumerate_rec(const Model& m,int j,double limit,double tau,State& s,std::vector<State>& out,int max_states) {
  if(j==static_cast<int>(m.margins.size())) {
    if(static_cast<int>(out.size())>=max_states)throw Error("The exact outcome support exceeds control$max_states. Increase the limit if memory permits; no priorities were discarded.");
    out.push_back(s);return;
  }
  const auto& t=m.margins[j].times;
  for(int k=0;k<static_cast<int>(t.size())&&t[k]<limit;++k){s.y[j]=t[k];s.lo[j]=k-1;s.hi[j]=k;enumerate_rec(m,j+1,limit,tau,s,out,max_states);}
  s.y[j]=tau;s.lo[j]=before(t,limit);s.hi[j]=-2;enumerate_rec(m,j+1,limit,tau,s,out,max_states);
}
std::vector<State> enumerate(const Model& m,double tau,int max_states) {
  const int d=m.margins.size();State s;s.y.resize(d);s.lo.resize(d);s.hi.resize(d);std::vector<State> out;
  const Vec& t=m.margins[0].times;
  for(int k=0;k<static_cast<int>(t.size())&&t[k]<tau;++k){s.y[0]=t[k];s.lo[0]=k-1;s.hi[0]=k;enumerate_rec(m,1,t[k],tau,s,out,max_states);}
  s.y[0]=tau;s.lo[0]=before(t,tau);s.hi[0]=-2;enumerate_rec(m,1,tau,tau,s,out,max_states);return out;
}
struct Rect {double value=0;Vec grad;};
Rect probability(const Model& m,const std::vector<Prediction>& p,const IVec& lo,const IVec& hi,bool derivative) {
  Vec a(p.size()),b(p.size());for(std::size_t j=0;j<p.size();++j){a[j]=sv(p[j],lo[j]);b[j]=sv(p[j],hi[j]);}
  auto r=drw::rectangle_survival(a,b,m.family,m.theta,derivative);
  if(!std::isfinite(r.value)||r.value < -1e-10)throw Error("The fitted copula gives an invalid outcome probability.");
  Rect z;z.value=std::max(0.0,r.value);z.grad=std::move(r.gradient);return z;
}
void add_margin_derivative(const Prediction& p,int idx,double value,const Vec& multipliers,std::vector<Vec>& gradient) {
  if(idx<0||value==0)return;
  for(std::size_t r=0;r<multipliers.size();++r)if(multipliers[r]!=0){
    const double w=value*multipliers[r];
    gradient[r][p.model->surv_offset+idx]+=w*p.dbase[idx];
    const double z=w*p.dlp[idx];
    for(int j=0;j<p.X->p;++j)gradient[r][p.model->beta_offset+j]+=z*p.X->at(p.row,j);
  }
}
void add_rect_derivative(const Model& m,const std::vector<Prediction>& p,const IVec& lo,const IVec& hi,const Rect& r,const Vec& mult,std::vector<Vec>& grad) {
  const int d=p.size();if(r.grad.size()!=static_cast<std::size_t>(2*d+1))throw Error("Internal rectangle derivative has the wrong dimension.");
  for(int j=0;j<d;++j){add_margin_derivative(p[j],lo[j],r.grad[j],mult,grad);add_margin_derivative(p[j],hi[j],r.grad[d+j],mult,grad);}
  for(std::size_t k=0;k<mult.size();++k)grad[k][m.theta_offset]+=mult[k]*r.grad[2*d];
}
struct Inverse {double value=1,derivative=0;int index=-1;};
struct Diagnostics {double minK=1,maxWeight=1,massError=0;int truncated=0;};
Inverse inverse(const Prediction& p,int idx,double floor,bool truncate,Diagnostics& diag) {
  double k=sv(p,idx);bool clipped=k<floor;
  diag.minK=std::min(diag.minK,k);
  if(clipped&&!truncate)throw Error("A required censoring survival probability is below control$weight_floor. Choose an earlier time, revise the censoring model, or explicitly request positivity='truncate'.");
  if(clipped){k=floor;++diag.truncated;}
  if(!(k>0)||!std::isfinite(k))throw Error("Invalid fitted censoring survival probability.");
  Inverse z;z.value=1/k;z.index=idx;z.derivative=clipped?0:-1/(k*k);
  diag.maxWeight=std::max(diag.maxWeight,z.value);return z;
}
struct Node {double t=0,w=0,den=1;IVec lo,hi;Inverse left,right;bool own=false;};
struct Group {int first=0,last=0;IVec observed;};
struct History {std::vector<Node> nodes;std::vector<Group> groups;Vec prefix,times;Inverse resolved_weight;bool resolved=false;};
IVec history_observed(const Vec& events,double t){IVec z(events.size(),0);for(std::size_t j=1;j<events.size();++j)z[j]=events[j]<=t;return z;}
void history_bounds(const Model& m,const Vec& events,double t,IVec& lo,IVec& hi) {
  int d=m.margins.size();lo.resize(d);hi.resize(d);
  for(int j=0;j<d;++j){const auto& times=m.margins[j].times;
    if(j>0&&events[j]<=t){auto it=std::lower_bound(times.begin(),times.end(),events[j]);if(it==times.end()||*it!=events[j])throw Error("An observed event is absent from its fitted outcome support.");int k=it-times.begin();lo[j]=k-1;hi[j]=k;}
    else{lo[j]=through(times,t);hi[j]=-2;}
  }
}
History make_history(const Model& model,const std::vector<Prediction>& p,const Prediction& c,const Vec& events,double stop,bool terminal,double tau,double floor,bool truncate,Diagnostics& diag) {
  History h;const double U=std::min(stop,tau);h.resolved=terminal||stop>=tau;
  if(h.resolved)h.resolved_weight=inverse(c,before(c.model->times,U),floor,truncate,diag);
  const bool censored=!terminal&&stop<tau;bool own_seen=false;
  for(int k=0;k<static_cast<int>(c.model->times.size());++k){double t=c.model->times[k];
    bool own=censored&&t==stop;if(t>U||(!own&&t==U))break;
    Node z;z.t=t;z.own=own;z.left=inverse(c,k-1,floor,truncate,diag);
    if(own){z.w=z.left.value;own_seen=true;}else{z.right=inverse(c,k,floor,truncate,diag);z.w=z.left.value-z.right.value;}
    history_bounds(model,events,t,z.lo,z.hi);z.den=probability(model,p,z.lo,z.hi,false).value;
    if(!(z.den>0))throw Error("The working outcome model assigns zero probability to an observed event history.");
    h.nodes.push_back(std::move(z));
  }
  if(censored&&!own_seen){Node z;z.t=stop;z.own=true;z.left=inverse(c,before(c.model->times,stop),floor,truncate,diag);z.w=z.left.value;history_bounds(model,events,stop,z.lo,z.hi);z.den=probability(model,p,z.lo,z.hi,false).value;if(!(z.den>0))throw Error("An observed event history has zero fitted probability.");h.nodes.push_back(std::move(z));}
  h.prefix.push_back(0);
  for(std::size_t j=0;j<h.nodes.size();++j){const auto& z=h.nodes[j];h.times.push_back(z.t);h.prefix.push_back(h.prefix.back()+z.w/z.den);IVec obs=history_observed(events,z.t);if(h.groups.empty()||h.groups.back().observed!=obs){Group g;g.first=j;g.last=j;g.observed=obs;h.groups.push_back(g);}else h.groups.back().last=j;}
  return h;
}
std::vector<std::pair<int,int>> compatible_ranges(const State& s,const Vec& events,const History& h) {
  std::vector<std::pair<int,int>> out;
  for(const auto& g:h.groups){bool good=true;double until=s.y[0];
    for(std::size_t j=1;j<s.y.size();++j){if(g.observed[j]){if(s.y[j]!=events[j]){good=false;break;}}else until=std::min(until,s.y[j]);}
    if(!good)continue;
    int end=std::min(g.last,static_cast<int>(std::lower_bound(h.times.begin(),h.times.end(),until)-h.times.begin())-1);
    if(end>=g.first)out.push_back({g.first,end});
  }
  return out;
}
double coefficient(const std::vector<std::pair<int,int>>& ranges,const History& h){double z=0;for(auto r:ranges)z+=h.prefix[r.second+1]-h.prefix[r.first];return z;}

struct TreeNode {int first=0,last=0,left=-1,right=-1;Vec lo,hi;};
class Comparator {
  const std::vector<Vec>& coords;Vec delta;IVec order;std::vector<TreeNode> nodes;bool lexical=true;
  int build(int first,int last,int depth){int id=nodes.size();TreeNode z;z.first=first;z.last=last;z.lo=Vec(delta.size(),std::numeric_limits<double>::infinity());z.hi=Vec(delta.size(),-std::numeric_limits<double>::infinity());
    for(int k=first;k<last;++k)for(std::size_t j=0;j<delta.size();++j){z.lo[j]=std::min(z.lo[j],coords[order[k]][j]);z.hi[j]=std::max(z.hi[j],coords[order[k]][j]);}
    nodes.push_back(z);if(last-first>24){int axis=depth%delta.size(),mid=(first+last)/2;std::nth_element(order.begin()+first,order.begin()+mid,order.begin()+last,[&](int a,int b){return coords[a][axis]<coords[b][axis];});int a=build(first,mid,depth+1),b=build(mid,last,depth+1);nodes[id].left=a;nodes[id].right=b;}return id;
  }
  double sum_node(int id,const Vec& w,Vec& sums)const{const auto& z=nodes[id];double s=0;if(z.left<0){for(int k=z.first;k<z.last;++k)s+=w[order[k]];}else s=sum_node(z.left,w,sums)+sum_node(z.right,w,sums);sums[id]=s;return s;}
  double query(int id,const Vec& lo,const Vec& hi,const Vec& w,const Vec& sums)const{const auto& z=nodes[id];bool all=true;for(std::size_t j=0;j<delta.size();++j){if(z.hi[j]<lo[j]||z.lo[j]>hi[j])return 0;all=all&&z.lo[j]>=lo[j]&&z.hi[j]<=hi[j];}if(all)return sums[id];if(z.left>=0)return query(z.left,lo,hi,w,sums)+query(z.right,lo,hi,w,sums);double s=0;for(int k=z.first;k<z.last;++k){int i=order[k];bool ok=true;for(std::size_t j=0;j<delta.size();++j)if(coords[i][j]<lo[j]||coords[i][j]>hi[j]){ok=false;break;}if(ok)s+=w[i];}return s;}
 public:
  Comparator(const std::vector<Vec>& x,const Vec& thresholds):coords(x),delta(thresholds){for(double d:delta)lexical=lexical&&d==0;if(!lexical){order.resize(x.size());std::iota(order.begin(),order.end(),0);build(0,x.size(),0);}}
  // p[0][s] = integral h_+(s,y) dw(y); p[1] is h_-(s,y).
  std::vector<Vec> project(const Vec& w)const{int K=w.size();std::vector<Vec> p(2,Vec(K));if(lexical){double total=std::accumulate(w.begin(),w.end(),0.0),running=0;for(int k=0;k<K;++k){p[0][k]=running;p[1][k]=total-running-w[k];running+=w[k];}return p;}
    Vec sums(nodes.size());sum_node(0,w,sums);double inf=std::numeric_limits<double>::infinity();
    for(int k=0;k<K;++k){Vec lo(delta.size(),-inf),hi(delta.size(),inf);for(std::size_t j=0;j<delta.size();++j){hi[j]=std::nextafter(coords[k][j]-delta[j],-inf);p[0][k]+=query(0,lo,hi,w,sums);hi[j]=inf;lo[j]=std::nextafter(coords[k][j]+delta[j],inf);p[1][k]+=query(0,lo,hi,w,sums);lo[j]=coords[k][j]-delta[j];hi[j]=coords[k][j]+delta[j];}}return p;
  }
};
double dot(const Vec& a,const Vec& b){double z=0;for(std::size_t j=0;j<a.size();++j)z+=a[j]*b[j];return z;}
void add(Vec& a,const Vec& b,double w){if(w==0)return;for(std::size_t j=0;j<a.size();++j)a[j]+=w*b[j];}
struct Target {
  Vec scale,S0,S1,T0,T1,A0,A1;std::vector<Vec> PS0,PS1,PT0,PT1,PA0,PA1;double self[2]={0,0};
  explicit Target(int m,int K):scale(m),S0(K),S1(K),T0(K),T1(K),A0(K),A1(K){}
};
SEXP named_list(const std::vector<std::string>& names){SEXP x=PROTECT(Rf_allocVector(VECSXP,names.size()));SEXP n=PROTECT(Rf_allocVector(STRSXP,names.size()));for(std::size_t j=0;j<names.size();++j)SET_STRING_ELT(n,j,Rf_mkChar(names[j].c_str()));Rf_setAttrib(x,R_NamesSymbol,n);UNPROTECT(2);return x;}

SEXP estimate(SEXP data,SEXP models,SEXP target_arg,SEXP control) {
  int n=Rf_asInteger(field(data,"n")),m=Rf_asInteger(field(data,"m"));
  Matrix events(field(data,"event_times")),Xc(field(data,"Xc"));int d=events.p;
  IVec cl=ivec(field(data,"cluster")),arm=ivec(field(data,"arm")),terminal=ivec(field(data,"terminal")),ca=ivec(field(data,"cluster_arm"));
  for(int& c:cl)--c;
  Vec stop=vec(field(data,"stop")),times=vec(field(data,"times")),delta=vec(field(data,"thresholds")),pi=vec(field(data,"cluster_pi"));
  IVec target=ivec(target_arg);int V=target.size(),E=2*V;
  if(m<2||n<2||d<1||delta.size()!=static_cast<std::size_t>(d))throw Error("Invalid prepared analysis dimensions.");
  std::vector<Matrix>X;SEXP rx=field(data,"X");for(int j=0;j<d;++j)X.emplace_back(VECTOR_ELT(rx,j));
  int L=0;SEXP ro=field(models,"outcome"),rc=field(models,"censoring");
  std::vector<Model>Q;for(int a=0;a<2;++a)Q.push_back(read_outcome(VECTOR_ELT(ro,a),d,L));
  std::vector<Margin>C;for(int a=0;a<2;++a)C.push_back(read_censor(VECTOR_ELT(rc,a),L));
  int max_states=Rf_asInteger(field(control,"max_states"));double max_bytes=number(field(control,"max_bytes"));
  double floor=number(field(control,"weight_floor"));bool truncate=string(field(control,"positivity"))=="truncate";
  bool finite_sample=boolean(field(control,"finite_sample")),need_derivative=boolean(field(models,"compute_if"));
  bool complete=string(field(field(data,"design"),"allocation"))=="complete";
  double gamma=complete?(m-1.0)/m:1,den=m*(m-1.0);
  Vec N(m),w0(m),w1(m);std::vector<IVec> members(m);
  for(int i=0;i<n;++i){++N[cl[i]];members[cl[i]].push_back(i);}
  for(int i=0;i<m;++i){w1[i]=ca[i]==1?1/pi[i]:0;w0[i]=ca[i]==0?1/(1-pi[i]):0;}
  int outRows=times.size()*V,outCols=2*outRows;
  SEXP ans=PROTECT(named_list({"estimate","projection","derivative","min_censor_survival","max_censor_weight","mass_error","truncated_weights","states"}));
  SEXP est=PROTECT(Rf_allocMatrix(REALSXP,outRows,2)),proj=PROTECT(Rf_allocMatrix(REALSXP,m,outCols)),der=PROTECT(Rf_allocMatrix(REALSXP,outCols,L));
  std::fill(REAL(est),REAL(est)+outRows*2,0);std::fill(REAL(proj),REAL(proj)+m*outCols,0);std::fill(REAL(der),REAL(der)+outCols*L,0);
  SET_VECTOR_ELT(ans,0,est);SET_VECTOR_ELT(ans,1,proj);SET_VECTOR_ELT(ans,2,der);
  for(int k=3;k<8;++k){SEXP z=PROTECT(Rf_allocVector(REALSXP,times.size()));SET_VECTOR_ELT(ans,k,z);UNPROTECT(1);}
  for(std::size_t ti=0;ti<times.size();++ti){R_CheckUserInterrupt();double tau=times[ti];Diagnostics diag;
    std::vector<State> support[2]={enumerate(Q[0],tau,max_states),enumerate(Q[1],tau,max_states)};
    std::vector<Vec> coords;for(int a=0;a<2;++a)for(auto& s:support[a])coords.push_back(s.y);
    std::vector<Vec> observed(n,Vec(d));IVec observed_index(n);
    for(int i=0;i<n;++i){double death=terminal[i]?std::min(stop[i],tau):tau;observed[i][0]=death;for(int j=1;j<d;++j){double h=events.at(i,j);observed[i][j]=h<death?h:tau;}coords.push_back(observed[i]);}
    std::sort(coords.begin(),coords.end());coords.erase(std::unique(coords.begin(),coords.end()),coords.end());int K=coords.size();
    double expected_bytes=8.0*K*(3.0*m+22.0*V+20)+8.0*L*E;
    if(expected_bytes>max_bytes)throw Error("The exact estimator exceeds control$max_bytes. Increase the memory budget if available, or reduce the fitted outcome support; no priorities were discarded.");
    for(int a=0;a<2;++a)for(auto& s:support[a])s.index=std::lower_bound(coords.begin(),coords.end(),s.y)-coords.begin();
    for(int i=0;i<n;++i)observed_index[i]=std::lower_bound(coords.begin(),coords.end(),observed[i])-coords.begin();
    std::vector<Vec>B0(m,Vec(K)),B1(m,Vec(K)),A(m,Vec(K));
    for(int i=0;i<n;++i){if(i%20==0)R_CheckUserInterrupt();Vec ev(d);for(int j=0;j<d;++j)ev[j]=events.at(i,j);int c=cl[i];
      for(int a=0;a<2;++a){std::vector<Prediction> p;for(int j=0;j<d;++j)p.push_back(predict(Q[a].margins[j],X[j],i));
        History hist;if(a==arm[i]){Prediction pc=predict(C[a],Xc,i);hist=make_history(Q[a],p,pc,ev,stop[i],terminal[i],tau,floor,truncate,diag);if(hist.resolved)A[c][observed_index[i]]+=hist.resolved_weight.value;}
        double sum=0,aug=hist.resolved?hist.resolved_weight.value:0;
        for(const auto& s:support[a]){double pr=probability(Q[a],p,s.lo,s.hi,false).value;sum+=pr;(a?B1[c]:B0[c])[s.index]+=pr;if(a==arm[i]){double coeff=coefficient(compatible_ranges(s,ev,hist),hist);A[c][s.index]+=pr*coeff;aug+=pr*coeff;}}
        if(std::fabs(sum-1)>1e-6)throw Error("The fitted joint outcome probabilities do not sum to one.");
        if(a==arm[i]){diag.massError=std::max(diag.massError,std::fabs(aug-1));if(std::fabs(aug-1)>1e-6)throw Error("The censoring augmentation failed its constant-outcome identity.");}
      }
    }
    Comparator cmp(coords,delta);std::vector<Target> targets;for(int v=0;v<V;++v){targets.emplace_back(m,K);auto& t=targets.back();for(int i=0;i<m;++i){t.scale[i]=target[v]==0?1/N[i]:1;double s=t.scale[i];add(t.S0,B0[i],s);add(t.S1,B1[i],s);add(t.T0,B0[i],s*w0[i]);add(t.T1,B1[i],s*w1[i]);add(t.A0,A[i],s*w0[i]);add(t.A1,A[i],s*w1[i]);}t.PS0=cmp.project(t.S0);t.PS1=cmp.project(t.S1);t.PT0=cmp.project(t.T0);t.PT1=cmp.project(t.T1);t.PA0=cmp.project(t.A0);t.PA1=cmp.project(t.A1);}
    std::vector<Vec> gradient(E,Vec(L)),raw(m,Vec(E));
    for(int ci=0;ci<m;++ci){R_CheckUserInterrupt();auto P0=cmp.project(B0[ci]),P1=cmp.project(B1[ci]);double self[2]={dot(B1[ci],P0[0]),dot(B1[ci],P0[1])};
      std::vector<Vec> fB0(E,Vec(K)),fB1(E,Vec(K)),fA(E,Vec(K));
      for(int v=0;v<V;++v){auto& t=targets[v];double s=t.scale[ci];for(int r=0;r<2;++r){int e=2*v+r;int rr=1-r;t.self[r]+=s*s*self[r];
          raw[ci][e]=s*(dot(B1[ci],t.PS0[r])+dot(B0[ci],t.PS1[rr]))-2*s*s*self[r]+gamma*s*(w1[ci]*(dot(A[ci],t.PA0[r])-dot(B1[ci],t.PT0[r]))+w0[ci]*(dot(A[ci],t.PA1[rr])-dot(B0[ci],t.PT1[rr])));
          if(need_derivative)for(int k=0;k<K;++k){fB1[e][k]=(s*t.PS0[r][k]-s*s*P0[r][k]-gamma*s*w1[ci]*t.PT0[r][k])/den;fB0[e][k]=(s*t.PS1[rr][k]-s*s*P1[rr][k]-gamma*s*w0[ci]*t.PT1[rr][k])/den;fA[e][k]=gamma*s*(w1[ci]*t.PA0[r][k]+w0[ci]*t.PA1[rr][k])/den;}
      }}
      if(!need_derivative)continue;
      for(int i:members[ci]){Vec ev(d);for(int j=0;j<d;++j)ev[j]=events.at(i,j);
        for(int a=0;a<2;++a){std::vector<Prediction> p;for(int j=0;j<d;++j)p.push_back(predict(Q[a].margins[j],X[j],i));Prediction pc=predict(C[arm[i]],Xc,i);History h;Diagnostics unused;
          if(a==arm[i])h=make_history(Q[a],p,pc,ev,stop[i],terminal[i],tau,floor,truncate,unused);
          std::vector<Vec> difference(E,Vec(h.nodes.size()+1));
          for(const auto& s:support[a]){auto pr=probability(Q[a],p,s.lo,s.hi,true);std::vector<std::pair<int,int>> ranges;double co=0;if(a==arm[i]){ranges=compatible_ranges(s,ev,h);co=coefficient(ranges,h);}
            Vec mult(E);for(int e=0;e<E;++e){mult[e]=(a?fB1[e][s.index]:fB0[e][s.index])+co*fA[e][s.index];if(a==arm[i]){double z=pr.value*fA[e][s.index];for(auto rr:ranges){difference[e][rr.first]+=z;difference[e][rr.second+1]-=z;}}}
            add_rect_derivative(Q[a],p,s.lo,s.hi,pr,mult,gradient);
          }
          if(a==arm[i]){for(int e=0;e<E;++e)for(std::size_t k=1;k<difference[e].size();++k)difference[e][k]+=difference[e][k-1];
            for(std::size_t k=0;k<h.nodes.size();++k){const auto& node=h.nodes[k];auto pr=probability(Q[a],p,node.lo,node.hi,true);Vec outcome_mult(E),censor_mult(E);for(int e=0;e<E;++e){censor_mult[e]=difference[e][k]/node.den;outcome_mult[e]=-node.w*censor_mult[e]/node.den;}
              add_rect_derivative(Q[a],p,node.lo,node.hi,pr,outcome_mult,gradient);add_margin_derivative(pc,node.left.index,node.left.derivative,censor_mult,gradient);if(!node.own)add_margin_derivative(pc,node.right.index,-node.right.derivative,censor_mult,gradient);
            }
            if(h.resolved){Vec mult(E);for(int e=0;e<E;++e)mult[e]=fA[e][observed_index[i]];add_margin_derivative(pc,h.resolved_weight.index,h.resolved_weight.derivative,mult,gradient);}
          }
        }
      }
    }
    double correction=finite_sample&&m>2?(m-1.0)/(m-2.0):1;
    for(int v=0;v<V;++v){auto& t=targets[v];int row=ti*V+v;for(int r=0;r<2;++r){int e=2*v+r,col=2*row+r;double point=(dot(t.S1,t.PS0[r])-t.self[r]+gamma*(dot(t.A1,t.PA0[r])-dot(t.T1,t.PT0[r])))/den;REAL(est)[row+outRows*r]=point;for(int i=0;i<m;++i)REAL(proj)[i+m*col]=correction*(raw[i][e]/(m-1.0)-2*point);for(int k=0;k<L;++k)REAL(der)[col+outCols*k]=gradient[e][k];}}
    REAL(VECTOR_ELT(ans,3))[ti]=diag.minK;REAL(VECTOR_ELT(ans,4))[ti]=diag.maxWeight;REAL(VECTOR_ELT(ans,5))[ti]=diag.massError;REAL(VECTOR_ELT(ans,6))[ti]=diag.truncated;REAL(VECTOR_ELT(ans,7))[ti]=K;
  }
  UNPROTECT(4);return ans;
}
} // namespace drwk

extern "C" SEXP drw_estimate(SEXP data,SEXP models,SEXP targets,SEXP control){
  try{return drwk::estimate(data,models,targets,control);}catch(const std::exception& e){Rf_error("DRwinStats: %s",e.what());}return R_NilValue;
}
