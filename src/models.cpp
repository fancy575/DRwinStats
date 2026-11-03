#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include "model_math.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <limits>
namespace {
SEXP field(SEXP x,const char* name){SEXP nm=Rf_getAttrib(x,R_NamesSymbol);for(R_xlen_t i=0;i<XLENGTH(x);++i)if(std::string(CHAR(STRING_ELT(nm,i)))==name)return VECTOR_ELT(x,i);Rf_error("Missing native model field: %s",name);return R_NilValue;}
SEXP optional_field(SEXP x,const char* name){SEXP nm=Rf_getAttrib(x,R_NamesSymbol);for(R_xlen_t i=0;i<XLENGTH(x);++i)if(std::string(CHAR(STRING_ELT(nm,i)))==name)return VECTOR_ELT(x,i);return R_NilValue;}
void names(SEXP x,const std::vector<const char*>&nm){SEXP n=PROTECT(Rf_allocVector(STRSXP,nm.size()));for(std::size_t i=0;i<nm.size();++i)SET_STRING_ELT(n,i,Rf_mkChar(nm[i]));Rf_setAttrib(x,R_NamesSymbol,n);UNPROTECT(1);}
int rows(SEXP x){return INTEGER(Rf_getAttrib(x,R_DimSymbol))[0];}int cols(SEXP x){return INTEGER(Rf_getAttrib(x,R_DimSymbol))[1];}
}
extern "C" SEXP drw_outcome_likelihood(SEXP par,SEXP problem,SEXP score_arg){
  using namespace drwins_model;
  SEXP xl=field(problem,"X"),lower=field(problem,"lower"),upper=field(problem,"upper"),counts=field(problem,"counts"),clusters=field(problem,"cluster");
  int d=XLENGTH(xl),n=rows(lower),m=Rf_asInteger(field(problem,"m")),family=Rf_asInteger(field(problem,"family"));
  bool association=Rf_asLogical(field(problem,"association")),scores=Rf_asLogical(score_arg);
  if(d>30)Rf_error("At most 30 outcome components are supported by the likelihood.");
  SEXP ridge_weights=optional_field(problem,"ridge_weights"),baseline_center=optional_field(problem,"baseline_center"),case_weights=optional_field(problem,"weights");
  double ridge=Rf_asReal(field(problem,"ridge")),aridge=Rf_asReal(field(problem,"association_ridge")),bpen=Rf_asReal(field(problem,"baseline_penalty"));
  std::vector<int>p(d),bo(d),ho(d);int q=0;
  for(int j=0;j<d;++j){p[j]=cols(VECTOR_ELT(xl,j));bo[j]=q;q+=p[j];}
  int ti=association?q++:-1;for(int j=0;j<d;++j){ho[j]=q;q+=INTEGER(counts)[j];}
  if(XLENGTH(par)!=q)Rf_error("Outcome parameter dimension mismatch.");
  SEXP ans=PROTECT(Rf_allocVector(VECSXP,4));SEXP value=PROTECT(Rf_allocVector(REALSXP,1)),grad=PROTECT(Rf_allocVector(REALSXP,q)),score=PROTECT(Rf_allocMatrix(REALSXP,scores?m:0,scores?q:0));
  std::fill(REAL(grad),REAL(grad)+q,0);if(scores)std::fill(REAL(score),REAL(score)+static_cast<std::size_t>(m)*q,0);
  std::vector<std::vector<real>>haz(d),cum(d),ch(d);std::vector<std::vector<std::vector<real>>>cch;
  if(scores)cch.resize(d);
  real theta=association?((family==1||family==4)?1:0)+std::exp(static_cast<real>(REAL(par)[ti])):Rf_asReal(field(problem,"theta"));
  int adim=2*d+1;real objective=0;
  bool valid=true;
  try {
    for(int j=0;j<d;++j){int k=INTEGER(counts)[j];haz[j].resize(k);cum[j].assign(k+1,0);ch[j].assign(k+1,0);if(scores)cch[j].assign(m,std::vector<real>(k+1,0));for(int l=0;l<k;++l){haz[j][l]=std::exp(static_cast<real>(REAL(par)[ho[j]+l]));cum[j][l+1]=cum[j][l]+haz[j][l];}}
    for(int i=0;i<n;++i){
      if(i%128==0)R_CheckUserInterrupt();
      std::vector<real>mult(d);std::vector<AD>lo(d),hi(d);std::vector<int>observed;
      for(int j=0;j<d;++j){SEXP x=VECTOR_ELT(xl,j);real lp=0;for(int k=0;k<p[j];++k)lp+=REAL(x)[i+n*k]*REAL(par)[bo[j]+k];mult[j]=std::exp(lp);int il=INTEGER(lower)[i+n*j],ih=INTEGER(upper)[i+n*j];real av0=cum[j][il]*mult[j],av1=cum[j][ih]*mult[j];lo[j]=AD::variable(av0,adim,2*j);hi[j]=AD::variable(av1,adim,2*j+1);if(ih>il)observed.push_back(j);}
      AD th=AD::variable(theta,adim,2*d),ll(0,adim);
      if(family==0||d==1){for(int j=0;j<d;++j){ll=ll-lo[j];if(INTEGER(upper)[i+n*j]>INTEGER(lower)[i+n*j])ll=ll+log(-expm1(-(hi[j]-lo[j])));}}
      else {
        AD ref=log_survival(lo,th,family),prob(0,adim);std::size_t combinations=std::size_t(1)<<observed.size();
        if(observed.size()>18)throw std::runtime_error("More than 18 observed components require too many likelihood rectangles.");
        for(std::size_t mask=0;mask<combinations;++mask){std::vector<AD>a=lo;int parity=0;for(std::size_t k=0;k<observed.size();++k)if(mask&(std::size_t(1)<<k)){a[observed[k]]=hi[observed[k]];++parity;}AD term=exp(log_survival(a,th,family)-ref);prob=(parity%2)?prob-term:prob+term;}
        ll=ref+log(prob);
      }
      if(case_weights!=R_NilValue)ll=ll*static_cast<real>(REAL(case_weights)[i]);
      if (!std::isfinite(ll.v)) {
        throw std::runtime_error("Nonfinite outcome likelihood");
      }
      objective-=ll.v;int ci=INTEGER(clusters)[i]-1;
      for(int j=0;j<d;++j){SEXP x=VECTOR_ELT(xl,j);real contribution=ll.g[2*j]*lo[j].v+ll.g[2*j+1]*hi[j].v;for(int k=0;k<p[j];++k){double z=static_cast<double>(contribution*REAL(x)[i+n*k]);REAL(grad)[bo[j]+k]-=z;if(scores)REAL(score)[ci+m*(bo[j]+k)]+=z;}int il=INTEGER(lower)[i+n*j],ih=INTEGER(upper)[i+n*j];ch[j][il]+=ll.g[2*j]*mult[j];ch[j][ih]+=ll.g[2*j+1]*mult[j];if(scores){cch[j][ci][il]+=ll.g[2*j]*mult[j];cch[j][ci][ih]+=ll.g[2*j+1]*mult[j];}}
      if(association){double z=static_cast<double>(ll.g[2*d]*(theta-((family==1||family==4)?1:0)));REAL(grad)[ti]-=z;if(scores)REAL(score)[ci+m*ti]+=z;}
    }
    for(int j=0;j<d;++j){
      real offset=0;SEXP center=baseline_center==R_NilValue?R_NilValue:VECTOR_ELT(baseline_center,j);
      if(center!=R_NilValue)for(int k=0;k<p[j];++k)offset+=REAL(center)[k]*REAL(par)[bo[j]+k];
      real baseline_factor=std::exp(-offset),acc=0;
      for(int k=INTEGER(counts)[j]-1;k>=0;--k){acc+=ch[j][k+1];REAL(grad)[ho[j]+k]-=static_cast<double>(haz[j][k]*acc);objective+=bpen*haz[j][k]*baseline_factor;REAL(grad)[ho[j]+k]+=static_cast<double>(bpen*haz[j][k]*baseline_factor);}
      if(center!=R_NilValue)for(int k=0;k<p[j];++k)REAL(grad)[bo[j]+k]-=static_cast<double>(bpen*cum[j].back()*baseline_factor*REAL(center)[k]);
      if(scores)for(int i=0;i<m;++i){acc=0;for(int k=INTEGER(counts)[j]-1;k>=0;--k){acc+=cch[j][i][k+1];REAL(score)[i+m*(ho[j]+k)]=static_cast<double>(haz[j][k]*acc);}}
      for(int k=0;k<p[j];++k){double b=REAL(par)[bo[j]+k],penalty=ridge*(ridge_weights==R_NilValue?1:REAL(ridge_weights)[bo[j]+k]);objective+=0.5L*penalty*b*b;REAL(grad)[bo[j]+k]+=penalty*b;}
    }
    if(association){real distance=theta-((family==1||family==4)?1:0);objective+=0.5L*aridge*distance*distance;REAL(grad)[ti]+=static_cast<double>(aridge*distance*distance);}
    for(int k=0;k<q;++k)if(!std::isfinite(REAL(grad)[k]))throw std::runtime_error("Nonfinite likelihood gradient");
  }catch(const std::exception&){valid=false;objective=1e100L;for(int k=0;k<q;++k){objective+=REAL(par)[k]*REAL(par)[k];REAL(grad)[k]=2*REAL(par)[k];}}
  REAL(value)[0]=static_cast<double>(objective);SET_VECTOR_ELT(ans,0,value);SET_VECTOR_ELT(ans,1,grad);SET_VECTOR_ELT(ans,2,score);SET_VECTOR_ELT(ans,3,Rf_ScalarLogical(valid));names(ans,{"value","gradient","score","valid"});UNPROTECT(4);return ans;
}

extern "C" SEXP drw_cox_likelihood(SEXP par,SEXP X,SEXP stop,SEXP event,SEXP cluster,SEXP m_arg,SEXP ridge_arg,SEXP censor_arg,SEXP scores_arg){
  int n=rows(X),p=cols(X),m=Rf_asInteger(m_arg);bool censor=Rf_asLogical(censor_arg),scores=Rf_asLogical(scores_arg);bool individual_order=XLENGTH(censor_arg)==n;
  std::vector<double>times;for(int i=0;i<n;++i)if(INTEGER(event)[i])times.push_back(REAL(stop)[i]);std::sort(times.begin(),times.end());times.erase(std::unique(times.begin(),times.end()),times.end());int K=times.size();
  SEXP ans=PROTECT(Rf_allocVector(VECSXP,8));SEXP val=PROTECT(Rf_allocVector(REALSXP,1)),grad=PROTECT(Rf_allocVector(REALSXP,p)),info=PROTECT(Rf_allocMatrix(REALSXP,p,p)),base=PROTECT(Rf_allocVector(REALSXP,K)),bt=PROTECT(Rf_allocVector(REALSXP,K)),sc=PROTECT(Rf_allocMatrix(REALSXP,scores?m:0,scores?p:0)),bs=PROTECT(Rf_allocMatrix(REALSXP,scores?m:0,scores?K:0)),means=PROTECT(Rf_allocMatrix(REALSXP,scores?p:0,scores?K:0));
  std::fill(REAL(grad),REAL(grad)+p,0);std::fill(REAL(info),REAL(info)+p*p,0);if(scores){std::fill(REAL(sc),REAL(sc)+m*p,0);std::fill(REAL(bs),REAL(bs)+m*K,0);std::fill(REAL(means),REAL(means)+p*K,0);}
  std::vector<long double>eta(n);long double high=-std::numeric_limits<long double>::infinity();for(int i=0;i<n;++i){eta[i]=0;for(int k=0;k<p;++k)eta[i]+=REAL(X)[i+n*k]*REAL(par)[k];high=std::max(high,eta[i]);}
  std::vector<long double>w(n);for(int i=0;i<n;++i)w[i]=std::exp(eta[i]-high);long double objective=0;
  for(int r=0;r<K;++r){if(r%64==0)R_CheckUserInterrupt();double t=times[r];long double s0=0;int ne=0;std::vector<long double>s1(p,0),s2(p*p,0),ev(p,0),riski(m,0),eventi(m,0);
    for(int i=0;i<n;++i){bool first=individual_order?LOGICAL(censor_arg)[i]:(censor&&!INTEGER(event)[i]);bool risk=REAL(stop)[i]>t || (REAL(stop)[i]==t && !first);if(risk){s0+=w[i];riski[INTEGER(cluster)[i]-1]+=w[i];for(int j=0;j<p;++j){s1[j]+=w[i]*REAL(X)[i+n*j];for(int k=0;k<p;++k)s2[j+p*k]+=w[i]*REAL(X)[i+n*j]*REAL(X)[i+n*k];}}if(INTEGER(event)[i]&&REAL(stop)[i]==t){++ne;eventi[INTEGER(cluster)[i]-1]+=1;objective-=eta[i];for(int j=0;j<p;++j)ev[j]+=REAL(X)[i+n*j];}}
    if (!(s0 > 0)) {
      Rf_error("Empty risk set at an observed event time.");
    }
    objective+=ne*(std::log(s0)+high);long double jump=ne*std::exp(-high)/s0;REAL(base)[r]=static_cast<double>(jump);REAL(bt)[r]=t;
    for(int j=0;j<p;++j){long double mean=s1[j]/s0;REAL(grad)[j]+=static_cast<double>(ne*mean-ev[j]);for(int k=0;k<p;++k)REAL(info)[j+p*k]+=static_cast<double>(ne*(s2[j+p*k]/s0-mean*s1[k]/s0));if(scores)REAL(means)[j+p*r]=static_cast<double>(mean);}
    if(scores){for(int i=0;i<m;++i)REAL(bs)[i+m*r]=static_cast<double>(m*std::exp(-high)*(eventi[i]-ne*riski[i]/s0)/s0);
      for(int i=0;i<n;++i){bool first=individual_order?LOGICAL(censor_arg)[i]:(censor&&!INTEGER(event)[i]);bool risk=REAL(stop)[i]>t||(REAL(stop)[i]==t&&!first);long double residual=(INTEGER(event)[i]&&REAL(stop)[i]==t?1:0)-(risk?ne*w[i]/s0:0);int ci=INTEGER(cluster)[i]-1;for(int j=0;j<p;++j)REAL(sc)[ci+m*j]+=static_cast<double>(residual*(REAL(X)[i+n*j]-s1[j]/s0));}}
  }
  for(int j=0;j<p;++j){double ridge=REAL(ridge_arg)[XLENGTH(ridge_arg)==1?0:j];objective+=0.5L*ridge*REAL(par)[j]*REAL(par)[j];REAL(grad)[j]+=ridge*REAL(par)[j];REAL(info)[j+p*j]+=ridge;}
  REAL(val)[0]=static_cast<double>(objective);SET_VECTOR_ELT(ans,0,val);SET_VECTOR_ELT(ans,1,grad);SET_VECTOR_ELT(ans,2,info);SET_VECTOR_ELT(ans,3,bt);SET_VECTOR_ELT(ans,4,base);SET_VECTOR_ELT(ans,5,sc);SET_VECTOR_ELT(ans,6,bs);SET_VECTOR_ELT(ans,7,means);names(ans,{"value","gradient","information","times","hazard","score","baseline_score","risk_mean"});UNPROTECT(9);return ans;
}
