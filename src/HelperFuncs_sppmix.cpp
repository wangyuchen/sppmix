#include "sppmix.h"
//helper functions only, used by cpp and R functions
//visible in the package only

// [[Rcpp::export]]
double Factorial_sppmix(int x)
{
  double num=1;
  while(x>1) num*=x--;
  return num;
}

// [[Rcpp::export]]
mat invmat2d_sppmix(mat const& A){
  mat B=zeros(2,2);
  double det=A(0,0)*A(1,1)-A(0,1)*A(1,0);
  //  if(det<=0){Rcout << "\n"<<"A is not pd " << std::endl;}
  B(0,0)=A(1,1)/det;
  B(0,1)=-A(0,1)/det;
  B(1,0)=-A(1,0)/det;
  B(1,1)=A(0,0)/det;
  return B;
}

// [[Rcpp::export]]
double densNormMixatx_sppmix(vec const& atx,List const& mix)
{
  //  densNormMix_sppmix(c(0,0),truemix)
  int j,m=mix.size();
  List mth_comp;
  double c1,val=0,psj;
  vec muk,mu1=zeros(2);
  mat invsigk,sigmak;
  //  Rcout << m<< std::endl ;
  //  Rcout << atx<< std::endl ;
  for(j=0;j<m;j++)
  {
    mth_comp = mix[j];
    psj=as<double>(mth_comp["p"]);
    muk = as<vec>(mth_comp["mu"]);
    //    Rcout << j<< std::endl ;
    //    Rcout << muk<< std::endl ;
    sigmak = as<mat>(mth_comp["sigma"]);
    //    Rcout << sigmak<< std::endl ;
    mu1(0)=atx(0)-muk(0);
    mu1(1)=atx(1)-muk(1);
    //    Rcout << "passed"<< std::endl ;
    invsigk=invmat2d_sppmix(sigmak);
    c1=1.0/sqrt(det(2*datum::pi*sigmak));
    val=val+psj*c1*
      exp(-.5*as_scalar(trans(mu1)*invsigk*mu1));
  }
  return val;
}

// [[Rcpp::export]]
mat dNormMix_sppmix(List const& mix, vec const& x,
                    vec const& y)
{
  int xnum=x.size(),ynum=y.size();
  //  int m = mix.size();
  //  List mth_comp;
  int i,j;
  mat z=zeros(xnum,ynum);
  vec atxy=zeros(2);
  for(i=0;i<xnum;i++)
    for(j=0;j<ynum;j++)
    {
      atxy(0)=x(i);
      atxy(1)=y(j);
      z(i,j)=densNormMixatx_sppmix(atxy,mix);
    }
    return z;
}

// [[Rcpp::export]]
vec Permute_vec_sppmix(vec const& oldvec,
                       vec const& perm)
{
  int p=perm.size();
  vec newvec(p);
  for (int i=0; i<p; i++)
    newvec(i)=oldvec(perm(i)-1);
  return newvec;
}


// [[Rcpp::export]]
mat Permute_mat_sppmix(mat const& oldmat,
                       vec const& perm)
{
  //permutes rows of a matrix
  int p=perm.size(),pp=oldmat.n_rows,q=oldmat.n_cols;
  if(p!=pp)
  {
    Rcout << "wrong dimensions" << std::endl ;
    return 0;
  }
  mat newmat(p,q);
  for (int i=0; i<p; i++)
    newmat.row(i)=oldmat.row(perm(i)-1);
  return newmat;
}

// [[Rcpp::export]]
mat GetAllPermutations_sppmix(int const& m)
{
  double permnum=Factorial_sppmix(m);
  mat allperms(permnum,m);
  int i,j;
  vec v(m);
  for(j=0;j<m;j++)v(j)=j+1;

  for (i=0;i<permnum;i++)
  {
    std::next_permutation(v.begin(),v.end());
    allperms.row(i)=v.t();
  }
  return allperms;
}

// [[Rcpp::export]]
vec GetAPermutation_sppmix(int const& m,int const& which)
{
  double permnum=Factorial_sppmix(m);
  vec aperm(m);
  int i,j;
  vec v(m);
  for(j=0;j<m;j++)v(j)=j+1;

  for (i=0;i<permnum;i++)
  {
    std::next_permutation(v.begin(),v.end());
    if (i==which)
    {
      aperm=v;
      break;
    }
  }
  return aperm;
}

//' @export
// [[Rcpp::export]]
List GetGrid_sppmix(int const& len,
            vec const& mins,vec const& maxs)
{
  //setup grid for truncation
int i,j;
vec ticsx=zeros(len),ticsy=zeros(len);
for (i=0;i<len;i++)
{
  ticsx(i)=mins(0)+i*(maxs(0)-mins(0))/(len-1);
  ticsy(i)=mins(1)+i*(maxs(1)-mins(1))/(len-1);
}
mat areas=zeros(len,len);
for(j=0;j<len-1;j++)
  for(i=0;i<len-1;i++)
    areas(i,j)=(ticsx(i+1)-ticsx(i))*(ticsy(j+1)-ticsy(j));

return List::create(
  Named("ticsx") = ticsx,
  Named("ticsy") = ticsy,
  Named("areas") = areas);
}

// [[Rcpp::export]]
bool EqVec_sppmix(vec const& v1,vec const& v2,
                  double const& tol)
{
  int len=v1.size();
  for(int i=0;i<len;i++)
  {
    if(abs(v1(i)-v2(i))>tol)
      return false;
  }
  return true;
}

// [[Rcpp::export]]
double logGammaFunc_sppmix(double const& x)
{
  // use approximation based on expansion of the Stirling series
  return .5*log(2*datum::pi)+(x-.5)*log(x)-x
  +1/(12*x)+1/(360*x*x*x)+1/(1260*x*x*x*x*x);
}

// [[Rcpp::export]]
double GammaFunc_sppmix(double const& x)
{
  return exp(logGammaFunc_sppmix(x));
}

// [[Rcpp::export]]
double dDirichlet_sppmix(vec const& ps,vec const& ds)
{
  int len=ps.size();
  double val=1;
  for(int i=0;i<len;i++)
    val=val*pow(ps(i),ds(i)-1)/GammaFunc_sppmix(ds(i));
  val=val*GammaFunc_sppmix(sum(ds));
  return val;
}



