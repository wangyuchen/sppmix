#include "sppmix.h"
//Operations on posterior realizations
//used by cpp and R functions
//visible in the package only

// [[Rcpp::export]]
List GetStats_sppmix(vec const& gens,
                     double const& alpha)
{
  double mu=mean(gens);
  int L=gens.size();
  vec CS(2);
  vec sortedgens=sort(gens);
  //  Rcout << sortedgens<< std::endl ;
  CS(1)=sortedgens(floor((1-alpha/2)*L));
  CS(0)=sortedgens(floor(alpha*L/2));
  return List::create(
    Named("Min") = gens.min(),
    Named("Mean") = mu,
    Named("Max") = gens.max(),
    Named("CredibleSet") = CS,
    Named("CredibleSetConfidence") = 100*(1-alpha));
}

// [[Rcpp::export]]
vec GetRealiz_ps_sppmix(List const& allgens,
                        int const& realiz)
{
  int j,L=allgens.size();
  if(realiz>=L)
  {
    Rcout << "index out of bounds" << std::endl ;
    return 0;
  }
  List mixcomp,gen_realiz=allgens[realiz];
  int m=gen_realiz.size();
  vec ps(m);
  for(j=0;j<m;j++)
  {
    mixcomp=gen_realiz[j];
    ps(j)=as<double>(mixcomp["p"]) ;
  }
  return ps;
}

// [[Rcpp::export]]
mat GetRealiz_mus_sppmix(List const& allgens,
                         int const& realiz)
{
  int j,L=allgens.size();
  if(realiz>=L)
  {
    Rcout << "index out of bounds" << std::endl ;
    return 0;
  }
  List mixcomp,gen_realiz=allgens[realiz];
  int m=gen_realiz.size();
  mat mus(m,2);
  for(j=0;j<m;j++)
  {
    mixcomp=gen_realiz[j];
    mus.row(j)= trans(as<vec>(mixcomp["mu"]));
  }
  return mus;
}

// [[Rcpp::export]]
mat GetRealiz_sigmas_sppmix(List const& allgens,
                            int const& realiz)
{
  int j,L=allgens.size();
  if(realiz>=L)
  {
    Rcout << "index out of bounds" << std::endl ;
    return 0;
  }
  List mixcomp,gen_realiz=allgens[realiz];
  int m=gen_realiz.size();
  //vectorize the matrix and return it
  mat sigmas(m,4);
  for(j=0;j<m;j++)
  {
    mixcomp=gen_realiz[j];
    mat sigmak = as<mat>(mixcomp["sigma"]);
    sigmas(j,0)=sigmak(0,0);
    sigmas(j,1)=sigmak(0,1);
    sigmas(j,2)=sigmak(1,0);
    sigmas(j,3)=sigmak(1,1);
    //  Rcout << muk<< std::endl ;
  }
  return sigmas;
}

// [[Rcpp::export]]
List PostGen_sppmix(List const& allgens)
{
  int i,L=allgens.size();
  List permuted_gens(L),mix1=allgens[0];
  int m=mix1.size();
  mat stateperm=zeros(L,m);
  mat oldindex=stateperm;
  int done=0;
  while (done==0)
  {
    //squared error loss, minimized at mean(thetas)
    vec sum_p=zeros(m);
    mat sum_mu=zeros(m,2);
    mat sum_sigma=zeros(m,4);

    for(i=0;i<L;i++)
    {
      vec v1=rPerm_sppmix(m);
      stateperm.row(i)=v1.t();
      //permute all generated p's,mu's and sigma's
      //and then get the average
      vec cur_ps=GetRealiz_ps_sppmix(allgens,i);
      //      Rcout << cur_ps << std::endl ;

      mat cur_mus=GetRealiz_mus_sppmix(allgens,i);
      mat cur_sigmas=GetRealiz_sigmas_sppmix(allgens,i);
      vec perm_ps=Permute_vec_sppmix(cur_ps,v1);
      mat perm_mus=Permute_mat_sppmix(cur_mus,v1);
      mat perm_sigmas=Permute_mat_sppmix(cur_sigmas,v1);
      sum_p=sum_p+perm_ps;
      sum_mu=sum_mu+perm_mus;
      sum_sigma=sum_sigma+perm_sigmas;
    }
    sum_p=sum_p/L;
    sum_mu=sum_mu/L;
    sum_sigma=sum_sigma/L;
    Rcout << sum_p << std::endl ;
    Rcout << sum_mu << std::endl ;
    Rcout << sum_sigma << std::endl ;
    done=1;
  }

  return List::create(
    Named("Permuted_gens") = permuted_gens);
}
