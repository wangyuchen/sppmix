#include "sppmix.h"
//Written by Sakis Micheas, 2015
//Operations on posterior realizations
//used by cpp and R functions
//visible in the package only

//' @export
// [[Rcpp::export]]
List GetStats_sppmix(vec const& gens,
                     double const& alpha)
{
  //apply burnin before calling this function
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

//' @export
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

//' @export
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

//' @export
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

//' @export
// [[Rcpp::export]]
List PostGenGetBestPerm_sppmix(List const& allgens)
{
//sppmix::PostGenGetBestPerm_sppmix(gens$allgens)
  int i,j,k,L=allgens.size();
  double mind,loss1;
  List permuted_gens(L),mix1=allgens[0];
  int m=mix1.size();
  double permnum=Factorial_sppmix(m);
  mat //stateperm=zeros(L,m),
    current_perm=zeros(L,m),//will contain the best perm
    previous_perm=ones(L,m);
//  vec oldindex(m),found(m);
  int minindex,done=0,count=0;
//  vec action(m);//estimating the ds hyperparam
//repeat until we find the L permutations
//that yield the smallest MC risk wrt choice of
//hyperparams ds for ps
//loss function is:-log(Dirichlet(ps,ds))
  while (done==0)
  {
    //squared error loss, minimized at mean(thetas)
    vec sum_p=zeros(m);
    mat sum_mu=zeros(m,2);
    mat sum_sigma=zeros(m,4);
    current_perm=previous_perm;
    for(i=0;i<L;i++)
    {
//      dDirichlet_sppmix(perm_ps,action);
      vec v1=current_perm.row(i).t(); //rPerm_sppmix(m);
//      current_perm.row(i)=v1.t();
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
//    Rcout << sum_p << std::endl ;
//    Rcout << sum_mu << std::endl ;
//    Rcout << sum_sigma << std::endl ;
    done=1;
//    Rcout << "passed1" << std::endl ;
    for(i=0;i<L;i++)
    {
//      vec v1=current_perm.row(i).t();
      vec cur_ps=GetRealiz_ps_sppmix(allgens,i);
      //      Rcout << cur_ps << std::endl ;
      mat cur_mus=GetRealiz_mus_sppmix(allgens,i);
      mat cur_sigmas=GetRealiz_sigmas_sppmix(allgens,i);
      //go through all permutations and find the one
      //that minimizes the loss for the given action
      //there are m! permutations...
      mind=1000000000000;
      for(k=0;k<permnum;k++)
      {
        vec v2=GetAPermutation_sppmix(m,k);
//        Rcout << v2 << std::endl ;
        vec perm_ps=Permute_vec_sppmix(cur_ps,v2);
        mat perm_mus=Permute_mat_sppmix(cur_mus,v2);
        mat perm_sigmas=Permute_mat_sppmix(cur_sigmas,v2);
        loss1=0;
        for(j=0;j<m;j++)
        {
          //        Rcout << as_scalar((sum_sigma.row(j)-perm_sigmas.row(j))*trans((sum_sigma.row(j)-perm_sigmas.row(j)))) << std::endl ;
          loss1=loss1+as_scalar((sum_sigma.row(j)-perm_sigmas.row(j))*trans((sum_sigma.row(j)-perm_sigmas.row(j))))
          +as_scalar((sum_mu.row(j)-perm_mus.row(j))*trans((sum_mu.row(j)-perm_mus.row(j))))
          +(sum_p(j)-perm_ps(j))*(sum_p(j)-perm_ps(j));
          //       sigmas(indthetas(i),stateperm(i,:),2,1))*(alphahat3(:,2,1)-sigmas(indthetas(i),stateperm(i,:),2,1))'+(alphahat3(:,1,1)-sigmas(indthetas(i),stateperm(i,:),1,1))*(alphahat3(:,1,1)-sigmas(indthetas(i),stateperm(i,:),1,1))'+(alphahat3(:,2,2)-sigmas(indthetas(i),stateperm(i,:),2,2))*(alphahat3(:,2,2)-sigmas(indthetas(i),stateperm(i,:),2,2))'+(alphahat3(:,1,2)-sigmas(indthetas(i),stateperm(i,:),1,2))*(alphahat3(:,1,2)-sigmas(indthetas(i),stateperm(i,:),1,2))'+(alphahat2(:,2)-mus(indthetas(i),stateperm(i,:),2))*(alphahat2(:,2)-mus(indthetas(i),stateperm(i,:),2))'+(alphahat2(:,1)-mus(indthetas(i),stateperm(i,:),1))*(alphahat2(:,1)-mus(indthetas(i),stateperm(i,:),1))'+(alphahat1-probs(indthetas(i),stateperm(i,:)))*(alphahat1-probs(indthetas(i),stateperm(i,:)))';
        }
        if(loss1<mind)
        {
          mind=loss1;
          minindex=k;
        }
      }
      current_perm.row(i)=GetAPermutation_sppmix(m,minindex).t();
    }

    for(i=0;i<L;i++)
    {
      for(j=0;j<m;j++)
      {
        if (current_perm(i,j)!=previous_perm(i,j))
        {
          done=0;
          break;
        }
      }
      if(done==0)break;
    }
//    oldindex=stateperm.row(minindex).t();
    previous_perm=current_perm;
    if(done)
    {
//      found=oldindex;
      break;
    }
    count++;
    printf("\rApplying Permutations, iteration %d",count);
    if (count>10000)
    {
      Rcout << "Didnt find it in 10000 iterations"<< std::endl ;
      break;
    }
  }
//  Rcout << "Best permutation:"<<current_perm//found.t()
//    << std::endl ;
//apply this permuation to all of the original
//realizations and return the permuted
// ps,mus and sigmas
  cube permgenmus=zeros(m,2,L);
  field<mat> permgensigmas(L,m);
  mat permgenps=zeros(L,m);

  for(i=0;i<L;i++)
  {
    vec cur_ps=GetRealiz_ps_sppmix(allgens,i);
    mat cur_mus=GetRealiz_mus_sppmix(allgens,i);
    mat cur_sigmas=GetRealiz_sigmas_sppmix(allgens,i);
    vec mu1(2),perm_ps=Permute_vec_sppmix(cur_ps,current_perm.row(i).t());
    mat perm_mus=Permute_mat_sppmix(cur_mus,current_perm.row(i).t());
    mat sig=zeros(2,2),perm_sigmas=Permute_mat_sppmix(cur_sigmas,current_perm.row(i).t());
//    vec mu1(2),perm_ps=Permute_vec_sppmix(cur_ps,found);
//    mat perm_mus=Permute_mat_sppmix(cur_mus,found);
//    mat sig=zeros(2,2),perm_sigmas=Permute_mat_sppmix(cur_sigmas,found);
    List mix2(m);
    for(j=0;j<m;j++)
    {
      sig(0,0)=perm_sigmas(j,0);
      sig(0,1)=perm_sigmas(j,1);
      sig(1,0)=perm_sigmas(j,2);
      sig(1,1)=perm_sigmas(j,3);
//     Rcout << perm_mus.row(j)<< std::endl ;
      mu1(0)=perm_mus(j,0);
      mu1(1)=perm_mus(j,1);
 //     Rcout << "passed2"<< std::endl ;
      mix2[j]=List::create(
        Named("p") = perm_ps(j),
        Named("mu") = mu1,
        Named("sigma") = sig);
      permgensigmas(i,j)=sig;
    }
    permuted_gens[i]=mix2;
    permgenps.row(i)=perm_ps.t();
    permgenmus.slice(i)=perm_mus;
  }
//apply the permutation to the realizations
//in the R code, just return the best one
  printf("\rDone                                                      \n");
  return List::create(
//    Named("best_perm") = found,
    Named("best_perm") = current_perm,
    Named("permuted_gens") = permuted_gens,
    Named("permuted_ps") = permgenps,
    Named("permuted_mus") = permgenmus,
    Named("permuted_sigmas") = permgensigmas);
//  return found;
}

//' @export
// [[Rcpp::export]]
List GetAllMeans_sppmix(List const& allgens,
                        int const& burnin)
{
//allgens is the list of lists of realizations
//  z=GetAllMeans_sppmix(gens$allgens,burnin)
  int i,j,L=allgens.size();
  List gen_realiz=allgens[0];
  int countgens=L-burnin,m=gen_realiz.size();
  vec sumps=zeros(m);
  mat summus=zeros(m,2);
  mat sumsigmas=zeros(m,4);

  for(i=burnin;i<L;i++)
  {
    sumps=sumps+GetRealiz_ps_sppmix(allgens,i);
    summus=summus+GetRealiz_mus_sppmix(allgens,i);
    sumsigmas=sumsigmas+GetRealiz_sigmas_sppmix(allgens,i);
//    countgens++;
  }
//  Rcout << countgens<< L-burnin<<std::endl ;
  mat sig(2,2);
  cube meansigmas(2,2,m);
  for(j=0;j<m;j++)
  {
    sig(0,0)=sumsigmas(j,0)/countgens;
    sig(0,1)=sumsigmas(j,1)/countgens;
    sig(1,0)=sumsigmas(j,2)/countgens;
    sig(1,1)=sumsigmas(j,3)/countgens;
    meansigmas.slice(j)=sig;
  }
  return List::create(
    Named("meanps") = sumps/countgens,
    Named("meanmus") = summus/countgens,
    Named("meansigmas") = meansigmas);
}

//' @export
// [[Rcpp::export]]
vec GetCompDistr_sppmix(vec const& numcomp,
    int const& maxnumcomp)
{
  //returns the distribution of the # of components
  //apply burnin before calling this function
  int L=numcomp.size();
  vec distr_numcomp(maxnumcomp);
//  vec newcomps=SubVec_sppmix(numcomp,burnin,L);
  //   numcomp.subvec(burnin,L);
  for(int j=0;j<maxnumcomp;j++)
  {
    uvec q1 = find(numcomp==j);
    //    Rcout << q1 << std::endl ;
    distr_numcomp(j)=1.0*q1.size()/L;
  }
  return distr_numcomp;
}

//' @export
// [[Rcpp::export]]
List GetBDCompRealiz_sppmix(List const& genBDmix,
    vec const& genlamdas,vec const& numcomp,
    int const& comp)
{
  //returns the gens for a specific # of components
  //apply burnin before calling this function
//  int L=numcomp.size();
  uvec indi=find(numcomp==comp);
  if(sum(indi)==0)//no realizations for this k
    return List::create();
  int nn1=indi.size();
  vec newlamdas=genlamdas(indi);
//  vec newnumcomp=numcomp(indi);
  List newgenmix(nn1);
  for(int i=0;i<nn1;i++)
    newgenmix[i]=genBDmix[indi(i)];

  return List::create(
    Named("newgenBD") = newgenmix,
    Named("newlamdas") = newlamdas);
    //,Named("newnumcomp") = newnumcomp);
}

//' @export
// [[Rcpp::export]]
mat GetAvgLabelsDiscrete2Multinomial_sppmix(mat
      const& genzs,int const& m)
{
  //returns the membership matrix from
  //the realizations matrix Lxn (the output
  //from DAMCMC)
  //apply burnin before calling this function
  int i,dat,iters=genzs.n_rows,n=genzs.n_cols;
  mat zmultinomial=zeros(n,m);
  vec ptavglabel=zeros(n);
  for(i=0;i<iters;i++)
  {
    //build the label matrix for this realization
    mat cur_z=zeros(n,m);
    //put 1 at location genzs(i,dat)
    for(dat=0;dat<n;dat++)
      cur_z(dat,genzs(i,dat))=1;
    zmultinomial+=cur_z;
//    ptavglabel+=genzs.row(i).t();
  }
  return zmultinomial/iters;
}
