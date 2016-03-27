#include "sppmix.h"
//Written by Sakis Micheas, 2015
//data in the form of points (x,y), or from an image
//hypers= 1:starting num of components, 2:tuning variable for kappa
//3:a=3, 4:g=.3, 5:gam=1, 6:l for kappa prior
//close all;[fitmus,fitsigmas,fitps,fitz,fitnumcomp,
//consts,invsigs]=BDMCMC2d(20,newdata2011,10000,2000,
//1,20,[15,.01,3,2,1,1],1);

//just get realizations no plotting here

//[[Rcpp::export]]
List BDMCMC2d_sppmix(int const& maxnumcomp,
                     mat const& data,
                     vec const& xlims,
                     vec const& ylims,
                     int const& L,
                     bool const& truncate,
                     double const& lamda,
                     double const& lamdab,
                     vec const& hyper)
{
  int n=data.n_rows;
  Rcout << "\nDataset has " << n <<" points, maximum # of components requested is "<<maxnumcomp << std::endl ;
  ivec numcomp(L);
  numcomp(0)=hyper(0);
  if (hyper(0)>maxnumcomp || hyper(0)<1)
    numcomp(0)=maxnumcomp;
  mat genz=zeros(L,n);
  cube genmus=zeros(maxnumcomp,2,L);
  //    Rcout << " passed" << std::endl ;
  //  genmus.slice(0), first realization
  //all elements in the field are 2x2 matrices
  field<mat> gensigmas(L,maxnumcomp),
    geninvsigmas(L,maxnumcomp);
  //  gensigmas(0,1), first realization for second comp
  mat genps=zeros(L,maxnumcomp),
    consts=zeros(L,maxnumcomp),approxmass=ones(L,maxnumcomp);
  vec meanps=zeros(maxnumcomp),
      propp=zeros(maxnumcomp);
  mat meanmus=zeros(maxnumcomp,2),
      propmus=zeros(maxnumcomp,2);
  cube meansigmas=zeros(2,2,maxnumcomp),
      propsig=zeros(2,2,maxnumcomp);
  //  Rcout << "meansigmas0="<<meansigmas.slice(0)<< std::endl ;

  int i,j,r,dat,comp,sum1;
  mat sumsig,propmu=zeros(1,2),
    genmutemp=zeros(1,2),
    mutemp1=zeros(1,2),propsigma,
    sumxmu,ps1,ps2,beta,cov1,betamat,betamatinv;
  vec approx=zeros(maxnumcomp),mu1=zeros(2),
    newmu=zeros(2);

  double birthprob,totdelta,sum11,sum22,Ly,
  //MHjump=0, approxcompj,
    ratio=1,val;//,quad;

  ivec which=randi(numcomp(0), distr_param(0,n-1));
  //process data for hyperparam values
  vec mins=zeros(2),maxs=zeros(2);
  if(xlims.size()==1)//passed 0, use data
  {
    mins(0)=min(data.col(0));
    mins(1)=min(data.col(1));
    maxs(0)=max(data.col(0));
    maxs(1)=max(data.col(1));
  }else
  {
    mins(0)=xlims(0);
    mins(1)=ylims(0);
    maxs(0)=xlims(1);
    maxs(1)=ylims(1);
  }
  //  Rcout << "mins="<<mins<< std::endl ;
  double Rx=max(data.col(0))-min(data.col(0)),
    Ry=max(data.col(1))-min(data.col(1));
  vec ksi=zeros(2,1);
  ksi(0)=sum(data.col(0))/n;
  ksi(1)=sum(data.col(1))/n;
//  Rcout << ksi<<"\n"<< std::endl ;
  mat kappa(2,2),kappainv(2,2);
  kappa(0,1)=0;
  kappa(1,0)=0;
  kappa(0,0)=hyper(1)/(Rx*Rx);
  kappa(1,1)=hyper(1)/(Ry*Ry);
  kappainv(0,1)=0;
  kappainv(1,0)=0;
  //  kappainv.eye(2,2);
  kappainv(0,0)=Rx*Rx/hyper(1);
  kappainv(1,1)=Ry*Ry/hyper(1);
//  Rcout << kappa<<"\n"<< std::endl ;
//  Rcout << kappainv<<"\n"<< std::endl ;
  //hypers 3:a=3, 4:g=.3, 5:gam=1
  double a=hyper(2),g=hyper(3),gam=hyper(4);
//  imat prevz,zmultinomial(n,maxnumcomp);
  mat hmat(2,2),hmatinv(2,2);
  hmat(0,1)=0;
  hmat(1,0)=0;
  //  hmat.eye(2,2);
  hmat(0,0)=100*g/(a*Rx*Rx);
  hmat(1,1)=100*g/(a*Ry*Ry);
  //  hmatinv.eye(2,2);
  hmatinv(0,1)=0;
  hmatinv(1,0)=0;
  hmatinv(0,0)=a*Rx*Rx/(100*g);
  hmatinv(1,1)=a*Ry*Ry/(100*g);
  //starting values
  for (i=0;i<numcomp(0);i++)
  {
    //    genmus.slice(0);
    genmus(i,0,0)=data(which(i),0);
    genmus(i,1,0)=data(which(i),1);
    gensigmas(0,i)=kappainv;
    geninvsigmas(0,i)=kappa;//invmat2d_sppmix(gensigmas(0,i));
    genps(0,i)=1.0/numcomp(0);
    consts(0,i)=1.0/sqrt(det(2*datum::pi*gensigmas(0,i)));
  }
//  Rcout << "passed1" << std::endl ;
  vec vec_of_ones=ones(numcomp(0));
  for (dat=0;dat<n;dat++)
  {
 //   Rcout << rDirichlet_sppmix(vec_of_ones)<< std::endl ;
    vec qij(numcomp(0));
    sum22=0;
    for(j=0;j<numcomp(0);j++)
    {
      mutemp1(0)=data(dat,0)-genmus(j,0,0);
      mutemp1(1)=data(dat,1)-genmus(j,1,0);
      sum11=1.0/sqrt(det(2*datum::pi*gensigmas(0,j)));
      val=sum11*exp(-.5*as_scalar(mutemp1*geninvsigmas(0,j)*trans(mutemp1)));
      qij(j)=genps(0,j)*val;
      sum22+=qij(j);
    }
//    sum11=sum(qij);//SumVec_sppmix(qij,0,numcomp(0)-1);
    if(sum22>0)
    {
//      qij=qij/sum11;
        genz(0,dat)=rDiscrete_sppmix(0,qij/sum22);
    }
    else
    {
      Rcout << "\nError: starting realization, a point not assigned a comp"<< std::endl ;
      return List::create();
    }
  }
  //  Rcout << sum(zmultinomial.col(0))<< std::endl ;
  //  Rcout << sum(zmultinomial.col(1))<< std::endl ;
  //Rcout << zmultinomial<< std::endl ;
  //    return List::create();
//  Rcout << "passed2" << std::endl ;
  //setup grid for truncation
/*  vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=mins(0)+i*(maxs(0)-mins(0))/(LL-1);
    ticsy(i)=mins(1)+i*(maxs(1)-mins(1))/(LL-1);
  }
  mat areas=zeros(LL,LL);
  for(j=0;j<LL-1;j++)
    for(i=0;i<LL-1;i++)
      areas(i,j)=(ticsx(i+1)-ticsx(i))*(ticsy(j+1)-ticsy(j));
*/
  int pt1,indInf=0;//,countiter=0;
  //start MCMC
  Rcout << "Preliminaries done. Starting Birth-Death MCMC" << std::endl ;
  for(i=0;i<L-1;i++)
  {
//    prevz=zmultinomial;
    printf("\rWorking: %3.1f%% complete",100.0*i/(L-2));
//compute death rates for each component
    for(j=0;j<numcomp(i);j++)
    {
      consts(i,j)=1.0/sqrt(det(2*datum::pi*gensigmas(i,j)));
      geninvsigmas(i,j)=invmat2d_sppmix(gensigmas(i,j));
    }
    vec deltas=zeros(numcomp(i));
    for(r=0;r<numcomp(i);r++)
    {
      Ly=1;
      for(dat=0;dat<n;dat++)
      {
        sum11=0;
        sum22=0;
        for(j=0;j<numcomp(i);j++)
        {
          mutemp1(0)=data(dat,0)-genmus(j,0,i);
          mutemp1(1)=data(dat,1)-genmus(j,1,i);
          val=consts(i,j)*exp(-.5*as_scalar(mutemp1*geninvsigmas(i,j)*trans(mutemp1)));
          sum11=sum11+genps(i,j)*val;
          if(j!=r)sum22=sum22+genps(i,j)*val;
        }
        if(sum11>0)
          Ly=Ly*sum22/sum11;
        else
        {
          indInf=r;
          Ly=datum::inf;
          break;
        }
      }
      deltas(r)=lamdab*Ly/lamda;
    }
    if(numcomp(i)>1)
    {
      totdelta=SumVec_sppmix(deltas,0,numcomp(i)-1);
      if(totdelta>0 && totdelta<datum::inf)
        deltas=deltas/totdelta;
      if(totdelta==datum::inf)
      {
        deltas=zeros(numcomp(i));
        deltas(indInf)=1;
      }
    }
    else
      totdelta=0;
///type of jump%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    birthprob=lamdab/(lamdab+totdelta);
//    Rcout << deltas << std::endl ;
    if(Rcpp::runif(1)[0]<birthprob
         && numcomp(i)<maxnumcomp)//do a birth'birth';
    {
 //     Rcout <<"\nbirth " <<i << std::endl ;
      numcomp(i+1)=numcomp(i)+1;
      propmus=genmus.slice(i);
      for(j=0;j<numcomp(i);j++)
        propsig.slice(j)=gensigmas(i,j);
//      Rcout << "passed1" << std::endl ;
      val=R::rbeta(1,numcomp(i));
      propp=(1-val)*genps.row(i).t();
      propp(numcomp(i+1)-1)=val;
      propmus.row(numcomp(i+1)-1)=rnorm2_sppmix(1,ksi,kappainv);
      betamat=rWishart_sppmix(2*g,.5*hmatinv);
      betamatinv=invmat2d_sppmix(betamat);
      ps1=rWishart_sppmix(2*a,.5*betamatinv);
      propsig.slice(numcomp(i+1)-1)=invmat2d_sppmix(ps1);
 //     Rcout << numcomp(i)<< " "<<numcomp(i+1) << std::endl ;
 //     Rcout << propp.subvec(0,numcomp(i+1)) << std::endl ;
 //     Rcout << propmus.submat(0,0,numcomp(i+1),1) << std::endl ;
    }
    else //do a death
    {
//      Rcout << "\ndeath "<<i << std::endl ;
      if(numcomp(i)>1)//'death';
      {
//        if(numcomp(i)==2)
 //         Rcout <<numcomp(i)<<" "<<maxnumcomp<< std::endl ;
        numcomp(i+1)=numcomp(i)-1;
//choose component to kill
//  Rcout <<numcomp(i)<<" "<<maxnumcomp<< std::endl ;
        if(numcomp(i)==maxnumcomp && totdelta<0.00001)
            //degenerate case choose one at random
        {
   //         Rcout <<"pt1"<< std::endl ;
            pt1=rDiscrete_sppmix(0,ones(numcomp(i))/numcomp(i));
   //         Rcout <<"after pt1"<< std::endl ;
            deltas(pt1)=1;
        }
   //     Rcout <<"numcomp(i) "<<numcomp(i)<<" deltas "<< deltas.size()<<std::endl ;
        comp=rDiscrete_sppmix(0,deltas);
   //     Rcout <<"after comp"<< std::endl ;
        propp=genps.row(i).t()/(1-genps(i,comp));
   //     Rcout <<"after propp "<<propp<<SubVec_sppmix(propp,comp,numcomp(i)-1)<< std::endl ;
//        propp=SubstituteVec_sppmix(propp,
//            SubVec_sppmix(propp,comp,numcomp(i)-1),comp);
//        Rcout <<"propp"<< propp << std::endl ;
  //       [ps(i,[1:1:numcomp(i)]~=comp)/(1-ps(i,comp))];
//        Rcout <<"comp="<<comp<<", numcomp(i+1)= "<< numcomp(i+1) << std::endl ;
        propmus=genmus.slice(i);
        for(j=0;j<numcomp(i);j++)
          propsig.slice(j)=gensigmas(i,j);
        for(j=comp;j<numcomp(i+1);j++)
        {
          propsig.slice(j)=propsig.slice(j+1);
          propmus.row(j)=propmus.row(j+1);
          propp(j)=propp(j+1);
        }
 //       Rcout <<"propmus"<< propmus << std::endl ;
      }
      else
      {
        Rcout << "Error: no components, death step fails"<< std::endl ;
        return List::create();
      }
    }
//    Rcout << "passed2" << std::endl ;
//Gibbs or MH steps for all variables
//using wishart for covariance matrices
//sample beta matrix
    sumsig=zeros(2,2);
    for (j=0;j<numcomp(i+1);j++)
    {
//      Rcout <<propsig.slice(j)<< std::endl ;
      geninvsigmas(i,j)=invmat2d_sppmix(propsig.slice(j));
//      Rcout <<geninvsigmas(i,j)<< std::endl ;
      consts(i,j)=1.0/sqrt(det(2*datum::pi*propsig.slice(j)));
//      Rcout <<consts(i,j)<< std::endl ;
      sumsig=sumsig+//invmat2d_sppmix(propsig.slice(j));
        geninvsigmas(i,j);
    }
  //  Rcout <<sumsig<< std::endl ;
    ps1=invmat2d_sppmix(2*hmat+2*sumsig);
    beta=rWishart_sppmix(2*g+2*numcomp(i+1)*a,ps1);
  //  Rcout <<beta<< std::endl ;
  //sample zi's:not multinomial, just discrete distribution much better
    for(dat=0;dat<n;dat++)
    {
      vec qij(numcomp(i+1));
      sum11=0;
      for(j=0;j<numcomp(i+1);j++)
      {
        mutemp1(0)=data(dat,0)-propmus(j,0);
        mutemp1(1)=data(dat,1)-propmus(j,1);
        qij(j)=propp(j)*consts(i,j)*exp(-.5*as_scalar(mutemp1*geninvsigmas(i,j)*trans(mutemp1)));
//        Rcout <<qij(j)<< std::endl ;
        sum11+=qij(j);
      }
  //    Rcout <<numcomp(i+1)<<" "<<dat<< std::endl ;
//      sum11=sum(qij);//SumVec_sppmix(qij,0,numcomp(i+1)-1);
      if(sum11>0)
      {
 //       qij=qij/sum11;
          genz(i+1,dat)=rDiscrete_sppmix(0,qij/sum11);
      }
      else
      {
        Rcout << "\nError: realiz "<<i<<", m="<<numcomp(i+1)<<", sum is "<<sum11<<", a point not assigned a comp\n"<< std::endl ;
        Rcout << "\npropp"<<propp<< std::endl ;
        Rcout << "\n"<<qij.t()<< std::endl ;
        return List::create();
      }
    }
// Rcout << "passed zs" << std::endl ;
//samples mus and sigmas
    vec ds=zeros(numcomp(i+1));
    for(j=0;j<numcomp(i+1);j++)
    {
//      Rcout << find(z.row(i+1)==j) << std::endl ;
      mat newdata=data.rows(find(genz.row(i+1)==j));
      sum1=newdata.n_rows;
//      Rcout << "passed newdata\n"<<newdata<<sum1 << std::endl ;
      //samples mus
      if(sum1>0)
      {
        newmu(0)=sum(newdata.col(0))/sum1;
        newmu(1)=sum(newdata.col(1))/sum1;
        //        Rcout << newmu << std::endl ;
      }
      else
      {
        newmu=ksi;
      }
//      Rcout << "newmu" <<newmu<< std::endl ;
//Rcout << "geninvsigmas(i,j)" <<geninvsigmas(i,j)<< std::endl ;
//      Rcout << "kappa" <<kappa<< std::endl ;
      cov1=invmat2d_sppmix(sum1*geninvsigmas(i,j)+kappa);
//      Rcout << "cov1" <<cov1<< std::endl ;
      mu1=cov1*(sum1*geninvsigmas(i,j)*newmu+kappa*ksi);
//      Rcout << "mu1" <<mu1<< std::endl ;
      genmutemp=rnorm2_sppmix(1,mu1,cov1);
//      Rcout << "genmu" <<genmutemp<< std::endl ;
      if(truncate)
      {
        mu1(0)=genmus(j,0,i);
        mu1(1)=genmus(j,1,i);
        ratio=ApproxMHRatiomu_sppmix(xlims,
          ylims,mu1,trans(genmutemp),
           gensigmas(i,j),sum1);
//          pow(ApproxMHRatiomu_sppmix(LL,xlims,ylims,
//           mu1,trans(genmutemp),gensigmas(i,j),
  //               geninvsigmas(i,j)),sum1);
      }
      else
        ratio=1;
      if(Rcpp::runif(1)[0]<ratio)
      {
        genmus(j,0,i+1)=genmutemp(0);
        genmus(j,1,i+1)=genmutemp(1);
      }
      else
      {
        genmutemp(0)=genmus(j,0,i);
        genmutemp(1)=genmus(j,1,i);
        genmus(j,0,i+1)=genmus(j,0,i);
        genmus(j,1,i+1)=genmus(j,1,i);
      }
//      Rcout << "passed mus" << std::endl ;
      //samples sigmas
      sumxmu=zeros(2,2);
      if (sum1>0)
        for(r=0;r<sum1;r++)
          sumxmu=sumxmu+(newdata.row(r)-genmutemp).t()*(newdata.row(r)-genmutemp);
//      Rcout << sum1<<"\n"<<std::endl;
//      Rcout << sumxmu<<"\n"<<std::endl;
      ps2=invmat2d_sppmix(2*beta+sumxmu);
      cov1=rWishart_sppmix(2*a+sum1,ps2);
      //      if(det(ps2)<=0){Rcout << "\n"<<        "ps2 not pd "<<det(ps2)<<"\n"<<        det(beta)<< " "<<det(sumxmu)<<std::endl;}
      propsigma=invmat2d_sppmix(cov1);
//     Rcout << "passed propsigma" <<propsigma<< std::endl ;
      if(truncate)
      {
        mu1(0)=genmus(j,0,i+1);
        mu1(1)=genmus(j,1,i+1);
        ratio=ApproxMHRatiosig_sppmix(xlims,
          ylims,mu1,gensigmas(i,j),
          propsigma,sum1);
 //         pow(ApproxMHRatiosig_sppmix(LL,xlims,ylims,mu1,
 //         propsigma,gensigmas(i,j),
  //            geninvsigmas(i,j)),sum1);
      }
      else
        ratio=1;
      if(Rcpp::runif(1)[0]<ratio)
      {
        geninvsigmas(i+1,j)=cov1;//invmat2d_sppmix(propsigma);
        gensigmas(i+1,j)=propsigma;
      }
      else
      {
        geninvsigmas(i+1,j)=geninvsigmas(i,j);
        gensigmas(i+1,j)=gensigmas(i,j);
      }
//      Rcout << "passed sigma" << std::endl ;
//      Rcout << gam<< std::endl ;
//      Rcout << sum1 << std::endl ;
      ds(j)=gam+sum1;
 //     Rcout << "passed ds "<<j << std::endl ;
    }
 //   Rcout << "passed sigmas" << std::endl ;
    //sample component probs
    if(numcomp(i+1)>1)
    {
      vec newps=rDirichlet_sppmix(ds);
  //    Rcout << newps<< std::endl ;
      for(j=0;j<numcomp(i+1);j++)
        genps(i+1,j)=newps(j);
    }
    else
      genps(i+1,0)=1;
//    Rcout << "passed ps" << std::endl ;
//    Rcout << gensigmas(i+1,0) << std::endl ;
    if(truncate)
      for(j=0;j<numcomp(i+1);j++)
      {
        mu1(0)=genmus(j,0,i+1);
        mu1(1)=genmus(j,1,i+1);
        approxmass(i,j)=
          ApproxBivNormProb_sppmix(xlims,
                                   ylims,mu1,gensigmas(i+1,j),2);
      }
  }
  printf("\rDone                                                      \n");
//  printf("\rMH acceptance %3.1f%%",100.0*MHjump/L);

  //create a list, with each element corresponding
  //to a single realization, which itself is a list
  //with each element containing the mixture ps, mus, sigs,
  //as a list of m elements
  List allgens(L);
  // ,mix(m);//list containing mixture ps,mus,sigs
  for(i=0;i<L;i++)
  {
    List mix(numcomp(i));
    for(j=0;j<numcomp(i);j++)
    {
      mu1(0)=genmus(j,0,i);
      mu1(1)=genmus(j,1,i);
      mix[j]=List::create(
        Named("p") = genps(i,j),
        Named("mu") = mu1,
        Named("sigma") = gensigmas(i,j));
    }
    allgens[i]=mix;
  }

  //sample lambdas
  double alamda=1,blamda=10000;
  vec lamdas=rgamma(L,n+alamda,1/(1+1/blamda));
  return List::create(
    Named("allgens_List") = allgens,
    Named("genps") = genps,
    Named("genmus") = genmus,
    Named("gensigmas") = gensigmas,
    Named("genzs") = genz,
    Named("genlamdas") = lamdas,
    Named("numcomp") = numcomp,
    Named("maxnumcomp") = maxnumcomp,
    Named("ApproxCompMass")=approxmass);
}
