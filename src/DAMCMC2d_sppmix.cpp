#include "../inst/include/sppmix.h"
//data, x-y limits, m=num of comps to fit
//L=num iter, burnin
//function [marginal,meanlamda,meanmus,meansigmas,
//meanp,mus,sigmas,ps,meanz,sumz,postmeanintensity,
//approx,lamdas,invsigs]=DAMCMC2dNormalMixture(
//data,xlimits,ylimits,m,L,burnin,trueintensity,
//trueps,truemus,truesigmas,truncate,skiplots)
//just get realizations no plotting here

//[[Rcpp::export]]
List DAMCMC2d_sppmix(mat const& data,
                     vec const& xlims,
                     vec const& ylims,
                     int const& m,int const& L,
                     int const& burnin,int const& LL,
                     bool const& truncate)
{
  int n=data.n_rows;
  Rcout << "Dataset has " << n <<" points" << std::endl ;
  //  cube(n_rows, n_cols, n_slices)
  cube //genz1=zeros(n,m,L),
    genmus=zeros(m,2,L);
  //    Rcout << " passed" << std::endl ;
  //  genmus.slice(0), first realization
  //all elements in the field are 2x2 matrices
  field<mat> gensigmas(L,m),geninvsigmas(L,m);
  //  gensigmas(0,1), first realization for second comp
  mat meanz=zeros(n,m),genz=zeros(n,m),
    genps=zeros(L,m),
    consts=zeros(L,m);
  vec meanps=zeros(m);
  mat meanmus=zeros(m,2);
  cube meansigmas=zeros(2,2,m);
  //  Rcout << "meansigmas0="<<meansigmas.slice(0)<< std::endl ;

  int i,j,r,dat;
  ivec which=randi(m, distr_param(0,n-1));
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
  //  Rcout << ksi<< std::endl ;
  mat kappa(2,2),kappainv(2,2);
  kappa(0,1)=0;
  kappa(1,0)=0;
  kappa(0,0)=1/(Rx*Rx);
  kappa(1,1)=1/(Ry*Ry);
  kappainv(0,1)=0;
  kappainv(1,0)=0;
  //  kappainv.eye(2,2);
  kappainv(0,0)=Rx*Rx;
  kappainv(1,1)=Ry*Ry;
  //hypers 3:a=3, 4:g=.3, 5:gam=1
  double a=3,g=1,gam=1;
  imat prevz,zmultinomial(n,m);
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
  for (i=0;i<m;i++)
  {
    //    genmus.slice(0);
    genmus(i,0,0)=data(which(i),0);
    genmus(i,1,0)=data(which(i),1);
    gensigmas(0,i)=kappainv;
    geninvsigmas(0,i)=kappa;//invmat2d_sppmix(gensigmas(0,i));
    genps(0,i)=1.0/m;
    consts(0,i)=1.0/sqrt(det(2*datum::pi*gensigmas(0,i)));
  }
  for (dat=0;dat<n;dat++)
    zmultinomial.row(dat)=reshape(rMultinomial_sppmix(1,rDirichlet_sppmix(ones(m))),1,m);
  //  Rcout << sum(zmultinomial.col(0))<< std::endl ;
  //  Rcout << sum(zmultinomial.col(1))<< std::endl ;
  //Rcout << zmultinomial<< std::endl ;
  //    return List::create();
  mat sumsig,newdata,propmu=zeros(1,2),
    genmutemp=zeros(1,2),
    mutemp1=zeros(1,2),propsigma,
    sumxmu,ps1,ps2,beta,cov1;
  vec qij(m),approx=zeros(m),
    ds,qs=zeros(m),mu1=zeros(2),
    newmu=zeros(2);

  double MHjump=0,ratio=1,sumd,approxcompj;//,quad;
  //setup grid for truncation
  int sum1;
  vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=mins(0)+i*(maxs(0)-mins(0))/(LL-1);
    ticsy(i)=mins(1)+i*(maxs(1)-mins(1))/(LL-1);
  }
  mat areas=zeros(LL,LL);
  for(j=0;j<LL-1;j++)
    for(i=0;i<LL-1;i++)
      areas(i,j)=(ticsx(i+1)-ticsx(i))*(ticsy(j+1)-ticsy(j));

  int countiter=0;
  //start MCMC
  Rcout << "Preliminaries done. Starting MCMC" << std::endl ;
  for(i=0;i<L-1;i++)
  {
    prevz=zmultinomial;
    printf("\rWorking: %3.1f%% complete",100.0*i/(L-2));
    //printf("\rWorking: %3.1f%% complete (iteration %d/%d)",100.0*i/(L-2),i,L);
    //  sample B matrix
    sumsig=zeros(2,2);
    for (j=0;j<m;j++)
      sumsig=sumsig+geninvsigmas(i,j);
    ps1=invmat2d_sppmix(2*hmat+2*sumsig);
    //   if(det(ps1)<=0){Rcout << "\n"<<"ps1 not pd "<<det(ps1)<<"\n"<<std::endl;}
    beta=rWishart_sppmix(2*g+2*m*a,ps1);
    //samples mus and sigmas
    ds=zeros(m,1);
    for(j=0;j<m;j++)
    {
      //samples mus
      sum1=sum(zmultinomial.col(j));
      //     indi=find(zmultinomial.col(j)==1);
      if(sum1>0)
      {
        newdata=data.rows(find(zmultinomial.col(j)==1));
        //        Rcout << newdata<< std::endl ;
        //    if(newdata.n_rows!=sum1){Rcout << "\n"<<"not reading data properly "<<newdata.n_rows<< " "<<sum1 << std::endl ;return List::create();}
        newmu(0)=sum(newdata.col(0))/sum1;
        newmu(1)=sum(newdata.col(1))/sum1;
        //        Rcout << newmu << std::endl ;
      }
      else
      {
        //       break;
        // Rcout <<"\n"<< i << " "<< j << " "<<sum1 << std::endl ;
        //         return List::create();
        newmu=zeros(2,1);
        //        Rcout << "\n"<<"Component with less than 2 points, exiting" << std::endl ;
        //        return List::create();
      }
      cov1=invmat2d_sppmix(sum1*geninvsigmas(i,j)+kappa);
      //    if(det(cov1)<=0){Rcout << "\n"<<"Cov1 not pd, "<<cov1 << std::endl ;}
      //      Rcout << newmu << sum1<< std::endl ;
      mu1=cov1*(sum1*geninvsigmas(i,j)*newmu+kappa*ksi);
      //      Rcout << mu1<<cov1<< std::endl ;
      //proposed mu
      genmutemp=rnorm2_sppmix(1,mu1,cov1);
      //      Rcout << "genmu" <<genmutemp<< std::endl ;
      if(truncate)
      {
        mu1(0)=genmus(j,0,i);
        mu1(1)=genmus(j,1,i);
        ratio=pow(ApproxMHRatiomu_sppmix(LL,ticsx,ticsy,areas,
                                         mu1,trans(genmutemp),gensigmas(i,j),
                                         geninvsigmas(i,j)),sum1);
      }
      else
        ratio=1;
      if(rUnif_sppmix()<ratio)
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
      //samples sigmas
      sumxmu=zeros(2,2);
      if (sum1>0)
        for(r=0;r<sum1;r++)
          sumxmu=sumxmu+(newdata.row(r)-genmutemp).t()*(newdata.row(r)-genmutemp);
      //      Rcout << sumxmu<<"\n"<<std::endl;
      ps2=invmat2d_sppmix(2*beta+sumxmu);
      cov1=rWishart_sppmix(2*a+sum1,ps2);
      //      if(det(ps2)<=0){Rcout << "\n"<<        "ps2 not pd "<<det(ps2)<<"\n"<<        det(beta)<< " "<<det(sumxmu)<<std::endl;}
      propsigma=invmat2d_sppmix(cov1);
      if(truncate)
      {
        mu1(0)=genmus(j,0,i+1);
        mu1(1)=genmus(j,1,i+1);
        ratio=pow(ApproxMHRatiosig_sppmix(LL,ticsx,ticsy,areas,mu1,
                                          propsigma,gensigmas(i,j),
                                          geninvsigmas(i,j)),sum1);
      }
      else
        ratio=1;
      if(rUnif_sppmix()<ratio)
      {
        geninvsigmas(i+1,j)=cov1;//invmat2d_sppmix(propsigma);
        gensigmas(i+1,j)=propsigma;
      }
      else
      {
        geninvsigmas(i+1,j)=geninvsigmas(i,j);
        gensigmas(i+1,j)=gensigmas(i,j);
      }
      ds(j)=gam+sum1;
      if (truncate)
      {
        mu1(0)=genmus(j,0,i+1);
        mu1(1)=genmus(j,1,i+1);
        approxcompj=ApproxCompMass_sppmix(LL,ticsx,ticsy,
                                          areas,mu1,gensigmas(i+1,j),
                                          geninvsigmas(i+1,j));
      }
      else
        approxcompj=1;
      consts(i+1,j)=1.0/(approxcompj*
        sqrt(det(2*datum::pi*gensigmas(i+1,j))));
    }
    //sample component probs
    genps.row(i+1)=rDirichlet_sppmix(ds).t();
    //    Rcout << genps.row(i+1)<< std::endl ;
    //    Rcout << ds<< std::endl ;
    //   Rcout << "1 "<<sum(zmultinomial.col(0))<< std::endl ;
    //   Rcout << "2 "<<sum(zmultinomial.col(1))<< std::endl ;
    //    Rcout << consts.row(i+1)<< std::endl ;

    //sample indicators zij
    for(dat=0;dat<n;dat++)
    {
      sumd=0;
      for(j=0;j<m;j++)
      {
        //        genmutemp(0)=genmus(j,0,i+1);
        //       genmutemp(1)=genmus(j,1,i+1);
        //        Rcout << genmutemp<< std::endl ;
        //        Rcout << data.row(dat)<< std::endl ;
        //        mutemp1=data.row(dat)-genmutemp;
        mutemp1(0)=data(dat,0)-genmus(j,0,i+1);
        mutemp1(1)=data(dat,1)-genmus(j,1,i+1);
        //        quad=as_scalar(mutemp1*geninvsigmas(i+1,j)*trans(mutemp1));
        //        Rcout << quad<< std::endl ;
        //        Rcout << mutemp1<< std::endl ;
        qij(j)=genps(i+1,j)*consts(i+1,j)*
          exp(-.5*as_scalar(mutemp1*geninvsigmas(i+1,j)*trans(mutemp1)));
        sumd=sumd+qij(j);
        //        Rcout << qij(dat,j)<< std::endl ;
      }
      qs=qij/sumd;
      //if(sum(qs)<1)
      //Rcout << sum(qs)<< std::endl ;
      //      if (sum(qij)==1)
      zmultinomial.row(dat)=reshape(rMultinomial_sppmix(1,qs),1,m);
      //      zmultinomial.row(dat)=reshape(rMultinomial_sppmix(1,qij/sum(qij)),1,m);
      //      Rcout << sum(zmultinomial.row(dat))<< std::endl ;
    }
    ratio=1;
    for(j=0;j<m;j++)
    {
      //component with less that 2 points is not accepted
      if (sum(zmultinomial.col(j))<2)
      {
        //        Rcout <<i<<zmultinomial<< std::endl ;
        //       Rcout <<"\nComponent "<<j+1<<" has less than 2 pts"<< std::endl ;
        //       Rcout << "1 "<<sum(zmultinomial.col(0))<< std::endl ;
        //       Rcout << "2 "<<sum(zmultinomial.col(1))<< std::endl ;
        //       return List::create();
        ratio=0;
        break;
      }
    }
    if(rUnif_sppmix()<ratio)
    {
      MHjump=MHjump+1;
      if(i+1>burnin)
      {
        countiter++;
        genz=genz+zmultinomial;
        meanps=meanps+trans(genps.row(i+1));
        meanmus=meanmus+genmus.slice(i+1);
        for(j=0;j<m;j++)
        {
          meansigmas.slice(j)=meansigmas.slice(j)
          +gensigmas(i+1,j);
        }
      }
    }
    else
    {
      genps.row(i+1)=genps.row(i);
      for(j=0;j<m;j++)
      {
        geninvsigmas(i+1,j)=geninvsigmas(i,j);
        gensigmas(i+1,j)=gensigmas(i,j);
        genmus(j,0,i+1)=genmus(j,0,i);
        genmus(j,1,i+1)=genmus(j,1,i);
      }
      if(i+1>burnin)
      {
        countiter++;
        genz=genz+prevz;
        meanps=meanps+trans(genps.row(i+1));
        meanmus=meanmus+genmus.slice(i+1);
        for(j=0;j<m;j++)
        {
          meansigmas.slice(j)=meansigmas.slice(j)
          +gensigmas(i+1,j);
          //         mat a(2,2); a=meansigmas.slice(j);
          //         a=a+gensigmas(i+1,j); meansigmas.slice(j)=a;
        }
      }
      zmultinomial=prevz;
    }
  }
  printf("\rDone                                                      \n");
  printf("\rMH acceptance %3.1f%%",1.0*MHjump/L);
  //compute means
  meanz=genz/countiter;
  meanps=meanps/countiter;
  meanmus=meanmus/countiter;
  meansigmas=meansigmas/countiter;

  //  Rcout <<gensigmas(1,0)<< std::endl ;
  //compute the average of the posterior surfaces
  //needs to be multiplied by the lamdas
  mat AvgPostIntensity = zeros(LL,LL),
    PostIntensityAvg = zeros(LL,LL);
  List mix1(m);//list containing mixture ps,mus,sigs
  vec xy(2);
  double intensityatxy;
  countiter=0;
  for(int x1=0;x1<LL;x1++)
    for(int y1=0;y1<LL;y1++)
    {
      printf("\rComputing intensity surfaces: %3.1f%% complete",100.0*countiter/(LL*LL));
      xy(0)=ticsx(x1);
      xy(1)=ticsy(y1);
      //for each realization, compute the intensity
      //surface at the point xy
      for(i=burnin;i<L;i++)
      {
        for(j=0;j<m;j++)
        {
          mu1(0)=genmus(j,0,i);
          mu1(1)=genmus(j,1,i);
          mix1[j]=List::create(
            Named("p") = genps(i,j),
            Named("mu") = mu1,
            Named("sigma") = gensigmas(i,j));
        }
        intensityatxy=densNormMixatx_sppmix(xy,mix1);
        AvgPostIntensity(x1,y1)=AvgPostIntensity(x1,y1)+intensityatxy/(L-burnin);
      }
      //      AvgPostIntensity(x1,y1)=AvgPostIntensity(x1,y1)/(L-burnin);
      countiter++;
    }
    //create a list, with each element corresponding
    //to a single realization, which itself is a list
    //with each element contaning the mixture ps, mus, sigs,
    //as a list of m elements
    List allgens(L),
    mix(m);//list containing mixture ps,mus,sigs
  for(i=0;i<L;i++)
  {
    for(j=0;j<m;j++)
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

  return List::create(
    Named("allgens") = allgens,
    Named("genps") = genps,
    Named("genmus") = genmus,
    Named("gensigmas") = gensigmas,
    Named("meanz") = meanz,
    Named("meanps") = meanps,
    Named("meanmus") = meanmus,
    Named("meansigmas") = meansigmas,
    Named("ticsnum") = LL,
    Named("ticsx") = ticsx,
    Named("ticsy") = ticsy,
    Named("ticsareas") = areas,
    Named("AvgofPostIntensity") = AvgPostIntensity,
    Named("PostIntensityofAvg") = PostIntensityAvg);
}
