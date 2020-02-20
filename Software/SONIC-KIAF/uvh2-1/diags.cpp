double anisospace()
{
  double diff=0,sum=0;
  
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	diff+=e[sx][sy]*((sy-Middle)*(sy-Middle)-(sx-Middle)*(sx-Middle));
	sum+=e[sx][sy]*((sx-Middle)*(sx-Middle)+(sy-Middle)*(sy-Middle));
      }

  return diff/sum;
}


/*
double anisomomentum(gsl_interp_accel *pacc,gsl_interp_accel *cs2acc)
{
  double diff=0,sum=0;
  
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	diff+=(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))*(u[0][sx][sy]*u[0][sx][sy]-u[1][sx][sy]*u[1][sx][sy]);
	diff+=pixx[sx][sy]-piyy[sx][sy];
	sum+=(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))*(u[0][sx][sy]*u[0][sx][sy]+u[1][sx][sy]*u[1][sx][sy]);
	sum+=pixx[sx][sy]+piyy[sx][sy];
	sum+=2*eos(e[sx][sy],pacc,cs2acc);
      }
  return diff/sum;
}
*/

void allanisomomentum(double& ex,double& ep, double& totpx, double& totpy, double& e1c, double& e1s, double& e2c, double& e2s, double& e3c, double& e3s, double& eps1c, double&eps1s, double& eps2c, double& eps2s, double& eps3c, double& eps3s, double& eps4c, double& eps4s, double& eps5c, double& eps5s, double& eps6c, double& eps6s, double& eps7c, double& eps7s,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc)
{
    double diffp=0,sump=0;
    double diffx=0,sumx=0;
    totpx=0.0;
    totpy=0.0;
    double den2=0,den3=0;
    double num1c=0,num2c=0,num3c=0;
    double num1s=0,num2s=0,num3s=0;
    double num2 = 0.0, numeps21 = 0.0, numeps22 = 0.0, eden = 0.0;
    double eden3 = 0.0, eden4 = 0.0, eden5 = 0.0, eden6 = 0.0, eden7 = 0.0;
    double numeps11 = 0.0, numeps12 = 0.0, numeps71 = 0.0, numeps72 = 0.0;
    double numeps31 = 0.0, numeps32 = 0.0,numeps41 = 0.0, numeps42 = 0.0;
    double numeps51 = 0.0, numeps52 = 0.0,numeps61 = 0.0, numeps62 = 0.0;
  
  double avgx = 0.0, avgy = 0.0, norm = 0.0;
    for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
           double x = (sx-Middle)*AT/fmtoGeV; //in fm
           double y = (sy-Middle)*AT/fmtoGeV; //in fm
           avgx += x*e[sx][sy];
           avgy += y*e[sx][sy];
           norm += e[sx][sy];
      }
    avgx /= norm;
    avgy /= norm;

  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	//standard momentum anisotropy
	diffp+=(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))*(u[0][sx][sy]*u[0][sx][sy]-u[1][sx][sy]*u[1][sx][sy]);
	diffp+=pixx[sx][sy]-piyy[sx][sy];
	sump+=(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))*(u[0][sx][sy]*u[0][sx][sy]+u[1][sx][sy]*u[1][sx][sy]);
	sump+=pixx[sx][sy]+piyy[sx][sy];
	sump+=2*eos(e[sx][sy],pacc,cs2acc)-2*pib[sx][sy];
	
	//standard (not participant) spatial anisotropy
	diffx+=e[sx][sy]*((sy-Middle)*(sy-Middle)-(sx-Middle)*(sx-Middle));
	sumx+=e[sx][sy]*((sx-Middle)*(sx-Middle)+(sy-Middle)*(sy-Middle));
	
	//net transverse momentum
	double diff=0,sum=0;
	double ux = u[0][sx][sy];
	double uy = u[1][sx][sy];
	double ut=sqrt(1+ux*ux+uy*uy);
	double vx=ux/ut;
	double vy=uy/ut;
	double pitx=vx*pixx[sx][sy]+vy*pixy[sx][sy];
	double pity=vx*pixy[sx][sy]+vy*piyy[sx][sy];
	totpx += eos(e[sx][sy],pacc,cs2acc)*ut*ux + pitx;
	totpy += eos(e[sx][sy],pacc,cs2acc)*ut*ux + pity;
	
	//Teaney's momenum anisotropies
	double phiu = atan2(uy,ux);
	double pitt=pixx[sx][sy]+piyy[sx][sy];
	num2c += ut*(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))*(ux*ux+uy*uy)*cos(2*phiu);
	num2s += ut*(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))*(ux*ux+uy*uy)*sin(2*phiu);
	den2 += ut*(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))*(ux*ux+uy*uy) + ut*e[sx][sy];
	den3 += ut*ut*(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))*ut*ut + ut*ut*pitt;
	num1c += ut*(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))*sqrt(ux*ux+uy*uy)*ux;
	num1s += ut*(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))*sqrt(ux*ux+uy*uy)*uy;
	num3c += ut*(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))*pow(ux*ux+uy*uy,1.5)*cos(3*phiu);
	num3s += ut*(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))*pow(ux*ux+uy*uy,1.5)*sin(3*phiu);

	//general spatial anisotropies
	double x = (sx-Middle)*AT/fmtoGeV - avgx; //in fm
	double y = (sy-Middle)*AT/fmtoGeV - avgy; //in fm
        double phi = atan2(y,x);
        double r2 = x*x + y*y; 
        double r3 = pow(r2,1.5);
        double r4 = r2*r2;
        double r5 = r2*r3;
        double r6 = r2*r2*r2;
        double r7 = r4*r3;

        eden += r2*e[sx][sy];
        eden3 += r3*e[sx][sy];
        eden4 += r4*e[sx][sy];
        eden5 += r5*e[sx][sy];
        eden6 += r6*e[sx][sy];
        eden7 += r7*e[sx][sy];

        num2 += (y*y - x*x)*e[sx][sy];// standard eccentricity
        numeps11 += r3*cos(phi)*e[sx][sy];
        numeps12 += r3*sin(phi)*e[sx][sy];
        numeps21 += r2*cos(2*phi)*e[sx][sy];// for participant eccentricity
        numeps22 += r2*sin(2*phi)*e[sx][sy];// for psrticipant eccentricity
        numeps31 += r3*cos(3*phi)*e[sx][sy];// for triangularity with r^n weighting
        numeps32 += r3*sin(3*phi)*e[sx][sy];// for triangularity with r^n weighting
        numeps41 += r4*cos(4*phi)*e[sx][sy];
        numeps42 += r4*sin(4*phi)*e[sx][sy];
        numeps51 += r5*cos(5*phi)*e[sx][sy];
        numeps52 += r5*sin(5*phi)*e[sx][sy];
        numeps61 += r6*cos(6*phi)*e[sx][sy];
        numeps62 += r6*sin(6*phi)*e[sx][sy];
        numeps71 += r7*cos(7*phi)*e[sx][sy];
        numeps72 += r7*sin(7*phi)*e[sx][sy];

      }
  ep = diffp/sump;
  ex = diffx/sumx;
  totpx *= t/AT;
  totpy *= t/AT;
  e1c = num1c/den3;
  e1s = num1s/den3;
  e2c = num2c/den2;
  e2s = num2s/den2;
  e3c = num3c/den3;
  e3s = num3s/den3;
  eps1c = numeps11/eden3;
  eps1s = numeps12/eden3;
  eps2c = numeps21/eden;
  eps2s = numeps22/eden;
  eps3c = numeps31/eden3;
  eps3s = numeps32/eden3;
  eps4c = numeps41/eden4;
  eps4s = numeps42/eden4;
  eps5c = numeps51/eden5;
  eps5s = numeps52/eden5;
  eps6c = numeps61/eden6;
  eps6s = numeps62/eden6;
  eps7c = numeps71/eden7;
  eps7s = numeps72/eden7;
  
//   return diff/sum;
}

/*
//total transverse momentum per unit rapidity
double totalmomentumx(gsl_interp_accel *pacc,gsl_interp_accel *cs2acc)
{
  double momx=0.0;
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	double diff=0,sum=0;
	double ut=sqrt(1+u[0][sx][sy]*u[0][sx][sy]+u[1][sx][sy]*u[1][sx][sy]);
	double vx=u[0][sx][sy]/ut;
	double vy=u[1][sx][sy]/ut;
	double pitx=vx*pixx[sx][sy]+vy*pixy[sx][sy];
	double pity=vx*pixy[sx][sy]+vy*piyy[sx][sy];
	momx += eos(e[sx][sy],pacc,cs2acc)*ut*u[0][sx][sy] + pitx;
      }
  momx *= t/AT;
  return momx;
}
double totalmomentumy(gsl_interp_accel *pacc,gsl_interp_accel *cs2acc)
{
  double momy=0.0;
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	double diff=0,sum=0;
	double ut=sqrt(1+u[0][sx][sy]*u[0][sx][sy]+u[1][sx][sy]*u[1][sx][sy]);
	double vx=u[0][sx][sy]/ut;
	double vy=u[1][sx][sy]/ut;
	double pitx=vx*pixx[sx][sy]+vy*pixy[sx][sy];
	double pity=vx*pixy[sx][sy]+vy*piyy[sx][sy];
	momy += eos(e[sx][sy],pacc,cs2acc)*ut*u[1][sx][sy] + pity;
      }
  momy *= t/AT;
  return momy;
}
*/



double uphi(int sx,int sy)
{
  double temp=(u[1][sx][sy]*sx-u[0][sx][sy]*sy)/sqrt(sx*sx+sy*sy);
  return temp;
}

double ur(int sx,int sy)
{
  double temp=(u[0][sx][sy]*sx+u[1][sx][sy]*sy)/sqrt(sx*sx+sy*sy);
  return temp;
}

int foundit(int sx, int sy,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc)
{
  globali=geti(e[sx][sy]);
  globalx=getx(globali,e[sx][sy]);

  freeze_out << (sx-Middle)/fmtoGeV*AT << "\t";
  freeze_out << (sy-Middle)/fmtoGeV*AT << "\t";
  freeze_out << u[0][sx][sy] << "\t";
  freeze_out << u[1][sx][sy] << "\t";
  freeze_out << pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) << "\t";
  freeze_out << pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) << "\t";
  freeze_out << piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) << "\n";
  return 0;
}

void fancyfoundx(double sx, int sy,streambuf* pbuf,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc)
{

  //fstream str;

  streambuf* important;
  
  important=cout.rdbuf();


  cout.rdbuf(pbuf);

  cout << (sx-Middle)/fmtoGeV*AT << "\t";
  cout << (sy-Middle)/fmtoGeV*AT << "\t";

  double workhorse[NUMT];
  double iarr[NUMT];
  gsl_spline * hspl=gsl_spline_alloc (gsl_interp_cspline, NUMT);
  gsl_interp_accel * hsacc=gsl_interp_accel_alloc ();


  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=u[0][i][sy];
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sx,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=u[1][i][sy];
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sx,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=pixx[i][sy]/(e[i][sy]+eos(e[i][sy],pacc,cs2acc));
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sx,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=pixy[i][sy]/(e[i][sy]+eos(e[i][sy],pacc,cs2acc));
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sx,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=piyy[i][sy]/(e[i][sy]+eos(e[i][sy],pacc,cs2acc));
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sx,hsacc) << "\n";

  /*
  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=T(i,sy)/AT;
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sx,hsacc) << "\n";
  */

  gsl_spline_free (hspl);
  gsl_interp_accel_free (hsacc);

  cout.rdbuf(important);

}

void fancyfoundy(int sx, double sy,streambuf* pbuf,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc)
{

  streambuf* important;
  
  important=cout.rdbuf();


  //fstream str;

  cout.rdbuf(pbuf);

  cout << (sx-Middle)/fmtoGeV*AT << "\t";
  cout << (sy-Middle)/fmtoGeV*AT << "\t";

  double workhorse[NUMT];
  double iarr[NUMT];
  gsl_spline * hspl=gsl_spline_alloc (gsl_interp_cspline, NUMT);
  gsl_interp_accel * hsacc=gsl_interp_accel_alloc ();
  

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=u[0][sx][i];
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sy,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=u[1][sx][i];
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sy,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=pixx[sx][i]/(e[sx][i]+eos(e[sx][i],pacc,cs2acc));
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sy,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=pixy[sx][i]/(e[sx][i]+eos(e[sx][i],pacc,cs2acc));
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sy,hsacc) << "\t";

  for (int i=1;i<=NUMT;i++)
    {
      workhorse[i-1]=piyy[sx][i]/(e[sx][i]+eos(e[sx][i],pacc,cs2acc));
      iarr[i-1]=i;
    }

  gsl_spline_init (hspl, iarr, workhorse, NUMT);
  
  cout << gsl_spline_eval(hspl,sy,hsacc) << "\n";

  gsl_spline_free (hspl);
  gsl_interp_accel_free (hsacc);

  cout.rdbuf(important);
}


double rootfunction(double x,void *params)
{
  double temp;
  temp=gsl_spline_eval (workspline, x, wac);
  return (temp-TF);
}
/*
//fancy freeze-out with interpolation
void fancyfreeze(gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{
  int dflag=1;
  int debug=0;


  if (T(Middle,Middle,Tacc)>TF)
    dflag=0;

  //if yes, use whole lattice to determine 
  //freeze-out surface as accurate as possible
  //(important for corse lattices)
  if (dflag==0)
    {
      
      if (debug==1)
	printf("Starting freeze-out\n");

      double warr[Middle];
      double iarr[Middle];

      int limiter;
      int status;
      int iter = 0, max_iter = 1000;
      const  gsl_root_fsolver_type *TT;
      gsl_root_fsolver *ss;
      double r=0;

      gsl_function F;
      double dummy=0;

      F.function=&rootfunction;
      F.params=&dummy;

      TT = gsl_root_fsolver_brent;
      ss = gsl_root_fsolver_alloc (TT);

      //endpoint
      int ssy=Middle;
      int ssx=Middle;
      for (ssx=Middle;ssx<=NUMT;ssx++)
	{
	  globali=geti(e[ssx][ssy]);
	  globalx=getx(globali,e[ssx][ssy]);
	  warr[ssx-Middle]=T(ssx,ssy);
	  iarr[ssx-Middle]=ssx;
	  if (debug==1)
	    printf("%i %f\n",ssx,warr[ssx-Middle]);
	}

      gsl_spline_init (workspline,iarr,warr,Middle);

      gsl_root_fsolver_set (ss, &F, Middle, NUMT);
      iter=0;
      do
	{
	  iter++;
	  status = gsl_root_fsolver_iterate (ss);
	  r = gsl_root_fsolver_root (ss);
	  double x_lo = gsl_root_fsolver_x_lower (ss);
          double x_hi = gsl_root_fsolver_x_upper (ss);
	  status = gsl_root_test_interval (x_lo, x_hi,AT/20, 0.001);
         }
       while (status == GSL_CONTINUE && iter < max_iter);
      fancyfoundx(r,ssy,freeze_out.rdbuf(),pacc,cs2acc);

      //correct slightly to
      //inhibit bad conversion properties
      limiter=(int) (r-(AT/20));

      printf("freezel=%f\n",(limiter-Middle)/fmtoGeV*AT);

      //if we are very close to the finish
      //decrease step size to capture
      //important final effects
      //note:does not work so great
      //would need more complicated adjustment
      //if (((limiter-Middle)/fmtoGeV*AT)<1.5)
      //	{
      //	  //SNAPUPDATE=(int) (SNAPUPDATE/2.);
      //	  //SNAPUPDATE++;
      //	  UPDATE=(int) (UPDATE/1.2);
      //	  UPDATE++;
      //	}

      if (debug==1)
	printf("E1 done limit=%i, %f\n",limiter,r);



      //by scanning in y:
      //upper right quarter
      for (int sx=limiter;sx>=Middle;sx--)
	{
	  for (int sy=Middle;sy<=NUMT;sy++)
	    {
	      globali=geti(e[sx][sy]);
	      globalx=getx(globali,e[sx][sy]);
	      warr[sy-Middle]=T(sx,sy);
	      iarr[sy-Middle]=sy;
	    }
	  
	  gsl_spline_init (workspline,iarr,warr,Middle);

	  iter = 0;
	  gsl_root_fsolver_set (ss, &F, Middle, NUMT);
	  do
	    {
	      iter++;
	      status = gsl_root_fsolver_iterate (ss);
	      r = gsl_root_fsolver_root (ss);
	      double x_lo = gsl_root_fsolver_x_lower (ss);
	      double x_hi = gsl_root_fsolver_x_upper (ss);
	      status = gsl_root_test_interval (x_lo, x_hi,AT/20, 0.001);
	    }
	  while (status == GSL_CONTINUE && iter < max_iter);
	  fancyfoundy(sx,r,freeze_out.rdbuf(),pacc,cs2acc);

	  if (debug==1)
	    printf("fq: %i done\n",sx);
	}



      //upper left quarter
      for (int sx=Middle-1;sx>=2*Middle-limiter;sx--)
	{
	  for (int sy=Middle;sy<=NUMT;sy++)
	    {
	      globali=geti(e[sx][sy]);
	      globalx=getx(globali,e[sx][sy]);
	      warr[sy-Middle]=T(sx,sy);
	      iarr[sy-Middle]=sy;
	    }
	  
	  gsl_spline_init (workspline,iarr,warr,Middle);
	  
	  gsl_root_fsolver_set (ss, &F, Middle, NUMT);
	  iter=0;
	  do
	    {
	      iter++;
	      status = gsl_root_fsolver_iterate (ss);
	      r = gsl_root_fsolver_root (ss);
	      double x_lo = gsl_root_fsolver_x_lower (ss);
	      double x_hi = gsl_root_fsolver_x_upper (ss);
	      status = gsl_root_test_interval (x_lo, x_hi,AT/20, 0.001);
	    }
	  while (status == GSL_CONTINUE && iter < max_iter);
	  fancyfoundy(sx,r,freeze_out.rdbuf(),pacc,cs2acc);
	  if (debug==1)
	    printf("sq: %i done\n",sx);
	}

      
      //endpoint
      ssy=Middle;
      ssx=Middle;
      for (ssx=1;ssx<=Middle;ssx++)
	{
	  globali=geti(e[ssx][ssy]);
	  globalx=getx(globali,e[ssx][ssy]);
	  warr[ssx-1]=T(ssx,ssy);
	  iarr[ssx-1]=ssx;
	}

      gsl_spline_init (workspline,iarr,warr,Middle);

      gsl_root_fsolver_set (ss, &F, 1,Middle);
      iter=0;
      do
	{
	  iter++;
	  status = gsl_root_fsolver_iterate (ss);
	  r = gsl_root_fsolver_root (ss);
	  double x_lo = gsl_root_fsolver_x_lower (ss);
          double x_hi = gsl_root_fsolver_x_upper (ss);
	  status = gsl_root_test_interval (x_lo, x_hi,AT/20, 0.001);
	}
      while (status == GSL_CONTINUE && iter < max_iter);      
      fancyfoundx(r,ssy,freeze_out.rdbuf(),pacc,cs2acc);

      if (debug==1)
	printf("E2 done %i, %f\n",2*Middle-limiter,r);


      //lower left quarter
      for (int sx=2*Middle-limiter;sx<=Middle;sx++)
	{
	  for (int sy=1;sy<=Middle;sy++)
	    {
	      globali=geti(e[sx][sy]);
	      globalx=getx(globali,e[sx][sy]);
	      warr[sy-1]=T(sx,sy);
	      iarr[sy-1]=sy;
	    }

	  gsl_spline_init (workspline,iarr,warr,Middle);
	  
	  gsl_root_fsolver_set (ss, &F, 1,Middle);
	  iter=0;
	  do
	    {
	      iter++;
	      status = gsl_root_fsolver_iterate (ss);
	      r = gsl_root_fsolver_root (ss);
	      double x_lo = gsl_root_fsolver_x_lower (ss);
	      double x_hi = gsl_root_fsolver_x_upper (ss);
	      status = gsl_root_test_interval (x_lo, x_hi,AT/20, 0.001);
	    }
	  while (status == GSL_CONTINUE && iter < max_iter);
	  fancyfoundy(sx,r,freeze_out.rdbuf(),pacc,cs2acc);
	  if (debug==1)
	    printf("tq: %i done\n",sx);
	}

      //lower right quarter
      for (int sx=Middle+1;sx<=limiter;sx++)
	{
	  for (int sy=1;sy<=Middle;sy++)
	    {
	      globali=geti(e[sx][sy]);
	      globalx=getx(globali,e[sx][sy]);
	      warr[sy-1]=T(sx,sy);
	      iarr[sy-1]=sy;
	    }

	  gsl_spline_init (workspline,iarr,warr,Middle);
	  
	  gsl_root_fsolver_set (ss, &F, 1,Middle);
	  iter=0;
	  do
	    {
	      iter++;
	      status = gsl_root_fsolver_iterate (ss);
	      r = gsl_root_fsolver_root (ss);
	      double x_lo = gsl_root_fsolver_x_lower (ss);
	      double x_hi = gsl_root_fsolver_x_upper (ss);
	      status = gsl_root_test_interval (x_lo, x_hi,AT/20, 0.001);
	    }
	  while (status == GSL_CONTINUE && iter < max_iter);
	  fancyfoundy(sx,r,freeze_out.rdbuf(),pacc,cs2acc);
	  if (debug==1)
	    printf("4q: %i done\n",sx);
	}

      
      //cout << " Found " << r << endl;


      
      gsl_root_fsolver_free (ss);
      
      freeze_out << "TIME \t" << t/fmtoGeV*AT << endl;
    }
  else
    reachedTf=1;
}
*/

void stupidfreeze(gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{

  if (T(Middle,Middle,Tacc)<TF)
    {
      reachedTf=1;
      
      for (int sx=1;sx<=NUMT;sx++)
	for (int sy=1;sy<=NUMT;sy++)
	  {
	    freeze_out << (sx-Middle)/fmtoGeV*AT << "\t";
	    freeze_out << (sy-Middle)/fmtoGeV*AT << "\t";
	    freeze_out << u[0][sx][sy] << "\t";
	    freeze_out << u[1][sx][sy] << "\t";
	    freeze_out << pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) << "\t";
	    freeze_out << pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) << "\t";
	    freeze_out << piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) << "\t";
	    freeze_out << T(sx,sy,Tacc)/AT << "\n";
	  }
      freeze_out << "TIME \t" << t/fmtoGeV*AT << endl;
    }
}



//freeze-out that doesn't assume a monotonically 
//decreasing temperature from the center out
//written by Matthew Luzum
void blockfreeze(gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{


  //All grid points that are nearest neighbors in x, y, or tau are checked pairwise.
  //A rectangular piece of the freezeout surface is defined to be half-way between
  //any pair whose temperatures straddle TF.
  
  reachedTf = 1; //hydro evolution is finished if everywhere in the time slice T < TF
  for (int sx=1;sx<=NUMT;sx++) {
    for (int sy=1;sy<=NUMT;sy++)
    {
      
      if (T(sx,sy,Tacc) > TF) reachedTf = 0; //if one cell has T>TF, freezeout not reached
      
      if (sy != NUMT)
      {
	if (((T(sx,sy,Tacc) > TF) && (T(sx,sy+1,Tacc) <= TF)) || ((T(sx,sy,Tacc) <= TF) && (T(sx,sy+1,Tacc) > TF)))
	{ 
	  reachedTf = 0;
	  int direction = 2;
	  if ((T(sx,sy,Tacc) <= TF) && (T(sx,sy+1,Tacc) > TF)) direction = -2;
	  freeze_out << (sx-Middle)/fmtoGeV*AT << "\t";
	  freeze_out << (sy-Middle+0.5)/fmtoGeV*AT << "\t";
	  freeze_out << t/fmtoGeV*AT << "\t";
	  freeze_out << direction << "\t";
	  freeze_out << 0.5 * (u[0][sx][sy] + u[0][sx][sy+1]) << "\t";
	  freeze_out << 0.5 * (u[1][sx][sy] + u[1][sx][sy+1]) << "\t";
	  freeze_out << 0.5 * (pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) 
			      + pixx[sx][sy+1]/(e[sx][sy+1]+eos(e[sx][sy+1],pacc,cs2acc))) << "\t";
	  freeze_out << 0.5 * (pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) 
			      + pixy[sx][sy+1]/(e[sx][sy+1]+eos(e[sx][sy+1],pacc,cs2acc)) ) << "\t";
	  freeze_out << 0.5 * (piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))
			      + piyy[sx][sy+1]/(e[sx][sy+1]+eos(e[sx][sy+1],pacc,cs2acc)))<< "\t";
	  freeze_out << 0.5 * (pib[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) 
			      + pib[sx][sy+1]/(e[sx][sy+1]+eos(e[sx][sy+1],pacc,cs2acc)) ) << "\t";
	  freeze_out << 0.5 * (T(sx,sy,Tacc) + T(sx,sy+1,Tacc))/AT << "\n";
	}
      }
      if (sx != NUMT)
      {
      	if (((T(sx,sy,Tacc) > TF) && (T(sx+1,sy,Tacc) <= TF)) || ((T(sx,sy,Tacc) <= TF) && (T(sx+1,sy,Tacc) > TF)))
	{ 
	  reachedTf = 0;
	  int direction = 1;
	  if ((T(sx,sy,Tacc) <= TF) && (T(sx+1,sy,Tacc) > TF)) direction = -1;
	  freeze_out << (sx-Middle+0.5)/fmtoGeV*AT << "\t";
	  freeze_out << (sy-Middle)/fmtoGeV*AT << "\t";
	  freeze_out << t/fmtoGeV*AT << "\t";
	  freeze_out << direction << "\t";
	  freeze_out << 0.5 * (u[0][sx][sy] + u[0][sx+1][sy]) << "\t";
	  freeze_out << 0.5 * (u[1][sx][sy] + u[1][sx+1][sy]) << "\t";
	  freeze_out << 0.5 * (pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) 
			      + pixx[sx+1][sy]/(e[sx+1][sy]+eos(e[sx+1][sy],pacc,cs2acc))) << "\t";
	  freeze_out << 0.5 * (pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) 
			      + pixy[sx+1][sy]/(e[sx+1][sy]+eos(e[sx+1][sy],pacc,cs2acc)) ) << "\t";
	  freeze_out << 0.5 * (piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))
			      + piyy[sx+1][sy]/(e[sx+1][sy]+eos(e[sx+1][sy],pacc,cs2acc)))<< "\t";
	  freeze_out << 0.5 * (pib[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))
			      + pib[sx+1][sy]/(e[sx+1][sy]+eos(e[sx+1][sy],pacc,cs2acc)))<< "\t";
	  freeze_out << 0.5 * (T(sx,sy,Tacc) + T(sx+1,sy,Tacc))/AT << "\n";
	}
      }

      if (t>TINIT*fmtoGeV/AT)
      {
	if ((T(sx,sy,Tacc) <= TF && (Tlast(sx,sy,Tacc) > TF)) || ((T(sx,sy,Tacc) > TF) && (Tlast(sx,sy,Tacc) <= TF)))
	{
	  int direction = 3;
	  if ((T(sx,sy,Tacc) > TF) && (Tlast(sx,sy,Tacc) <= TF)) direction = -3;
	  freeze_out << (sx-Middle)/fmtoGeV*AT << "\t";
	  freeze_out << (sy-Middle)/fmtoGeV*AT << "\t";
	  freeze_out << t/fmtoGeV*AT - 0.5 * UPDATE*EPS*AT/fmtoGeV << "\t";
	  freeze_out << direction << "\t";
	  freeze_out << 0.5 * (u[0][sx][sy] + ulast[0][sx][sy]) << "\t";
	  freeze_out << 0.5 * (u[1][sx][sy] + ulast[1][sx][sy]) << "\t";
	  freeze_out << 0.5 * (pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) 	
			      + pixxlast[sx][sy]/(elast[sx][sy]+eos(elast[sx][sy],pacc,cs2acc))) << "\t";
	  freeze_out << 0.5 * (pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) 
			      + pixylast[sx][sy]/(elast[sx][sy]+eos(elast[sx][sy],pacc,cs2acc)) ) << "\t";
	  freeze_out << 0.5 * (piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))
			      + piyylast[sx][sy]/(elast[sx][sy]+eos(elast[sx][sy],pacc,cs2acc)))<< "\t";
	  freeze_out << 0.5 * (pib[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))
			      + pilast[sx][sy]/(elast[sx][sy]+eos(elast[sx][sy],pacc,cs2acc)))<< "\t";
	  freeze_out << 0.5 * (T(sx,sy,Tacc) + Tlast(sx,sy,Tacc))/AT << "\n";
	}
      }
    }
  }
}

void outputMeasurements(double t,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc) 
{
  cout.precision(5);
  int dwidth = 13;
  cout.width(dwidth); cout << t/fmtoGeV*AT;
  globali=geti(e[Middle][Middle]);
  globalx=getx(globali,e[Middle][Middle]);
//   double ex=anisospace();
//   double ep=anisomomentum();
//   double px = totalmomentumx();
//   double py = totalmomentumy();
  double ex, ep, totpx, totpy, e1c, e1s, e2c, e2s, e3c, e3s;
  double eps1c, eps1s, eps2c, eps2s, eps3c, eps3s, eps4c, eps4s;
  double eps5c, eps5s, eps6c, eps6s, eps7c, eps7s;
  allanisomomentum(ex,ep, totpx, totpy, e1c, e1s, e2c, e2s, e3c, e3s,
		  eps1c, eps1s, eps2c, eps2s, eps3c, eps3s, eps4c, 
		   eps4s, eps5c, eps5s, eps6c, eps6s, eps7c, eps7s,pacc,cs2acc);

  cout.width(dwidth); cout << T(Middle,Middle,Tacc)/AT;
  cout << "\t" << ex << "\t" << ep;
  //cout << "\t" << overlapS()/AT/AT;
  //cout << "\t" << l1coeff(T(Middle,Middle,Tacc),Tacc);  
  cout << endl;

  meta << t/fmtoGeV*AT <<"\t";
  meta << T(Middle,Middle,Tacc)/AT << "\t";
  meta << e[Middle][Middle]/AT/AT/AT/AT << "\t";
  meta << -pi(3,3,Middle,Middle)*t*t/AT/AT/AT/AT << "\t";
  meta << pib[Middle][Middle]/AT/AT/AT/AT << "\t";
  meta << "\n";
  
  ecces << t/fmtoGeV*AT <<"\t";
  ecces << ex << "\t";
  ecces << ep << "\t";
  ecces << totpx << "\t";
  ecces << totpy << "\t";
  ecces << e1c << "\t";
  ecces << e1s << "\t";
  ecces << e2c << "\t";
  ecces << e2s << "\t";
  ecces << e3c << "\t";
  ecces << e3s << "\t";
  ecces << eps1c << "\t";
  ecces << eps1s << "\t";
  ecces << eps2c << "\t";
  ecces << eps2s << "\t";
  ecces << eps3c << "\t";
  ecces << eps3s << "\t";
  ecces << eps4c << "\t";
  ecces << eps4s << "\t";
  ecces << eps5c << "\t";
  ecces << eps5s << "\t";
  ecces << eps6c << "\t";
  ecces << eps6s << "\t";
  ecces << eps7c << "\t";
  ecces << eps7s << "\n";

    //freezeout();
//   if (FREEZE==1)
//     fancyfreeze();
//   else
//     stupidfreeze();
    
  switch (FREEZE){
    case 0:
      stupidfreeze(pacc,cs2acc,Tacc);
      break;
      /*    case 1:
      fancyfreeze(pacc,cs2acc);
      break;*/
    case 2:
      blockfreeze(pacc,cs2acc,Tacc);
      break;
    default:
      blockfreeze(pacc,cs2acc,Tacc);
  }
}

void snapTprofile(double time,gsl_interp_accel *Tacc)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/Tprofile_%.2f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    globali=geti(e[Middle][s]);
    globalx=getx(globali,e[Middle][s]);
    out << (s-Middle)/fmtoGeV*AT << "\t";
    out << T(Middle,s,Tacc)/AT << endl;
  }
  out.close();
}

void snapEDprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/EDprofile_%.2f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    
    out << (s-Middle)/fmtoGeV*AT << "\t";
    out << e[Middle][s]/AT/AT/AT/AT << endl;
  }
  out.close();
}

void snapTcontour(double time,gsl_interp_accel *Tacc)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/Tcontour_%.3f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);

   out << "#x [fm] \t y [fm] \t T [GeV] \n";
  for (int sy=1;sy<=NUMT;sy++)
    {
      for (int sx=1;sx<=NUMT;sx++)
	{
	  globali=geti(e[sx][sy]);
	  globalx=getx(globali,e[sx][sy]);
	  out << (sx-Middle)/fmtoGeV*AT;
	  out << "\t";
	  out << (sy-Middle)/fmtoGeV*AT;
	  out << "\t";
	  out << T(sx,sy,Tacc)/AT << "\n";
	}
      out << endl;
    }
  out.close();
}

void snapVcontour(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/Vcontour_%.3f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  out << "#x [fm] \t y [fm] \t ux \t uy \t gamma \n";
  for (int sy=1;sy<=NUMT;sy++)
    {
      for (int sx=1;sx<=NUMT;sx++)
	{
	  globali=geti(e[sx][sy]);
	  globalx=getx(globali,e[sx][sy]);
	  out << (sx-Middle)/fmtoGeV*AT << "\t";
	  out << (sy-Middle)/fmtoGeV*AT << "\t";
	  out << u[0][sx][sy] << "\t";
	  out << u[1][sx][sy] << "\t";
	  out << sqrt(u[0][sx][sy]*u[0][sx][sy]+u[1][sx][sy]*u[1][sx][sy]) <<"\n";
	}
    }
  out.close();
}

void snapFOdata(double time,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/FOdata_%.3f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int sy=1;sy<=NUMT;sy++)
    {
      for (int sx=1;sx<=NUMT;sx++)
	{
	  out << u[0][sx][sy] << "\t";
	  out << u[1][sx][sy] << "\t";
	  out << pixx[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) << "\t";
	  out << pixy[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) << "\t";
	  out << piyy[sx][sy]/(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)) << "\t";
	}
      out << endl;
    }
  out.close();
  //putting names into 
  //file to facilitate later 
  //freeze-out calc
  meta << fname << "\n";
}

void snapVprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/Vprofile_%.2f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    globali=geti(e[Middle][s]);
    globalx=getx(globali,e[Middle][s]);
    out << (s-Middle)/fmtoGeV*AT << "\t";
    out << u[1][Middle][s]/ut(Middle,s) << endl;
  }
  out.close();
}

void snapVxprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/Vxprofile_%.2f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    globali=geti(e[s][Middle]);
    globalx=getx(globali,e[s][Middle]);
    out << (s-Middle)/fmtoGeV*AT << "\t";
    out << u[0][s][Middle]/ut(s,Middle) << endl;
  }
  out.close();
}

void snapV2profile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/V2profile_%.2f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int s=0;s<NUMT/2./sqrt(2.);s++)
  {
    globali=geti(e[Middle+s][Middle+s]);
    globalx=getx(globali,e[Middle+s][Middle+s]);
    out << s*sqrt(2)/fmtoGeV*AT << "\t";
    out << (u[1][Middle+s][Middle+s]+u[0][Middle+s][Middle+s])/sqrt(2.)/ut(Middle+s,Middle+s) << endl;
  }
  out.close();
}

void snappieeprofile(double time)
{
  fstream out;
  char fname[255];
  double vr=0;
  sprintf(fname,"data/snapshot/Piprofile_%.2f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    globali=geti(e[Middle][s]);
    globalx=getx(globali,e[Middle][s]);
    out << (s-Middle)/fmtoGeV*AT << "\t";
    out << -pi(3,3,Middle,s)*t*t/(e[Middle][s]*4./3.)<< endl;
  }
  out.close();
}

void snappirrprofile(double time)
{
  fstream out;
  char fname[255];
  double vr=0;
  sprintf(fname,"data/snapshot/PiRprofile_%.2f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    globali=geti(e[Middle][s]);
    globalx=getx(globali,e[Middle][s]);
    out << (s-Middle)/fmtoGeV*AT << "\t";
    out << -pi(1,1,Middle,s)/(4./3.*e[Middle][s])<< endl;
  }
  out.close();
}

void snappibulkprofile(double time)
{
  fstream out;
  char fname[255];
  double vr=0;
  sprintf(fname,"data/snapshot/PiBulkprofile_%.2f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    globali=geti(e[Middle][s]);
    globalx=getx(globali,e[Middle][s]);
    out << (s-Middle)/fmtoGeV*AT << "\t";
    out << pib[Middle][s]/pow(AT,4)<< endl;
  }
  out.close();
}


void snapuphiprofile(double time)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/uphiprofile_%.2f.dat",time/fmtoGeV*AT);
  out.open(fname, ios::out);
  for (int s=Middle;s<=NUMT;s++)
  {
    globali=geti(e[Middle][s]);
    globalx=getx(globali,e[Middle][s]);
    out << (s-Middle)/fmtoGeV*AT << "\t";
    out << u[0][Middle][s]/ut(Middle,s) << endl;
  }
  out.close();
}


#ifdef BG
void betzgyulassy(double tt,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{
  for (int sy=1;sy<=NUMT;sy++)
    {
      for (int sx=1;sx<=NUMT;sx++)
	{
	  globali=geti(e[sx][sy]);
	  globalx=getx(globali,e[sx][sy]);
	  bgout << t/fmtoGeV*AT;
	  bgout << "\t";
	  bgout << (sx-Middle)/fmtoGeV*AT;
	  bgout << "\t";
	  bgout << (sy-Middle)/fmtoGeV*AT;
	  bgout << "\t";
	  bgout << T(sx,sy,Tacc)/AT;
	  bgout << "\t";
	  bgout << u[0][sx][sy];
	  bgout << "\t";
	  bgout << u[1][sx][sy];
	  bgout << "\t";
	  bgout << e[sx][sy]/(AT*AT*AT*AT);
	  bgout << "\n";
	}      
    }
}
#endif

#ifdef AMNR
void amnr(double tt,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/Tdata_%05.2f.dat",tt/fmtoGeV*AT);
  out.open(fname, ios::out);

  for (int sy=1;sy<=NUMT;sy++)
    {
      for (int sx=1;sx<=NUMT;sx++)
	{
	  globali=geti(e[sx][sy]);
	  globalx=getx(globali,e[sx][sy]);
	  out << T(sx,sy,Tacc)/AT << "\t";
	}
      out << endl;
    }
  out.close();

  fstream fout;
  sprintf(fname,"data/snapshot/FOdata_%05.2f.dat",tt/fmtoGeV*AT);
  fout.open(fname, ios::out);

  for (int sy=1;sy<=NUMT;sy++)
    {
      for (int sx=1;sx<=NUMT;sx++)
	{
	  globali=geti(e[sx][sy]);
	  globalx=getx(globali,e[sx][sy]);
	  fout << T(sx,sy,Tacc)/AT << "\t";
	  fout << u[0][sx][sy] << "\t";
	  fout << u[1][sx][sy] << "\t";
	  fout << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t";
	}
      fout << endl;
    }
  fout.close();
}
#endif


void snapshot(double tt,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{
  double hel2[4];
  double dt[4];
  double tb[4],td[4];
 
  //fordebug(10,10);
  
  snapTcontour(tt,Tacc);
  snapVcontour(tt);
  //snapFOdata(tt,pacc,cs2acc);

  snapTprofile(tt,Tacc);   
  snapEDprofile(tt);
  snapVprofile(tt);
  //snapVxprofile(tt);
  //snapV2profile(tt);
  //snapV2profile(tt);
  //snapuphiprofile(tt);
  snappieeprofile(tt);
  snappirrprofile(tt);
  snappibulkprofile(tt);
}
