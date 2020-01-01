//vapour pressure
double pvap(int sx, int sy, gsl_interp_accel *Tacc)
{
	double temp=0.0; 
	// temp=0.0143941 *T(sx,sy,Tacc) - 0.0018451 ; //no influence, deactivated
	return temp;
	
}


double upastt(int sx,int sy)
{
  double temp=1.;
  temp+=upast[0][sx][sy]*upast[0][sx][sy]; 
  temp+=upast[1][sx][sy]*upast[1][sx][sy]; //these three lines create (upast_t)^2 = 1+(upast_x)^2+(upast_y)^2
  return sqrt(temp);
}


double upastmu(int mu, int sx, int sy)
{
  if(mu==2)
    return upastt(sx,sy);
  else
    return upast[mu][sx][sy];
}


double Ut(int sx,int sy)
{
  double temp=1.;
  temp+=U[0][sx][sy]*U[0][sx][sy]; 
  temp+=U[1][sx][sy]*U[1][sx][sy]; //these three lines create (U_t)^2 = 1+(U_x)^2+(U_y)^2
  return sqrt(temp);
}


double Umu(int mu, int sx, int sy)
{
  if(mu==2)
    return Ut(sx,sy);
  else
    return U[mu][sx][sy];
}


//
// 1st order derivatives
//

double dtu(int mu, int sx, int sy)
{
  double temp=0;
  temp=(Umu(mu,sx,sy)-upastmu(mu,sx,sy))/2.0/EPS;
  return temp;
}

//this provides dx U[i]
double dxU(int i, int sx, int sy)
{
  double temp=0;
  if(sx==1)
    temp=(-3.0*Umu(i,sx,sy) + 4.0*Umu(i,sx+1,sy) - Umu(i,sx+2,sy))/2.0;
  else if(sx==NUMT)
    temp=( 3.0*Umu(i,sx,sy) - 4.0*Umu(i,sx-1,sy) + Umu(i,sx-2,sy))/2.0;
  else
    temp=(Umu(i,sx+1,sy) - Umu(i,sx-1,sy))/2.0;
  return temp;
}

//this provides dy U[i]
double dyU(int mu,int sx,int sy)
{
  double temp=0;
  if(sy==1)
    temp=(-3.0*Umu(mu,sx,sy) + 4.0*Umu(mu,sx,sy+1) - Umu(mu,sx,sy+2))/2.0;
  //    temp=(-Umu(mu,sx,sy) + Umu(mu,sx,sy+1); //primitive derivative
  else if(sy==NUMT)
    temp=( 3.0*Umu(mu,sx,sy) - 4.0*Umu(mu,sx,sy-1) + Umu(mu,sx,sy-2))/2.0;
  //  temp= Umu(mu,sx,sy) - Umu(mu,sx,sy-1); //primitive derivative
  else
    temp=(Umu(mu,sx,sy+1) - Umu(mu,sx,sy-1))/2.0;
  return temp;
}


//this provides dx upast[i]


double dxupast(int i, int sx, int sy)
{
  double temp=0;
  if(sx==1)
    temp=(-3.0*upastmu(i,sx,sy) + 4.0*upastmu(i,sx+1,sy) - upastmu(i,sx+2,sy))/2.0;
  else if(sx==NUMT)
    temp=( 3.0*upastmu(i,sx,sy) - 4.0*upastmu(i,sx-1,sy) + upastmu(i,sx-2,sy))/2.0;
  else
    temp=(upastmu(i,sx+1,sy) - upastmu(i,sx-1,sy))/2.0;
  return temp;
}

//this provides dy upast[i]
double dyupast(int mu,int sx,int sy)
{
  double temp=0;
  if(sy==1)
    temp=(-3.0*upastmu(mu,sx,sy) + 4.0*upastmu(mu,sx,sy+1) - upastmu(mu,sx,sy+2))/2.0;
  else if(sy==NUMT)
    temp=( 3.0*upastmu(mu,sx,sy) - 4.0*upastmu(mu,sx,sy-1) + upastmu(mu,sx,sy-2))/2.0;
  else
    temp=(upastmu(mu,sx,sy+1) - upastmu(mu,sx,sy-1))/2.0;
  return temp;
}




//
//begin 2nd deriv.
//


double dx2u(int mu, int sx, int sy)
{
  double temp=0;
  if(sx==1)
    temp=umu(mu,sx,sy) - 2.0*umu(mu,sx+1,sy) + umu(mu,sx+2,sy);  
  else if(sx==NUMT)
    temp=umu(mu,sx,sy) - 2.0*umu(mu,sx-1,sy) + umu(mu,sx-2,sy);
  else
    temp=umu(mu,sx+1,sy) - 2.0*umu(mu,sx,sy) + umu(mu,sx-1,sy);
  return temp;
}

double dy2u(int mu, int sx, int sy)
{
  double temp=0;
  if(sy==1)
    temp=umu(mu,sx,sy) - 2.0*umu(mu,sx,sy+1) + umu(mu,sx,sy+2);  
  else if(sy==NUMT)
    temp=umu(mu,sx,sy) - 2.0*umu(mu,sx,sy-1) + umu(mu,sx,sy-2);
  else
    temp=umu(mu,sx,sy+1) - 2.0*umu(mu,sx,sy) + umu(mu,sx,sy-1);
  return temp;
}

double dt2u(int mu, int sx, int sy)
{
  double temp=0;
  temp=(upastmu(mu,sx,sy) - 2.0*umu(mu,sx,sy) + Umu(mu,sx,sy))/EPS/EPS;
  return temp;
}


//this provides dydx u[i]
double dydxu(int i, int sx, int sy)
{
  double temp=0;
  if(sy==1)
    temp=(-3.0*dxu(i,sx,sy) + 4.0*dxu(i,sx,sy+1) - dxu(i,sx,sy+2))/2.0;
  else if(sy==NUMT)
    temp=( 3.0*dxu(i,sx,sy) - 4.0*dxu(i,sx,sy-1) + dxu(i,sx,sy-2))/2.0;
  else
    temp=(dxu(i,sx,sy+1) - dxu(i,sx,sy-1))/2.0;
  return temp;
}

//this provides dydx u[i]
double dxdyu(int i, int sx, int sy)
{
  double temp=0;
  if(sx==1)
    temp=(-3.0*dyu(i,sx,sy) + 4.0*dyu(i,sx+1,sy) - dyu(i,sx+2,sy))/2.0;
  else if(sx==NUMT)
    temp=( 3.0*dyu(i,sx,sy) - 4.0*dyu(i,sx-1,sy) + dyu(i,sx-2,sy))/2.0;
  else
    temp=(dyu(i,sx+1,sy) - dyu(i,sx-1,sy))/2.0;
  return temp;
}

double DmuUmu(int sx, int sy, double time) //summed components
{
  double temp=0;
  temp =dxu(0,sx,sy) + dyu(1,sx,sy) + dtu(2,sx,sy) + umu(2,sx, sy)/time;
  return temp;	
}

// this provides d_mu u^nu; all 1st order partial derivatives; not summed, though!
double dmuUnu(int mu, int nu, int sx, int sy)
{
  double temp=0;
  if(mu==0) //d_x
    temp=dxu(nu,sx,sy);
  else if(mu==1) //d_y
    temp=dyu(nu,sx,sy);
  else if (mu==2) //d_t
    temp=dtu(nu,sx,sy);
  else
    temp=0;
  return temp;
}

//
// Building Xi terms
//


double H60term(int sx, int sy, double time)
{
  double temp=0;

  // arising from christoffelsymbol
  temp= - umu(2,sx,sy)*umu(2,sx,sy)/time 
	+ umu(2,sx,sy)*dmuUnu(2,2,sx,sy)  
	+ umu(0,sx,sy)*dmuUnu(0,2,sx,sy) + umu(1,sx,sy)*dmuUnu(1,2,sx,sy); 
  temp/=time;
  //ux dx2 ux + uy dy2 uy + ut dt2 ut
  temp+= umu(0,sx,sy)*dx2u(0,sx,sy) + umu(1,sx,sy)*dy2u(1,sx,sy) +umu(2,sx,sy)*dt2u(2,sx,sy);
  return temp;
}


double H61term(int sx, int sy)
{
  double temp=0;
  temp+= umu(0,sx,sy) * (dxU(2,sx,sy) - dxupast(2,sx,sy))/2.0/EPS; //ux * dt * dx ut
  temp+= umu(0,sx,sy) * dxdyu(1,sx,sy); //ux * dx * dy uy
  return temp;
}

double H62term(int sx, int sy)
{
  double temp=0;
  temp+= umu(1,sx,sy) * (dyU(2,sx,sy) - dyupast(2,sx,sy))/2.0/EPS; //uy * dt * dy ut
  temp+= umu(1,sx,sy) * dxdyu(0,sx,sy) ; //uy * dy * dx ux
  return temp;
}

double H63term(int sx, int sy)
{
  double temp=0;
  temp+= umu(2,sx,sy) * (dxU(0,sx,sy)-dxupast(0,sx,sy))/2.0/EPS; //ut * dt * dx ux
  temp+= umu(2,sx,sy) * (dyU(1,sx,sy)-dyupast(1,sx,sy))/2.0/EPS; //ut * dt * dy uy
  return temp;
}


double H6term(int sx, int sy, double time)
{
  double temp=0;
  temp+=H60term(sx,sy,time);
  temp+=H61term(sx,sy);
  temp+=H62term(sx,sy);
  temp+=H63term(sx,sy);
  temp*=(4.0-log(4.0));
  return temp;
}


double H7term(int sx, int sy)
{
	double temp=0.;
	for(int mu=0;mu<4;mu++)
	  {
	    for(int nu=0;nu<4;nu++)
	      {
		for(int alpha=0;alpha<4;alpha++)
		  {
		    for(int beta=0;beta<4;beta++)
		      {
			temp=pishell(nu,mu,sx,sy)*gdown(mu,alpha)*pishell(alpha,beta, sx,sy)*gdown(beta,nu);
		      }
		
		  }
	      }
	  }
	return temp;	
}

double H8term(int sx, int sy, double time,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc)
{
  double temp=0.;
  temp = (4.0-log(4.0))*cs2(sx,sy,pacc,cs2acc)*DmuUmu(sx,sy,time)*DmuUmu(sx,sy,time);
  return temp;
}

//building coefficients using relation between eta & zeta in arXiv:0906.4787v2 
// zeta/s stays linear
// Xi term, 2nd order gradient
double Xi(int sx, int sy, double time,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc)
{
	double temp=0.0;
	temp+=H6term(sx,sy,time);
	temp+=H7term(sx,sy);
	temp+=H8term(sx,sy,time, pacc, cs2acc);
	return temp;
}


double bulkzeta1st(int sx, int sy, double time, gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{
  double temp=0;
  temp = (eos(e[sx][sy],pacc,cs2acc)-pvap(sx, sy,Tacc)) * T(sx, sy,Tacc);
  temp/= DmuUmu(sx, sy, time) ;
  temp/= ( eos(e[sx][sy],pacc,cs2acc) + e[sx][sy] );
	return temp;
}


double bulkzeta2nd(int sx, int sy, double time, gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{
  double temp=0;
  temp = (eos(e[sx][sy],pacc,cs2acc)-pvap(sx, sy,Tacc)) * T(sx, sy,Tacc);
  temp/= DmuUmu(sx, sy, time) - etaos(T(sx,sy,Tacc))*Xi(sx,sy,time,pacc,cs2acc)/T(sx,sy,Tacc);
  temp/= eos(e[sx][sy],pacc,cs2acc) + e[sx][sy] ;
  return temp;
}
  

void bulkvisc(double tt,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{
  string tab="\t";
  int sy=Middle;  
  for (int sx=Middle;sx<NUMT+1;sx++) //safe of boundary effects
    {
      //Tbulkvisc << (sx-Middle)/5.06842*AT << tab; //1
      Tzetaos << T(sx,sy,Tacc)/AT << tab; //2
      Tzetaos << bulkzeta1st(sx,sy,tt,pacc,cs2acc,Tacc) << tab;//3
      Tzetaos << bulkzeta2nd(sx,sy,tt,pacc,cs2acc,Tacc) << tab;//4
      Tzetaos << endl;
    }
}

//routine to calculate the effective pressure (averaged components of the EMT)
double peff(int sx, int sy, double tt, gsl_interp_accel *pacc, gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc, bool sndorder)
{
  double temp=0;
  double zetaos=0;
  zetaos=zetaoeta(T(sx,sy,Tacc)) * etaos(T(sx,sy,Tacc));
  //Switch for second order terns
  if (sndorder)
    {
      temp+= Xi(sx,sy,tt,pacc,cs2acc)/AT/AT;
      temp*= zetaos *  etaos(T(sx,sy,Tacc)) / T(sx,sy,Tacc) * AT;
    }
  temp+= eos(e[sx][sy],pacc,cs2acc)/pow(AT,4);
  temp-= zetaos * DmuUmu(sx,sy,tt) /AT;
  return temp;
}


void peffout(double tt, gsl_interp_accel *pacc, gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{
  string tab="\t";
  fstream out;
  char fname[255];
  sprintf(fname,"data/snapshot/Peff_%.2f.dat",tt/fmtoGeV*AT);
  out.open(fname, ios::out);
  /*
    if ( ZETANAME != zetaoverS-zero.dat )
    cout << Non-zero bulkviscosity << endl;
  */
  out << ZETANAME << endl;
  for (int sy=Middle;sy<NUMT+1;sy++)  
    {
      for (int sx=Middle;sx<NUMT+1;sx++) //safe of boundary effects
	{
	  //	  if( peff(sx,sy,tt,pacc,cs2acc,Tacc,false) <= 0 || peff(sx,sy,tt,pacc,cs2acc,Tacc,true) <= 0 ) //this is cavitation
	  if(true)
	    {
	      globali=geti(e[sx][sy]);
	      globalx=getx(globali,e[sx][sy]);
	      out << (sx-Middle)/fmtoGeV*AT << tab;
	      out << (sy-Middle)/fmtoGeV*AT << tab;
	      out << T(sx,sy,Tacc)/AT << tab; //3
	      out << peff(sx,sy,tt,pacc,cs2acc,Tacc,false) << tab; //4
	      out << peff(sx,sy,tt,pacc,cs2acc,Tacc,true) << tab; //5
	      out << eos(e[sx][sy],pacc,cs2acc)/pow(AT,4) - pib[sx][sy]/pow(AT,4) << tab; //6
	      out << eos(e[sx][sy],pacc,cs2acc)/pow(AT,4) << endl; //7
	    }
	}
    }
  out.close();
}

void snapshotBulkvisc(double tt,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{
  bulkvisc(tt, pacc, cs2acc, Tacc); 
  peffout(tt,pacc,cs2acc,Tacc);
}
