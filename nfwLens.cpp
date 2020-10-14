#include "Cosmology.h"
#include "nfwLens.h"
#include <gsl/gsl_sf_bessel.h>

// properties of the linear matter power spectrum.............
std:: string method_Pk = "EisensteinHu";
std:: string method_bias = "Tinker";
const bool do_nonlinear = false;
// const double zs = 1.5;
//............................................................

const size_t nw1 = 2048;
const size_t nw2 = 2048;
const size_t nw3 = 2048;
const size_t nw4 = 2048;

double epsabs = 1e-3; // 1e-5
double epsrel = 1e-3; // 1e-3

double epsabsq = 0; // 0
double epsrelq = 0.5e-2; // 1e-3

nfwLens:: nfwLens(cbl:: cosmology:: Cosmology* _cosmo,
		  const double& _redshift,
		  const double& _m200,
		  const double& _conc, 
		  const double& _trunc_fact_,
		  const double& _mis_frac,
		  const double& _mis_scale,
		  const bool& _use_2halo)
  : cosmo(_cosmo), redshift(_redshift), m200(_m200), conc(_conc), trunc_fact(_trunc_fact_),
    mis_frac(_mis_frac), mis_scale(_mis_scale), use_2halo(_use_2halo){
  
  // calcualte r200 for given mass and cosmology
  double H = cosmo->HH(redshift)/cosmo->HH(0.)*100./3.0857e+19; // in sec^-1
  double rho_crit = 3*H*H/8/M_PI/6.6732e-8/1.98892e+33*3.0857e+24*3.0857e+24*3.0857e+24; // in h^2*M_sun/Mpc^3
  r200 = pow( 3*m200/(4*M_PI*200.*rho_crit), 1./3. );
  
  rscale = r200/conc;
  trunc = trunc_fact*r200;
  tau = trunc/rscale;
  
  // here we calculate rhos, such that the integral inside r200 is M200
  calc_rhos();
  // allocate gsl workspaces to save time
  w1  = gsl_integration_workspace_alloc (nw1);  
  w2  = gsl_integration_workspace_alloc (nw2);  
  w3  = gsl_integration_workspace_alloc (nw3);
  w4  = gsl_integration_workspace_alloc (nw4);  

  /** ... not necessary for now, but sigma_crit can be computed in this way!
      double Ds = cosmo->D_C(zs);
      double Dl  = cosmo->D_C(redshift);
      double Dls = (Ds - Dl)/(1.+zs); // angularDiameterDistance
      Ds = Ds/(1+zs); // angularDiameterDistance
      Dl = Dl/(1+redshift); // angularDiameterDistance all in [Mpc/h]
      double newton = 4.3011e-9;
      double light = 2.99792e5;
      sigma_crit = light*light/(4.0*M_PI*newton)*Ds/Dl/Dls;  
  */
}

struct nfw_struct
{
  nfwLens* nfw;
};

struct nfw_struct_and_rad
{
  nfwLens* nfw;
  double rad;
};

struct data_for_meantau_factor
{
  double tau;
  double mis_scale;
  gsl_integration_workspace * w;
};

struct data_for_tau_factor
{
  double a;
  double b;
  double inv_mis_scale;
};

double integrand_tau_factor1 (double s, void * params)
{
  double a = ((struct data_for_tau_factor *)params)->a;
  double b = ((struct data_for_tau_factor *)params)->b;
  double inv_mis_scale = ((struct data_for_tau_factor *)params)->inv_mis_scale;

  double tau_factor = (exp(-0.5*(s*s+a)*inv_mis_scale))/sqrt(b-a-s*s);

  return tau_factor;
}

double integrand_tau_factor2 (double s, void * params)
{
  double a = ((struct data_for_tau_factor *)params)->a;
  double b = ((struct data_for_tau_factor *)params)->b;
  double inv_mis_scale = ((struct data_for_tau_factor *)params)->inv_mis_scale;
  
  double tau_factor = (exp(-0.5*(b-s*s)*inv_mis_scale))/sqrt(b-a-s*s);

  return tau_factor;
}

double nfwLens::deltasigma_1h(const double& rad){
  return meansigma_cen(rad) - sigma_cen(rad);
}

double nfwLens:: deltasigma_1h_mis(const double& rad){
  double Ds1h;
  if (mis_frac <= 1.e-5 || mis_scale <= 1.e-5) 
    Ds1h = deltasigma_1h(rad);    
  else{
    double a = deltasigma_1h(rad);
    double b = deltasigma_off(rad);
    // std:: cout << "  " << std:: endl;
    // std:: cout << a << "  " << b << std:: endl;
    Ds1h =  (1-mis_frac)*a+mis_frac*b;
  }
  return Ds1h;
}

double nfwLens:: deltasigma(const double& rad){
  double Ds1h;
  if (mis_frac <= 1.e-5 || mis_scale <= 1.e-5) 
    Ds1h = deltasigma_1h(rad);    
  else{
    double a = deltasigma_1h(rad);
    double b = deltasigma_off(rad);
    // std:: cout << "  " << std:: endl;
    // std:: cout << a << "  " << b << std:: endl;
    Ds1h =  (1-mis_frac)*a+mis_frac*b;
  }
  if (use_2halo == false)
    return Ds1h;
  else
    {
      double Ds2h = deltasigma_2h(rad);
      return Ds1h+Ds2h;
    }  
}

double nfwLens::deltasigma_off(const double& rad)
{
  return meansigma_off(rad) - sigma_off(rad);
}

double integrand_sigma_off(double tau, void * params)
{
  nfwLens* nfw = ((struct nfw_struct_and_rad *)params)->nfw;
  double rad = ((struct nfw_struct_and_rad *)params)->rad;  
  double mis_scale = nfw->get_mis_scale();
  gsl_integration_workspace * w2 = nfw->w2;

  double sigma_off_int = 0.;

  double a = (tau-rad)*(tau-rad);
  double b = (tau+rad)*(tau+rad);
  double c = (a+b)*0.5;
  double inv_mis_scale = 1./mis_scale/mis_scale;

  double result, error;

  gsl_function F;
  F.function = &integrand_tau_factor1;

  data_for_tau_factor params_for_tau_factor = {a,b,inv_mis_scale};
  F.params = &params_for_tau_factor;

  double int_left = 0.;
  double int_right = std::min(sqrt(c-a),sqrt(25*mis_scale*mis_scale-a));

  gsl_integration_qag(&F,int_left,int_right,epsabs,epsrel,nw2,1,w2,&result,&error);
  
  sigma_off_int += result*rad*2.;
  
  if (c < 25*mis_scale*mis_scale)
    {
      F.function = &integrand_tau_factor2;

      int_left = std::max(0.,sqrt(b-25*mis_scale*mis_scale));
      int_right = sqrt(b-c);
      gsl_integration_qag(&F,int_left,int_right,epsabs,epsrel,nw2,1,w2,&result,&error);
      
      sigma_off_int += result*rad*2.;
    }
  
  sigma_off_int *= tau*nfw->sigma_cen(tau);

  return sigma_off_int;
}

double nfwLens:: sigma_off(const double& rad)
{
  double result, error;

  gsl_function F;
  F.function = &integrand_sigma_off;

  nfw_struct_and_rad nfw_and_rad = {this,rad};
  F.params = &nfw_and_rad;
  gsl_integration_qags(&F,std::max(0.,rad-5*mis_scale),rad+5.*mis_scale,epsabs,epsrel,nw1,w1,&result,&error);
  
  double sigma_off = result/rad/M_PI/mis_scale/mis_scale;

  return sigma_off;
}

double integrand_meantau_factor(double rad, void * params)
{
  double tau = ((struct data_for_meantau_factor *)params)->tau;
  double mis_scale = ((struct data_for_meantau_factor *)params)->mis_scale;
  gsl_integration_workspace * w3 = ((struct data_for_meantau_factor *)params)->w;
  double inv_mis_scale = 1./mis_scale/mis_scale;

  double meantau_factor = 0.;

  double a = (tau-rad)*(tau-rad);
  double b = (tau+rad)*(tau+rad);
  double c = (a+b)*0.5;

  double result, error;
  
  gsl_function F;
  F.function = &integrand_tau_factor1;
  
  data_for_tau_factor params_for_tau_factor = {a,b,inv_mis_scale};
  F.params = &params_for_tau_factor;

  double int_left = 0.;
  double int_right = std::min(sqrt(c-a),sqrt(25*mis_scale*mis_scale-a));

  gsl_integration_qag(&F,int_left,int_right,epsabs,epsrel,nw3,1,w3,&result,&error);
  
  meantau_factor += result*rad*2.;

  if (c < 25*mis_scale*mis_scale)
    {
      F.function = &integrand_tau_factor2;

      int_left = std::max(0.,sqrt(b-25*mis_scale*mis_scale));
      int_right = sqrt(b-c);
      gsl_integration_qag(&F,int_left,int_right,epsabs,epsrel,nw3,1,w3,&result,&error);
      
      meantau_factor += result*rad*2.;
    }
  return meantau_factor;
}

double integrand_meansigma_off(double tau, void * params)
{
  double rad = ((struct nfw_struct_and_rad *)params)->rad;
  nfwLens* nfw = ((struct nfw_struct_and_rad *)params)->nfw;
  double mis_scale = nfw->get_mis_scale();
  gsl_integration_workspace * w2 = nfw->w2;

  double result, error;

  gsl_function F;
  F.function = &integrand_meantau_factor;

  data_for_meantau_factor params_for_meantau_factor = {tau,mis_scale,nfw->w3};
  F.params = &params_for_meantau_factor;
  gsl_integration_qags(&F,std::max(0.,tau-5*mis_scale),std::min(rad,tau+5*mis_scale),epsabs,epsrel,nw2,w2,&result,&error);

  double sigma_off_int = result*tau*nfw->sigma_cen(tau);

  return sigma_off_int;  
}

double nfwLens::meansigma_off(const double& rad)
{
  double result, error;

  gsl_function F;
  F.function = &integrand_meansigma_off;

  nfw_struct_and_rad nfw_and_rad = {this,rad};
  F.params = &nfw_and_rad;
  gsl_integration_qags(&F,0.,rad+5.*mis_scale,epsabs,epsrel,nw1,w1,&result,&error);
  
  double meansigma_off = 2.*result/M_PI/mis_scale/mis_scale/rad/rad;

  return meansigma_off;
}

double integrand_2h(double l, void * params){
  nfwLens* nfw = ((struct nfw_struct_and_rad *)params)->nfw;
  double rad = ((struct nfw_struct_and_rad *)params)->rad;  
  double Dl = nfw->get_cosmo()->D_C(nfw->get_redshift()); // comoving
  double theta = rad/Dl;
  double x = theta*l;
  double J2 = gsl_sf_bessel_Jn(2,x);
  double k1 = l/Dl;  
  double Pklinz = nfw->get_cosmo()->Pk(k1,method_Pk,do_nonlinear,nfw->get_redshift());
  return l*J2*Pklinz;
}

double integrand_2hl(double l, void * params){
  l = pow(10.,l);
  nfwLens* nfw = ((struct nfw_struct_and_rad *)params)->nfw;
  double rad = ((struct nfw_struct_and_rad *)params)->rad;  
  double Dl = nfw->get_cosmo()->D_C(nfw->get_redshift()); // comoving
  double theta = rad/Dl;
  double x = theta*l;
  double J2 = gsl_sf_bessel_Jn(2,x);
  double k1 = l/Dl;  
  double Pklinz = nfw->get_cosmo()->Pk(k1,method_Pk,do_nonlinear,nfw->get_redshift());
  return l*J2*Pklinz*l*log(10.);
}

double nfwLens::deltasigma_2h(const double& rad){
  double rhobz = cosmo->rho_m(redshift);
  double Dl = cosmo->D_C(redshift);
  double halo_bias = cosmo->bias_halo(m200,redshift,method_bias,method_Pk);
  double factor = halo_bias*rhobz/(pow((1+redshift),1.0)*pow(Dl,2.0)*2.0*M_PI);
  gsl_function func;
  func.function = &integrand_2hl;
  nfw_struct_and_rad cosmo_param_and_rad = {this,rad};
  func.params = &cosmo_param_and_rad;
  double result, error;  
  // gsl_integration_qagiu(&func,0.0,epsabsq,epsrelq,nw4,w4,&result,&error);
  gsl_integration_qag(&func,-5.0,5.0,epsabs,epsrel,nw4,6,w4,&result,&error);  
  return result*factor/1e12; // convert in pc
}

double integrand_func_sigma_cen (double rad, void * params)
{
  nfwLens* nfw = ((struct nfw_struct *)params)->nfw;
  double sigma = nfw->sigma_cen(rad)*2.*rad;
  return sigma;
}

double nfwLens::meansigma_cen(const double& rad){
  double result, error;

  gsl_function F;
  F.function = &integrand_func_sigma_cen;

  nfw_struct params_for_sigma_cen = {this};
  F.params = &params_for_sigma_cen;
  gsl_integration_qags(&F,0.,rad,epsabs,epsrel,nw1,w1,&result,&error);
  
  double meansigma_cen = result/rad/rad;

  return meansigma_cen;
}

double nfwLens::sigma_cen(const double& rad)
{
  double x = rad/rscale;
  double tausq = tau*tau;
  double tausq_xsq = tausq+x*x; 
  double tau4th = tausq*tausq;
  double prefact = rhos*rscale*tau4th/(tausq+1.)/(tausq+1.)/(tausq+1.);
  double a = 2*(tausq+1.)/(x*x-1.)*(1-F(x));
  double b = 8*F(x);
  double c = (tau4th-1.)/tausq/tausq_xsq;
  double d = -M_PI*(4*tausq_xsq+tausq+1.)/pow(tausq_xsq,1.5);
  double e = (tausq*(tau4th-1.)+(tausq_xsq)*(3*tau4th-6*tausq-1))*L(x)/(tau*tau*tau)/pow(tausq_xsq,1.5);
  return prefact*(a+b+c+d+e)/1.e+12; // in (M_sun/h)/((pc/h)^2) = h*M_sun/pc^2
}

double nfwLens::density_nfw(const double& rad)
{
  double x = rad/rscale;
  return rhos/x/(1+x)/(1+x); // in (M_sun/h)/((Mpc/h)^3)  
}

double nfwLens:: density(const double& rad)
{
  double x = rad/rscale;
  return rhos/x/(1+x)/(1+x)*(tau*tau/(tau*tau+x*x))*(tau*tau/(tau*tau+x*x)); // in (M_sun/h)/((Mpc/h)^3)
}

double integrand_mass(double rad, void * params)
{
  nfwLens* nfw = ((struct nfw_struct *)params)->nfw;
  return 4.*M_PI*rad*rad*nfw->density(rad);
}

void nfwLens::calc_rhos(){
  rhos = 1.;
  double result, error;
  size_t nw = 2048;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (nw);
  
  gsl_function F;
  F.function = integrand_mass;
  nfw_struct params_for_integrate_mass = {this};

  F.params = &params_for_integrate_mass;

  gsl_integration_qags(&F,0.,r200,epsabs,epsrel,nw,w,&result,&error);

  double mass_int = result;
  gsl_integration_workspace_free (w);

  rhos = m200/mass_int;
}

double nfwLens::F(const double & x)
{
  // both operations must be generalised to complex numbers when x<1
  // in that case, both num and den are imaginary, and their ratio is real again
  std::complex<double> num = std::acos(std::complex<double>(1./x,0.));
  std::complex<double> den = std::sqrt(std::complex<double>(x*x-1.,0.));  
  std::complex<double> res = num/den;
  return std::abs(real(res));
}

double nfwLens::L(const double & x)
{
  return log(x/(sqrt(x*x+tau*tau)+tau));
}

double nfwLens::get_M200m()
{
  double H = cosmo->HH(redshift)/cosmo->HH(0.)*100./3.0857e+19; // in sec^-1                                            
  double rho_crit = 3*H*H/8/M_PI/6.6732e-8/1.98892e+33*3.0857e+24*3.0857e+24*3.0857e+24; // in h^2*M_sun/Mpc^3
  double mass_int = m200;
  double rad = r200;  
  double enclosed_density = mass_int*3/(4*M_PI*rad*rad*rad);
  // need to be checked!
  double mean_density = rho_crit*cosmo->OmegaM(redshift);
  mean_density =  cosmo->rho_m(redshift); // with the CBL
  while (enclosed_density > 200*mean_density)
    {
      mass_int += 4*M_PI*rad*rad*density(rad);

      rad += 0.01*r200;

      enclosed_density = mass_int*3/(4*M_PI*rad*rad*rad);
    }
  return rad;
}

