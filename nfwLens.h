#ifndef LIBCLUSTERWL_NFWLENS_H
#define LIBCLUSTERWL_NFWLENS_H

#include <fstream>
#include <iostream>
#include <vector>
#include "Cosmology.h"
#include <gsl/gsl_integration.h>

// Class defining a truncated NFW lens
// (Baltz, Marshall & Oguri, 2009)
// + mis-centering term (Sereno, Veropalumbo et al. 2016)
// + 2-halo term (from Mauro S. tables)
class nfwLens
{
 private:
  cbl:: cosmology:: Cosmology* cosmo;
  double sigma_crit;
  double redshift;
  double m200;
  double r200;
  double rscale;
  double conc;

  double rhos; // scale density (=M0/(4 pi rs^3))

  double trunc_fact;
  
  double trunc; // truncation radius
  double tau; // trunc/rscale

  double mis_frac; // fraction of miscentered halos in the stacking
  double mis_scale; // scale length of off-sets (Mpc/h)

  double sigma_fac; // c^2/(4piG) in Msun/pc  
  double clu_dist;

  bool use_2halo; // true if 2-halo term is used
  double bias;

  // std::vector<double> gamma_arr;
  // std::string bias_file = "/home/lorenzo/KIDS/cluster_wl_lib/data_halo_bias_TI+10_LCDM.dat";
  // std::string gamma_file = "/home/lorenzo/KIDS/cluster_wl_lib/data_2-halo_gamma_t_LCDM.dat";

  void calc_rhos();

  double L(const double& x);
  double F(const double& x);

 public:  
  nfwLens(cbl:: cosmology:: Cosmology* _cosmo, const double& _redshift, const double& _m200,
	  const double& _conc, const double& trunc_fact, const double& _mis_frac = 0.,
	  const double& _mis_scale = 0., const bool& _use_2halo = false);
  double deltasigma(const double& rad); 
  double deltasigma_1h(const double& rad);
  double deltasigma_1h_mis(const double& rad);
  double deltasigma_2h(const double& rad);   
  double deltasigma_off(const double& rad); 
  double meansigma_cen(const double& rad); 
  double meansigma_off(const double& rad); 
  double sigma_cen(const double& rad); 
  double sigma_off(const double& rad); 
  double sigma_tot(const double& rad); 
  // double 2halo(const double& rad);

  double density(const double& rad); // 3d density
  double density_nfw(const double& rad); // 3d density, ignoring truncation

  double get_redshift() {return redshift;};
  double get_r200() {return r200;};  

  double get_M200m();

  double get_mis_scale(){return mis_scale;};

  cbl:: cosmology:: Cosmology* get_cosmo(){return cosmo;};

  gsl_integration_workspace * w1, * w2, * w3, *w4;
  
};
    
#endif


