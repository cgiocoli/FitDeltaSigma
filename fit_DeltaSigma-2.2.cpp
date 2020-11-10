// ========================================================================================
// Perform a Bayesian fit to a set of data points with a generic model
// ========================================================================================

#include "Cosmology.h"
#include "Data1D.h"
#include "Posterior.h"
#include "nfwLens.h"
#include "TaperedCovarianceMatrix.h"
#include <numeric>

using namespace std;
double epsilon = 8e-3;
std:: string fepsilon = "change_epsilon.txt";
double average_redshift;

double prec = 1e-10;
std:: string fprec = "change_prec.txt";

// int max_iter=10000000;
// double tol = 1.e-1;

// 3.0 reference by Bellagamba+19
double truncation_radius = 3.0;

// halo bias model
std:: string halo_bias_model;

// min gal per bin - project like 300th
double min_ngal_in_bin=10;
std:: string fngal_in_bin = "ngals_in_bin.txt";

// these two variables contain the name of the CosmoBolognaLib
// directory and the name of the current directory (useful when
// launching the code on remote systems)
string cbl::par::DirCosmo = DIRCOSMO, cbl::par::DirLoc = DIRL;
std:: vector<double> mean_redshift,d_c_referenceLCDM;
std:: vector<double> scrit,z_sources;
std:: vector<double> one_over_scrit;

// this function marginalise over 2 parameters:
// 1 - Mass (M_200c)
// 2 - Concentration (c_200c)
vector<double> model_function2(const vector<double> x, const shared_ptr<void> modelInput, std::vector<double> &parameter){
  // the object Cosmology, used to create the NfwLens
  cbl::cosmology::Cosmology cosm = *static_pointer_cast<cbl::cosmology::Cosmology>(modelInput);
  
  // initialize an NfwLens, parameters
  // 1: cosmology
  // 2: redshift
  // 3: m200
  // 4: concentration
  // 5: truncation factor in unit of r200
  // 6: miscentering factor f_off
  // 7: miscentering scale  \sigma_off
  // 8: bool use or not the two halo term
  // parameter 1
  double c200 = parameter[1];

  // parameter 0
  double m200 = pow(10., parameter[0]);
  double f_off = 0;
  double sigma_f_off = 1e-6;
  
  nfwLens lens(&cosm, average_redshift, m200, c200, truncation_radius, f_off, sigma_f_off, true);
  
  vector<double> model(x.size(), 0.);
  for (size_t i = 0; i < x.size(); ++i)
    // model[i] = lens.deltasigma(x[i]);
    model[i] = lens.deltasigma_1h(x[i]);  
  return model;
}
// =====================================================================

// this function marginalise over 4 parameters:
// 1 - Mass (M_200c)
// 2 - Concentration (c_200c)
// 3 - Offset Parameter (f_off)
// 4 - Dispersion of the Offest (sigma_off)
vector<double> model_function4 (const vector<double> x, const shared_ptr<void> modelInput, std::vector<double> &parameter){
  // the object Cosmology, used to create the NfwLens
  cbl::cosmology::Cosmology cosm = *static_pointer_cast<cbl::cosmology::Cosmology>(modelInput);
  
  // initialize an NfwLens, parameters
  // 1: cosmology
  // 2: redshift
  // 3: m200
  // 4: concentration
  // 5: truncation factor in unit of r200
  // 6: miscentering factor f_off
  // 7: miscentering scale  \sigma_off
  // 8: bool use or not the two halo term
  // parameter 1
  double c200 = parameter[1];

  // parameter 0
  double m200 = pow(10.,parameter[0]);
  // parameter 2
  double f_off = parameter[2];
  // parameter 3
  double sigma_f_off = parameter[3];
  
  nfwLens lens(&cosm, average_redshift, m200, c200, truncation_radius, f_off, sigma_f_off, true, halo_bias_model);
  
  vector<double> model(x.size(), 0.);
  for (size_t i=0; i<x.size(); ++i)
    model[i] = lens.deltasigma(x[i]);
  return model;
}
// =====================================================================

// this function marginalise over 5 parameters:
// 1 - Mass (M_200c)
// 2 - Concentration (c_200c)
// 3 - Offset Parameter (f_off)
// 4 - Dispersion of the Offest (sigma_off)
// 5 - \Omega_m, \sigma_8 is a derived quantity
vector<double> model_function5 (const vector<double> x, const shared_ptr<void> modelInput, std::vector<double> &parameter){
  // the object Cosmology, used to create the NfwLens
  cbl::cosmology::Cosmology cosm = *static_pointer_cast<cbl::cosmology::Cosmology>(modelInput);

  // initialize an NfwLens, parameters
  // 1: cosmology
  // 2: redshift
  // 3: m200
  // 4: concentration
  // 5: truncation factor in unit of r200
  // 6: miscentering factor f_off
  // 7: miscentering scale  \sigma_off
  // 8: bool use or not the two halo term
  // parameter 1
  double c200 = parameter[1]; // as fixed by Bellagamba et al. 2019 to 4 --- maybe we can add a c-M relation

  // parameter 0
  double m200 = pow(10.,parameter[0]);
  // parameter 2
  double f_off = parameter[2];
  // parameter 3
  double sigma_f_off = parameter[3];
  // parameter 4  
  double O_matter = parameter[4];  
  cosm.set_Omega(O_matter);    
  nfwLens lens(&cosm, average_redshift, m200, c200, truncation_radius, f_off, sigma_f_off, true, halo_bias_model);

  vector<double> model(x.size(), 0.);

  std:: vector<double> d_c_cl,d_c_sources,dl,ds,dls;
  // rescale the radius to the new cosmological model
  for(long unsigned int i=0;i<x.size();i++){
    d_c_cl.push_back(cosm.D_C(mean_redshift[i]));
    d_c_sources.push_back(cosm.D_C(z_sources[i]));
  }
  for(long unsigned int i=0;i<x.size();i++){
    dl.push_back(d_c_cl[i]/(1.0+mean_redshift[i]));
    ds.push_back(d_c_sources[i]/(1.0+z_sources[i]));
    dls.push_back((d_c_sources[i]-d_c_cl[i])/(1.0+z_sources[i]));
  }
  for (size_t i=0; i<x.size(); ++i)
    // here we do change it for the cosmological model .... . ... !!!
    model[i] = lens.deltasigma(x[i]/d_c_referenceLCDM[i]*cosm.D_C(mean_redshift[i]))/scrit[i]*(ds[i]/(dl[i]*dls[i]));
  return model;
}

// this function marginalise over 6 parameters:
// 1 - Mass (M_200c)
// 2 - Concentration (c_200c)
// 3 - Offset Parameter (f_off)
// 4 - Dispersion of the Offest (sigma_off)
// 5 - \Omega_m
// 6 - \sigma_8
vector<double> model_function6 (const vector<double> x, const shared_ptr<void> modelInput, std::vector<double> &parameter){
  // the object Cosmology, used to create the NfwLens
  cbl::cosmology::Cosmology cosm = *static_pointer_cast<cbl::cosmology::Cosmology>(modelInput);

  // initialize an NfwLens, parameters
  // 1: cosmology
  // 2: redshift
  // 3: m200
  // 4: concentration
  // 5: truncation factor in unit of r200
  // 6: miscentering factor f_off
  // 7: miscentering scale  \sigma_off
  // 8: bool use or not the two halo term
  // parameter 1
  double c200 = parameter[1]; // as fixed by Bellagamba et al. 2019 to 4 --- maybe we can add a c-M relation

  // parameter 0
  double m200 = pow(10.,parameter[0]);
  // parameter 2
  double f_off = parameter[2];
  // parameter 3
  double sigma_f_off = parameter[3];
  // parameter 4
  double O_matter = parameter[4];
  // parameter 5
  double s_8 = parameter[5];
  cosm.set_Omega(O_matter);
  cosm.set_sigma8(s_8);    
  nfwLens lens(&cosm, average_redshift, m200, c200, truncation_radius, f_off, sigma_f_off, true, halo_bias_model);

  vector<double> model(x.size(), 0.);

  std:: vector<double> d_c_cl,d_c_sources,dl,ds,dls;
  // rescale the radius to the new cosmological model
  for(long unsigned int i=0;i<x.size();i++){
    d_c_cl.push_back(cosm.D_C(mean_redshift[i]));
    d_c_sources.push_back(cosm.D_C(z_sources[i]));
  }
  for(long unsigned int i=0;i<x.size();i++){
    dl.push_back(d_c_cl[i]/(1.0+mean_redshift[i]));
    ds.push_back(d_c_sources[i]/(1.0+z_sources[i]));
    dls.push_back((d_c_sources[i]-d_c_cl[i])/(1.0+z_sources[i]));
  }
  for (size_t i=0; i<x.size(); ++i)
    // here we do change it for the cosmological model .... . ... !!!
    model[i] = lens.deltasigma(x[i]/d_c_referenceLCDM[i]*cosm.D_C(mean_redshift[i]))/scrit[i]*(ds[i]/(dl[i]*dls[i]));
  return model;
}
// =====================================================================

int main () {
  try {
    average_redshift = 0;
    /*************************************************************************************************************************/
    /*************************************************************************************************************************/
    // 300th Cosmology
    double Omega_matter = 0.307;
    double Omega_baryon = 0.048;
    const double Omega_nu = 0.;
    double massless_neutrinos = 2.04;
    const double Omega_radiation = 0.;
    int massive_neutrinos = 1.0;
    double Omega_DE = 0.693;
    double hh = 0.678;
    double scalar_amp = 2.139e-9;
    double scalar_pivot = 0.05;
    double n_spec = 0.96;
    // double tau = 0.066;
    const double wa = 0.;
    const double w0 = -1.;
    //
    /***** COSMOLOGY *****/  /* BEGIN */
    // cbl::cosmology::Cosmology cosmology {cbl::cosmology::CosmologicalModel::_Planck15_};
    cbl::cosmology::Cosmology cosmology {Omega_matter, Omega_baryon,
                                        Omega_nu, massless_neutrinos, massive_neutrinos, Omega_DE, Omega_radiation,
                                        hh, scalar_amp, scalar_pivot, n_spec, w0, wa};
    // as used to comupute the stacked shear profiles
    std:: cout << "s8 EH " << cosmology.sigma8_Pk("EisensteinHu",0.0) << std:: endl;
    // print cosmological parameters
    cosmology.print_parameters();
    // cosmology.set_Omega(0.3);
    // cosmology.set_H0(70.0);
    /***** COSMOLOGY *****/  /* END */
    /*************************************************************************************************************************/
    /*************************************************************************************************************************/
    /**** READING INPUT FILE PROVIDED BY THE USER as:
	  echo INPUT_FILE | ./fit_DeltaSigma-2.2
    ****/

    std:: string filinput;
    std:: cin >> filinput;

    int npars;
    double varA, varB, varC, varD, varE, varF;
    double varAmin, varAmax;
    double varBmin, varBmax;
    double varCmin, varCmax;
    double varDmin, varDmax;
    double varEmin, varEmax;
    double varFmin, varFmax;
    std:: string file_input,file_cov,outputdir;
    std:: string str;
    int iseedlikelihood, inwalkers, ichain_size;
    int Z_bin; // Redshift bin ... first column
    int A_bin; // Amplitude bin ... second column
    double rmin;

    std:: ifstream ifilinput;
    ifilinput.open(filinput.c_str());
    if(ifilinput.is_open()){

      ifilinput >> str;
      ifilinput >> file_input;

      ifilinput >> str;
      ifilinput >> file_cov;

      ifilinput >> str;
      ifilinput >> Z_bin;

      ifilinput >> str;
      ifilinput >> A_bin;

      ifilinput >> str;
      ifilinput >> rmin;

      ifilinput >> str;
      ifilinput >> truncation_radius;

      ifilinput >> str;        
      ifilinput >> halo_bias_model;
      
      // seed and chain parameters
      ifilinput >> str;
      ifilinput >> iseedlikelihood >> inwalkers >> ichain_size;

      ifilinput >> str;
      ifilinput >> outputdir;

      std:: cout << " -------------------------------------------------------------------------------------- " << std:: endl;
      std:: cout << " considering the input file       = " << file_input << std:: endl;
      std:: cout << " considering the covariance file  = " << file_cov << std:: endl;
      std:: cout << " considering the Z bin            = " << Z_bin << std:: endl;
      std:: cout << " considering the A bin            = " << A_bin << std:: endl;
      std:: cout << " considering data points above    = " << rmin << "[Mpc/h]" << std:: endl;
      std:: cout << " considering r_t in unit of R_200 = " << truncation_radius  << std:: endl;
      std:: cout << " halo bias model                  = " << halo_bias_model  << std:: endl;      
      std:: cout << " writing output files in          = " << outputdir << std:: endl;        
           
      ifilinput >> str;
      ifilinput >> npars;          //  number of parameters

      std:: cout << " considering npars = " << npars << std:: endl;

      if(npars < 2){
        std:: cout << " npars smaller than 2 " << std:: endl;
	std:: cout << " I will STOP here!!! " << std:: endl;
	exit(0);
      }
      if(npars > 6){
	std:: cout << " npars larger than 6 " << std:: endl;
	std:: cout << " I will STOP here!!! " << std:: endl;
	exit(0);
      }

      ifilinput >> str;
      ifilinput >> varA;          //  LogM200

      ifilinput >> str;
      ifilinput >> varB;          //  c200

      if(npars>2){

	ifilinput >> str;
	ifilinput >> varC;          //  f_off

	ifilinput >> str;
	ifilinput >> varD;          //  sigma_off
	
	if(npars == 5){
	  ifilinput >> str;
	  ifilinput >> varE;          //  OmegaM
	}
	
	if(npars == 6){
	  ifilinput >> str;
	  ifilinput >> varE;          //  OmegaM
	  
	  ifilinput >> str;
	  ifilinput >> varF;          //  sigma8
	}
      }
      
      ifilinput >> str;
      ifilinput >> varAmin >> varAmax;
      
      ifilinput >> str;
      ifilinput >> varBmin >> varBmax;
      
      if(npars>2){
	ifilinput >> str;
	ifilinput >> varCmin >> varCmax;
	
	ifilinput >> str;
	ifilinput >> varDmin >> varDmax;
	
	if(npars == 5){
	  ifilinput >> str;
	  ifilinput >> varEmin >> varEmax;
	}
	
	if(npars == 6){
	  ifilinput >> str;
	  ifilinput >> varEmin >> varEmax;
	  
	  ifilinput >> str;
	  ifilinput >> varFmin >> varFmax;
	}
      }
      
      std:: cout << " *** Running with:    " << std:: endl;
      std:: cout << "     LogM200 = " << varA << "  ranging from " << varAmin << " to " << varAmax << std:: endl;
      std:: cout << "        c200 = " << varB << "  ranging from " << varBmin << " to " << varBmax << std:: endl;
      if(npars>2){
	std:: cout << "       f_off = " << varC << "  ranging from " << varCmin << " to " << varCmax << std:: endl;
	std:: cout << "   sigma_off = " << varD << "  ranging from " << varDmin << " to " << varDmax << std:: endl;
      }
      std:: cout << " iseedlikelihood = " << iseedlikelihood << std:: endl;
      std:: cout << " inwalkers       = " << inwalkers << std:: endl;
      std:: cout << " ichain_size     = " << ichain_size << std:: endl;

      if(npars == 5){
	std:: cout << "     OmegaM  = " << varE << "  ranging from " << varEmin << " to " << varEmax << std:: endl;
      }

      if(npars == 6){
	std:: cout << "     OmegaM  = " << varE << "  ranging from " << varEmin << " to " << varEmax << std:: endl;
	std:: cout << "     sigma8  = " << varF << "  ranging from " << varFmin << " to " << varFmax << std:: endl;
      }
    }else{
      std:: cout << " input file provided " << filinput << std:: endl;
      std:: cout << " does not exist " << std:: endl;
      std:: cout << " I will STOP here!!! " << std:: endl;
      exit(0);
    }
    /**** DONE READING ****/

    // *** CHECK THE EXISTENCE OF FILES:
    std:: cout << " ...   checking EXTRA PARAMETER FILES ... ... ... " << std:: endl;
    std:: ifstream ifepsilon;
    ifepsilon.open(fepsilon.c_str());
    if(ifepsilon.is_open()){
      std:: cout << " change_epsilon.txt exists I will read it! " << std:: endl;
      ifepsilon >> epsilon;
      std:: cout << " new epsilon value = " << epsilon << std:: endl;
      ifepsilon.close();
    }
    std:: ifstream ifprec;
    ifprec.open(fprec.c_str());
    if(ifprec.is_open()){
      std:: cout << " change_prec.txt exists I will read it! " << std:: endl;
      ifprec >> prec;
      std:: cout << " new prec value = " << prec << std:: endl;
      ifprec.close();
    }
    std:: ifstream ingal_in_bin;
    ingal_in_bin.open(fngal_in_bin.c_str());
    if(ingal_in_bin.is_open()){
      std:: cout << " ngals_in_bin.txt exists I will read it! " << std:: endl;
      ifprec >> min_ngal_in_bin;
      std:: cout << " new min_ngal_in_bin value = " << min_ngal_in_bin << std:: endl;
      ingal_in_bin.close();
    }
    std:: cout << " ...   checking EXTRA PARAMETER FILES ... ... ... " << std:: endl;    
    std:: cout << "    ... ... ... ... CHECK DONE ... ... ... ...    " << std:: endl;
    // *** DONE WITH CHECK *** //
    
    // int seedlikelihood = 12345;
    int seedlikelihood = iseedlikelihood;
    
    std:: ifstream ifile_input;// (file_input);
    ifile_input.open(file_input.c_str());
    
    int i1, i2, i13, i18;
    double a3,a4, a5, a6, a7, a8, a9, a10, a11, a12, a14, a15, a16, a17, a19, a20;
    
    std:: vector<double> redshift_sigma;
    std:: vector<double> mean_amplitude, amplitude_sigma;
    std:: vector<double> radius, radius_sigma;
    std:: vector<double> weight,weighted_amplitude;
    std:: vector<double> DeltaSigma, DeltaSigma_err, DeltaSigma_err_bootstrap,eta;
    std:: vector<int> nclusters_per_bin;
    int nl=0;
    int ndata=0;
    
    if(Z_bin>=0){
      std:: string line;
      int nlines_to_skip = 24;
      for(int i=0;i<nlines_to_skip;i++){
        getline(ifile_input,line);
        nl++;
      }
      
      if(npars == 4){
        while (true){
	  ifile_input >> i1  >> i2  >> a3  >> a4
		      >> a5  >> a6  >> a7  >> a8
		      >> a9  >> a10 >> a11 >> a12
		      >> i13 >> a14 >> a15 >> a16
		      >> a17 >> i18;
	  if(ifile_input.eof()) break;
	  
	  // amplitude bin and redshift bin to be considered
	  
	  if(i1 == Z_bin && i2 == A_bin){
	    std:: cout << i1 << "  " << i2 << "  " << a7 << "   " << a9 << std:: endl;
	    // exclude DeltaSigma<=0
	    //if(a9>0){
	    // exclude values below rmin
	    if(a7>rmin){
	      mean_redshift.push_back(a3);
	      redshift_sigma.push_back(a4);
	      mean_amplitude.push_back(a5);
	      amplitude_sigma.push_back(a6);
	      radius.push_back(a7);
	      radius_sigma.push_back(a8);
	      DeltaSigma.push_back(a9);
	      DeltaSigma_err.push_back(a10);
	      weight.push_back(a14);
	      weighted_amplitude.push_back(a14*a5);
	      DeltaSigma_err_bootstrap.push_back(a15);
	      nclusters_per_bin.push_back(i13);
	    }
	    //}
	    ndata++;
	  }
	  nl++;
        }
      }
      if (npars >= 5){
        while (true){
	  ifile_input >> i1  >> i2  >> a3  >> a4
		      >> a5  >> a6  >> a7  >> a8
		      >> a9  >> a10 >> a11 >> a12
		      >> i13 >> a14 >> a15 >> a16
		      >> a17 >> i18 >> a19 >> a20;
	  if(ifile_input.eof()) break;
	  
      	  // amplitude bin and redshift bin to be considered
	  if(i1 == Z_bin && i2 == A_bin){
	    std:: cout << i1 << "  " << i2 << "  " << a7 << "   " << a9 << std:: endl;
	    // exclude DeltaSigma<=0
	    //if(a9>0){
	    // exclude values below rmin
	    if(a7>rmin){
	      mean_redshift.push_back(a3);
	      redshift_sigma.push_back(a4);
	      mean_amplitude.push_back(a5);
	      amplitude_sigma.push_back(a6);
	      radius.push_back(a7);
	      radius_sigma.push_back(a8);
	      DeltaSigma.push_back(a9);
	      DeltaSigma_err.push_back(a10);
	      weight.push_back(a14);
	      weighted_amplitude.push_back(a14*a5);
	      DeltaSigma_err_bootstrap.push_back(a15);
	      nclusters_per_bin.push_back(i13);
	      z_sources.push_back(a19);
	      eta.push_back(a20);
	    }
	    //}
	    ndata++;
	  }
	  nl++;
        }
      }
      ifile_input.close();
      
      double weighted_amplitude_N=0.;
      double weighted_amplitude_D=0.;
      
      int n = nclusters_per_bin.size();
      std:: cout << " I have found " << nl << " lines in the data file " << std:: endl;
      std:: cout << " I have read and stored " << n << " lines " << "("<<ndata<<")"<< std:: endl;
      
      if ( n != 0) average_redshift = accumulate( mean_redshift.begin(), mean_redshift.end(), 0.0) / n;
      std:: cout << " the average redshift for the clusters is = " << average_redshift << std:: endl;
      double new_average_redshift;
      if(npars >= 5){
        std:: vector<double> redshift_1h;
        std:: string file_redshift = "/home/PERSONALE/carlo.giocoli/Codes/KiDS_ClusterProf/REDSHIFT_CL_1H_and_w2H.txt";
        std:: ifstream ifile_redshift;
        ifile_redshift.open(file_redshift.c_str());
        if(ifile_redshift.is_open()){
	  double z1h, z2h;
      	  while(ifile_redshift >> z1h >> z2h){
	    redshift_1h.push_back(z1h);
	  }
        }else{
	  std:: cout << " redshift file does not exist " << file_redshift << std:: endl;
	  std:: cout << " I will STOP here!!! " << std:: endl;
	  exit(1);
        }
        if(Z_bin == 0 && A_bin == 0){
	  new_average_redshift = redshift_1h[0];
        }
        if(Z_bin == 0 && A_bin == 1){
	  new_average_redshift = redshift_1h[1];
        }
        if(Z_bin == 0 && A_bin == 2){
	  new_average_redshift = redshift_1h[2];
        }
        if(Z_bin == 0 && A_bin == 3){
	  new_average_redshift = redshift_1h[3];
        }
        if(Z_bin == 0 && A_bin == 4){
	  new_average_redshift = redshift_1h[4];
        }
        if(Z_bin == 1 && A_bin == 0){
	  new_average_redshift = redshift_1h[5];
        }
        if(Z_bin == 1 && A_bin == 1){
	  new_average_redshift = redshift_1h[6];
        }
        if(Z_bin == 1 && A_bin == 2){
	  new_average_redshift = redshift_1h[7];
        }
        if(Z_bin == 1 && A_bin == 3){
	  new_average_redshift = redshift_1h[8];
        }
        if(Z_bin == 1 && A_bin == 4){
	  new_average_redshift = redshift_1h[9];
        }
        if(Z_bin == 2 && A_bin == 0){
	  new_average_redshift = redshift_1h[10];
        }
        if(Z_bin == 2 && A_bin == 1){
	  new_average_redshift = redshift_1h[11];
        }
        if(Z_bin == 2 && A_bin == 2){
	  new_average_redshift = redshift_1h[12];
        }
        if(Z_bin == 2 && A_bin == 3){
	  new_average_redshift = redshift_1h[13];
        }
        std:: cout << " the average redshift for the clusters from 1 Halo term is = " << new_average_redshift << std:: endl;
        average_redshift = new_average_redshift;
      }
      
      double average_amplitude;
      if ( n != 0) average_amplitude = accumulate( mean_amplitude.begin(), mean_amplitude.end(), 0.0) / n;
      std:: cout << " the average amplitude for the clusters is = " << average_amplitude << std:: endl;
      
      for(std::vector<double>::iterator it = weighted_amplitude.begin(); it != weighted_amplitude.end(); ++it)
	weighted_amplitude_N += *it;
      
      for(std::vector<double>::iterator it = weight.begin(); it != weight.end(); ++it)
	weighted_amplitude_D += *it;
      
      std:: cout << " the weighted average amplitude for the clusters is = " << weighted_amplitude_N/weighted_amplitude_D << std:: endl;
      std:: vector<double> ds, dls,dl,d_c_sources;
      if(npars >= 5){
        // rescale the radius to the new cosmological model
        for(long unsigned int i=0;i<radius.size();i++){
	  d_c_referenceLCDM.push_back(cosmology.D_C(mean_redshift[i]));
	  d_c_sources.push_back(cosmology.D_C(z_sources[i]));
        }
        for(long unsigned int i=0;i<radius.size();i++){
	  dl.push_back(d_c_referenceLCDM[i]/(1.0+mean_redshift[i]));
	  ds.push_back(d_c_sources[i]/(1.0+z_sources[i]));
	  dls.push_back((d_c_sources[i]-d_c_referenceLCDM[i])/(1.0+z_sources[i]));
        }
        for(long unsigned int i=0;i<radius.size();i++){
	  scrit.push_back(ds[i]/(dl[i]*dls[i]));
	  one_over_scrit.push_back(dl[i]*dls[i]/ds[i]);
	  std:: cout << eta[i]/1e6 << "   " << dl[i]*dls[i]/ds[i] << std:: endl;
        }
        double average_eta = accumulate(eta.begin(), eta.end(), 0.0) / 1e6 / eta.size();
        double average_eta2 = accumulate(one_over_scrit.begin(), one_over_scrit.end(), 0.0) / one_over_scrit.size();
        std:: cout << " average eta (Fabio's code)" << "  " << " average eta " << "   " <<  "1 / averagege scrit " << std:: endl;
        std:: cout << " eventually the mean is = " << average_eta << "   " << average_eta2
		   << "   " << 1.0/(accumulate(scrit.begin(), scrit.end(), 0.0) / scrit.size()) << std:: endl;
      }
      
      // write the values
      std:: cout << "  E(z_average) / E(0.35) = " <<  cosmology.HH(average_redshift) / cosmology.HH(0.35) << std:: endl;
    }else{
      if (average_redshift<1e-8){
        std:: cout << " the average redshift of the clusters is " << average_redshift << std:: endl;
        std:: cout << " please provide a realistic value " << std:: endl;
        std::cin >> average_redshift;
      }
      std::cout << " average_redshift = " <<  average_redshift << std::endl;
      std::cout << "  E(z_average) / E(0.35) = " << cosmology.HH(average_redshift) / cosmology.HH(0.35) << std::endl;
      // read file with only three comumn
      nl = 0;
      double but;
      double ngal_in_bin;
      while (true){
	// ifile_input >> a3  >> a4 >> a5 >> but >> but >> but >> but >> but >> ngal_in_bin;
	ifile_input >> a3  >> but >> but >> but >> but >> but >> a4 >> a5 >> ngal_in_bin;	 
	if(ifile_input.eof()) break;
	if (a3 > rmin){
	  // in this case we do not account for the SL region in simulations
	  //if(a4 > 0){
	  if(ngal_in_bin >= min_ngal_in_bin){
	    DeltaSigma.push_back(a4);
	    DeltaSigma_err.push_back(a5);
	    radius.push_back(a3);
	    // std::cout << a3 << "  " << a4 << "  " << a5 << std::endl;
	    //}
	  }
	  nl++;
	}
      }
      ifile_input.close();
      std:: cout << " I have read " << nl << " lines in the file" << std:: endl;
      ndata = nl;
    }
    // --- set the input/output file/directories ---
    const string dir_output = cbl::par::DirLoc+outputdir+"/";
    
    // --- construct the dataset copying the vectors read
    cbl::data::Data1D data(radius,DeltaSigma,DeltaSigma_err);
    std:: cout << "  " << std:: endl;
    std:: cout << " size of the data " << data.xx().size() << std:: endl;
    int nn = data.xx().size();
    for(int i=0;i<nn;i++){
      std:: cout << data.xx(i) << std:: endl;
    }
    std:: vector<vector<double>> covariance(nn, vector<double>(nn));
    
    int ni = ndata-nn;
    // we read the covariance file if provided and different from NO
    if(file_cov != "NO"){
      std:: ifstream ifile_cov;
      ifile_cov.open(file_cov.c_str());
      if(ifile_cov.is_open()){
	double icov, prec;
	int a,z,aa,bb;
	while(ifile_cov >> z >> a >> aa >> bb >> icov >> prec){
	  if(z == Z_bin && a == A_bin){
	    if(aa>=ni && bb>=ni){
	      // std:: cout << aa << "  " << bb << "  " << aa-ni << "  " << bb-ni << std:: endl;
	      // we are skipping r<rmin in the covariance too
	      covariance[aa-ni][bb-ni]=icov;
	    }
	  }
	}
	data.set_covariance(covariance);
	std:: cout << "  " << std:: endl;
	std:: cout << " covariance read and set " << std:: endl;
	std:: cout << "  " << std:: endl;
      }else{
	std:: cout << file_cov << std:: endl;
	std:: cout << " Covariance file provided does not exist, I will STOP here!!! " << std:: endl;
	exit(1);
      }
      ifile_cov.close();
    }
    
    std:: cout << "  " << std:: endl;
    shared_ptr<cbl::data::Data> ptr_data = make_shared<cbl::data::Data1D>(data);
    
    /*** CONSTANT PARAMETERS ***/
    const int nbins = 50;
    const bool show_mode = false;
    // parameters to sample the posterior
    // const int nwalkers = 32;
    // const int chain_size = 100;
    // long_chain
    // const int nwalkers = 96;
    // const int chain_size = 5000;
    const int nwalkers = inwalkers;
    const int chain_size = ichain_size;
    
    // parameter to show all the MCMC statistics on screen
    const int burn_in = 0;
    const int thin = 1;    
    /**************************/
    
    if(npars == 2){
      // --- set the model to construct the likelihood ---
      // number of model parameters
      const int nparameters = npars;
      
      // names of the model parameters
      const vector<string> parNames = {"LogM200", "c200"};
      
      // starting values ... as READ IN THE INPUT FILE
      double valA = varA, valB = varB;
      
      // vector containing the model parameters
      vector<cbl::statistics::ParameterType> parType(nparameters);
      
      auto ptr_modelInput = make_shared<cbl::cosmology::Cosmology>(cosmology);
      
      // construct the model
      const cbl::statistics::Model1D model(&model_function2, nparameters, parType, parNames, ptr_modelInput);
      auto ptr_model = make_shared<cbl::statistics::Model1D>(model);
      
      // --- construct and maximize the likelihood ---
      const std::vector< size_t> x_index = {0,2};
      const int w_index = -1;
      const std:: shared_ptr<cbl::statistics::ModelParameters> model_parameters = NULL;
      // show_results parameters
      
      std:: cout << "  " << std:: endl;
      
      const vector<double> start = {valA, valB};
      
      // limits for the parameters as READ IN THE INPUT FILE
      const double minA = varAmin, maxA = varAmax;
      const double minB = varBmin, maxB = varBmax;
      const vector<vector<double>> limits = {{minA, maxA}, {minB, maxB}};
      
      // --- construct the priors ---
      auto prior_A = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minA, maxA));
      auto prior_B = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minB, maxB));
      const vector<shared_ptr<cbl::statistics::PriorDistribution>> prior_distributions = {prior_A, prior_B};
      
      if(file_cov == "NO"){
	std:: cout << "  " << std:: endl;
	std:: cout << " using the Gaussian Error  "   << std:: endl;
	cbl::statistics::Likelihood likelihood(ptr_data, ptr_model, cbl::statistics::LikelihoodType::_Gaussian_Error_,
					       x_index,w_index,model_parameters,prec);      // uses diagonal
      	std:: cout << "  " << std:: endl;
        std::cout << "   " << std::endl;
	
        // maximize the likelihood and write the output
	likelihood.maximize(start, limits);
	const string file_output_bestfit = "model_bestfit.dat";
	likelihood.write_model_at_bestfit(cbl::par::DirLoc+outputdir+"/",file_output_bestfit);
	
	std:: cout << "   " << std:: endl;
	
	// --- construct, maximize and sample the posterior ---
	cbl::statistics::Posterior posterior(prior_distributions, likelihood, seedlikelihood);
	
	// here is starting from the maximum of the posterior
	posterior.initialize_chains(chain_size, nwalkers, 1.e-5, {valA, valB},10000,1e-6,epsilon);
        //posterior.initialize_chains(chain_size, nwalkers);
        posterior.sample_stretch_move(2);
	
        // show the median MCMC values of the four parameters on screen
	cout << endl;
	for (size_t i=0; i<posterior.parameters()->nparameters(); ++i)
	  cout << setprecision(4) << "Posterior median of " << posterior.parameters()->name(i) << " = " << posterior.parameters()->bestfit_value(i) << endl;
	
	// posterior.show_results(burn_in, thin,nbins,show_mode,ns,nb);
	posterior.show_results(burn_in, thin,nbins,show_mode);
	
	// store the chain ouputs
	posterior.write_results(cbl::par::DirLoc+outputdir+"/", "chains_linear_relation", burn_in, thin);
	
	// store the best-fit model
	posterior.write_model_from_chain(cbl::par::DirLoc+outputdir+"/", "model_from_chain.dat", {}, {}, burn_in, thin);
      }else{
	
	std:: cout << "  " << std:: endl;
	std:: cout << " using the Gaussian Covariance  "   << std:: endl;
	cbl::statistics::Likelihood likelihood(ptr_data, ptr_model, cbl::statistics::LikelihoodType::_Gaussian_Covariance_,
					       x_index,w_index,model_parameters,prec,10000);   // uses full covariance
	
	std:: cout << "  " << std:: endl;
	
	// maximize the likelihood and write the output
	likelihood.maximize(start, limits);
	const string file_output_bestfit = "model_bestfit.dat";
	likelihood.write_model_at_bestfit(cbl::par::DirLoc+outputdir+"/",file_output_bestfit);
	
	std:: cout << "  " << std:: endl;
	
	// --- construct, maximize and sample the posterior ---
	cbl::statistics::Posterior posterior(prior_distributions, likelihood, seedlikelihood);
	
	// here is starting from the maximum of the posterior
	posterior.initialize_chains(chain_size, nwalkers, 1.e-5, {valA, valB},10000,1e-6,epsilon);
	posterior.sample_stretch_move(2);
	
	// show the median MCMC values of the four parameters on screen
	cout << endl;
	for (size_t i=0; i<posterior.parameters()->nparameters(); ++i)
	  cout << setprecision(4) << "Posterior median of " << posterior.parameters()->name(i) << " = " << posterior.parameters()->bestfit_value(i) << endl;
	
	posterior.show_results(burn_in, thin,nbins,show_mode,10000,nn);
	
	// store the chain ouputs
	posterior.write_results(cbl::par::DirLoc+outputdir+"/", "chains_linear_relation", burn_in, thin);
	
	// store the best-fit model
	posterior.write_model_from_chain(cbl::par::DirLoc+outputdir+"/", "model_from_chain.dat", {}, {}, burn_in, thin);
      }
    }
    if(npars == 4){
      // --- set the model to construct the likelihood ---
      // number of model parameters
      const int nparameters = npars;
      
      // names of the model parameters
      const vector<string> parNames = {"LogM200", "c200", "f_off", "sigma_off"};
      
      // starting values ... as READ IN THE INPUT FILE
      double valA = varA, valB = varB;
      double valC = varC, valD = varD;
      
      // vector containing the model parameters
      vector<cbl::statistics::ParameterType> parType(nparameters);
      
      auto ptr_modelInput = make_shared<cbl::cosmology::Cosmology>(cosmology);
      
      // construct the model
      const cbl::statistics::Model1D model(&model_function4, nparameters, parType, parNames, ptr_modelInput);
      auto ptr_model = make_shared<cbl::statistics::Model1D>(model);
      
      // --- construct and maximize the likelihood ---
      const std::vector< size_t> x_index = {0,2};
      const int w_index = -1;
      const std:: shared_ptr<cbl::statistics::ModelParameters> model_parameters = NULL;
      // show_results parameters
      
      std:: cout << "  " << std:: endl;
      
      const vector<double> start = {valA, valB, valC, valD};
      
      // limits for the parameters as READ IN THE INPUT FILE
      const double minA = varAmin, maxA = varAmax;
      const double minB = varBmin, maxB = varBmax;
      const double minC = varCmin, maxC = varCmax;
      const double minD = varDmin, maxD = varDmax;
      const vector<vector<double>> limits = {{minA, maxA}, {minB, maxB}, {minC, maxC}, {minD, maxD}};
      
      // --- construct the priors ---
      auto prior_A = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minA, maxA));
      auto prior_B = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minB, maxB));
      auto prior_C = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minC, maxC));
      auto prior_D = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minD, maxD));
      const vector<shared_ptr<cbl::statistics::PriorDistribution>> prior_distributions = {prior_A, prior_B, prior_C, prior_D};
      
      if(file_cov == "NO"){
	std:: cout << "  " << std:: endl;
	std:: cout << " using the Gaussian Error  "   << std:: endl;
	cbl::statistics::Likelihood likelihood(ptr_data, ptr_model, cbl::statistics::LikelihoodType::_Gaussian_Error_,
					       x_index,w_index,model_parameters,prec);      // uses diagonal

        std::cout << "  " << std::endl;
	
        // maximize the likelihood and write the output
	likelihood.maximize(start, limits);
	const string file_output_bestfit = "model_bestfit.dat";
	likelihood.write_model_at_bestfit(cbl::par::DirLoc+outputdir+"/",file_output_bestfit);
	
	std:: cout << "   " << std:: endl;
	
	// --- construct, maximize and sample the posterior ---
	cbl::statistics::Posterior posterior(prior_distributions, likelihood, seedlikelihood);
	
	// here is starting from the maximum of the posterior
	posterior.initialize_chains(chain_size, nwalkers, 1.e-5, {valA, valB, valC, valD},10000,1e-6,epsilon);
        // posterior.initialize_chains(chain_size, nwalkers);
        posterior.sample_stretch_move(2);
	
        // show the median MCMC values of the four parameters on screen
	cout << endl;
	for (size_t i=0; i<posterior.parameters()->nparameters(); ++i)
	  cout << setprecision(4) << "Posterior median of " << posterior.parameters()->name(i) << " = " << posterior.parameters()->bestfit_value(i) << endl;
	
	// posterior.show_results(burn_in, thin,nbins,show_mode,ns,nb);
	posterior.show_results(burn_in, thin,nbins,show_mode);
	
	// store the chain ouputs
	posterior.write_results(cbl::par::DirLoc+outputdir+"/", "chains_linear_relation", burn_in, thin);
	
	// store the best-fit model
	// posterior.write_model_from_chain(cbl::par::DirLoc+outputdir+"/", "model_from_chain.dat", {}, {}, burn_in, thin);
      }else{
	
	std:: cout << "  " << std:: endl;
	std:: cout << " using the Gaussian Covariance  "   << std:: endl;
	cbl::statistics::Likelihood likelihood(ptr_data, ptr_model, cbl::statistics::LikelihoodType::_Gaussian_Covariance_,
					       x_index,w_index,model_parameters,prec,10000);   // uses full covariance
	
	std:: cout << "  " << std:: endl;
	
	// maximize the likelihood and write the output
	likelihood.maximize(start, limits);
	const string file_output_bestfit = "model_bestfit.dat";
	likelihood.write_model_at_bestfit(cbl::par::DirLoc+outputdir+"/",file_output_bestfit);
	
	std:: cout << "  " << std:: endl;
	
	// --- construct, maximize and sample the posterior ---
	cbl::statistics::Posterior posterior(prior_distributions, likelihood, seedlikelihood);
	
	// here is starting from the maximum of the posterior
	posterior.initialize_chains(chain_size, nwalkers, 1.e-5, {valA, valB, valC, valD},10000,1e-6,epsilon);
	posterior.sample_stretch_move(2);
	
	// show the median MCMC values of the four parameters on screen
	cout << endl;
	for (size_t i=0; i<posterior.parameters()->nparameters(); ++i)
	  cout << setprecision(4) << "Posterior median of " << posterior.parameters()->name(i) << " = " << posterior.parameters()->bestfit_value(i) << endl;
	
	posterior.show_results(burn_in, thin,nbins,show_mode,10000,nn);
	
	// store the chain ouputs
	posterior.write_results(cbl::par::DirLoc+outputdir+"/", "chains_linear_relation", burn_in, thin);
	
	// store the best-fit model
	// posterior.write_model_from_chain(cbl::par::DirLoc+outputdir+"/", "model_from_chain.dat", {}, {}, burn_in, thin);
      }
    }
    if(npars == 5){
      // --- set the model to construct the likelihood ---
      // number of model parameters
      const int nparameters = npars;
      
      // names of the model parameters
      const vector<string> parNames = {"LogM200", "c200", "f_off", "sigma_off", "OmegaM"};
      
      // starting values ... as READ IN THE INPUT FILE
      double valA = varA, valB = varB;
      double valC = varC, valD = varD;
      double valE = varE;
      
      // vector containing the model parameters
      vector<cbl::statistics::ParameterType> parType(nparameters);
      
      auto ptr_modelInput = make_shared<cbl::cosmology::Cosmology>(cosmology);
      
      // construct the model
      const cbl::statistics::Model1D model(&model_function5, nparameters, parType, parNames, ptr_modelInput);
      auto ptr_model = make_shared<cbl::statistics::Model1D>(model);
      
      // --- construct and maximize the likelihood ---
      const std::vector< size_t> x_index = {0,2};
      const int w_index = -1;
      const std:: shared_ptr<cbl::statistics::ModelParameters> model_parameters = NULL;
      // show_results parameters
      std:: cout << "  " << std:: endl;
      
      const vector<double> start = {valA, valB, valC, valD, valE};
      
      // limits for the parameters as READ IN THE INPUT FILE
      const double minA = varAmin, maxA = varAmax;
      const double minB = varBmin, maxB = varBmax;
      const double minC = varCmin, maxC = varCmax;
      const double minD = varDmin, maxD = varDmax;
      const double minE = varEmin, maxE = varEmax;
      const vector<vector<double>> limits = {{minA, maxA}, {minB, maxB}, {minC, maxC}, {minD, maxD}, {minE, maxE}};
      
      // --- construct the priors ---
      auto prior_A = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minA, maxA));
      auto prior_B = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minB, maxB));
      auto prior_C = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minC, maxC));
      auto prior_D = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minD, maxD));
      auto prior_E = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minE, maxE));
      const vector<shared_ptr<cbl::statistics::PriorDistribution>> prior_distributions = {prior_A, prior_B, prior_C, prior_D, prior_E};
      
      if(file_cov == "NO"){
	std:: cout << "  " << std:: endl;
	std:: cout << " using the Gaussian Error  "   << std:: endl;
	cbl::statistics::Likelihood likelihood(ptr_data, ptr_model, cbl::statistics::LikelihoodType::_Gaussian_Error_,
					       x_index,w_index,model_parameters,prec);      // uses diagonal
	std:: cout << "  " << std:: endl;
	
	// maximize the likelihood and write the output
	likelihood.maximize(start, limits);
	const string file_output_bestfit = "model_bestfit.dat";
	likelihood.write_model_at_bestfit(cbl::par::DirLoc+outputdir+"/",file_output_bestfit);
	
	std:: cout << "  " << std:: endl;
	
	// --- construct, maximize and sample the posterior ---
	cbl::statistics::Posterior posterior(prior_distributions, likelihood, seedlikelihood);
	
	// here is starting from the maximum of the posterior
	posterior.initialize_chains(chain_size, nwalkers, 1.e-5, {valA, valB, valC, valD, valE},10000,1e-6,epsilon);
	posterior.sample_stretch_move(2);
	
	// show the median MCMC values of the four parameters on screen
	cout << endl;
	for (size_t i=0; i<posterior.parameters()->nparameters(); ++i)
	  cout << setprecision(4) << "Posterior median of " << posterior.parameters()->name(i) << " = " << posterior.parameters()->bestfit_value(i) << endl;
	
	// posterior.show_results(burn_in, thin,nbins,show_mode,ns,nb);
	posterior.show_results(burn_in, thin,nbins,show_mode);
	
	// store the chain ouputs
	posterior.write_results(cbl::par::DirLoc+outputdir+"/", "chains_linear_relation", burn_in, thin);
	
	// store the best-fit model
	// posterior.write_model_from_chain(cbl::par::DirLoc+outputdir+"/", "model_from_chain.dat", {}, {}, burn_in, thin);
      }else{
	
	std:: cout << "  " << std:: endl;
	std:: cout << " using the Gaussian Covariance  "   << std:: endl;
	cbl::statistics::Likelihood likelihood(ptr_data, ptr_model, cbl::statistics::LikelihoodType::_Gaussian_Covariance_,
					       x_index,w_index,model_parameters,prec,10000);   // uses full covariance
	std:: cout << "  " << std:: endl;
	
	// maximize the likelihood and write the output
	likelihood.maximize(start, limits);
	// likelihood.maximize(start, limits, max_iter, tol);
	const string file_output_bestfit = "model_bestfit.dat";
	likelihood.write_model_at_bestfit(cbl::par::DirLoc+outputdir+"/",file_output_bestfit);
	
	std:: cout << "  " << std:: endl;
	
	// --- construct, maximize and sample the posterior ---
	cbl::statistics::Posterior posterior(prior_distributions, likelihood, seedlikelihood);
	
	// here is starting from the maximum of the posterior
	posterior.initialize_chains(chain_size, nwalkers, 1.e-5, {valA, valB, valC, valD, valE},10000,1e-6,epsilon);
	posterior.sample_stretch_move(2);
	
	// show the median MCMC values of the four parameters on screen
	cout << endl;
	for (size_t i=0; i<posterior.parameters()->nparameters(); ++i)
	  cout << setprecision(4) << "Posterior median of " << posterior.parameters()->name(i) << " = " << posterior.parameters()->bestfit_value(i) << endl;
	
	posterior.show_results(burn_in, thin,nbins,show_mode,10000,nn);
	
	// store the chain ouputs
	posterior.write_results(cbl::par::DirLoc+outputdir+"/", "chains_linear_relation", burn_in, thin);
	
	// store the best-fit model
	// posterior.write_model_from_chain(cbl::par::DirLoc+outputdir+"/", "model_from_chain.dat", {}, {}, burn_in, thin);
      }
    }
    if(npars == 6){
      // --- set the model to construct the likelihood ---
      // number of model parameters
      const int nparameters = npars;
      
      // names of the model parameters
      const vector<string> parNames = {"LogM200", "c200", "f_off", "sigma_off", "OmegaM", "sigma8"};
      
      // starting values ... as READ IN THE INPUT FILE
      double valA = varA, valB = varB;
      double valC = varC, valD = varD;
      double valE = varE, valF = varF;
      
      // vector containing the model parameters
      vector<cbl::statistics::ParameterType> parType(nparameters);
      
      auto ptr_modelInput = make_shared<cbl::cosmology::Cosmology>(cosmology);
      
      // construct the model
      const cbl::statistics::Model1D model(&model_function6, nparameters, parType, parNames, ptr_modelInput);
      auto ptr_model = make_shared<cbl::statistics::Model1D>(model);
      
      // --- construct and maximize the likelihood ---
      const std::vector< size_t> x_index = {0,2};
      const int w_index = -1;
      const std:: shared_ptr<cbl::statistics::ModelParameters> model_parameters = NULL;
      // show_results parameters
      std:: cout << "  " << std:: endl;
      
      const vector<double> start = {valA, valB, valC, valD, valE, valF};
      
      // limits for the parameters as READ IN THE INPUT FILE
      const double minA = varAmin, maxA = varAmax;
      const double minB = varBmin, maxB = varBmax;
      const double minC = varCmin, maxC = varCmax;
      const double minD = varDmin, maxD = varDmax;
      const double minE = varEmin, maxE = varEmax;
      const double minF = varFmin, maxF = varFmax;
      const vector<vector<double>> limits = {{minA, maxA}, {minB, maxB}, {minC, maxC}, {minD, maxD}, {minE, maxE}, {minF, maxF}};
      
      // --- construct the priors ---
      auto prior_A = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minA, maxA));
      auto prior_B = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minB, maxB));
      auto prior_C = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minC, maxC));
      auto prior_D = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minD, maxD));
      auto prior_E = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minE, maxE));
      auto prior_F = make_shared<cbl::statistics::PriorDistribution>(cbl::statistics::PriorDistribution(cbl::glob::DistributionType::_Uniform_, minF, maxF));
      const vector<shared_ptr<cbl::statistics::PriorDistribution>> prior_distributions = {prior_A, prior_B, prior_C, prior_D, prior_E, prior_F};
      
      if(file_cov == "NO"){
	std:: cout << "  " << std:: endl;
	std:: cout << " using the Gaussian Error  "   << std:: endl;
	cbl::statistics::Likelihood likelihood(ptr_data, ptr_model, cbl::statistics::LikelihoodType::_Gaussian_Error_,
					       x_index,w_index,model_parameters,prec);      // uses diagonal
	std:: cout << "  " << std:: endl;
	
	// maximize the likelihood and write the output
	likelihood.maximize(start, limits);
	const string file_output_bestfit = "model_bestfit.dat";
	likelihood.write_model_at_bestfit(cbl::par::DirLoc+outputdir+"/",file_output_bestfit);
	
	std:: cout << "  " << std:: endl;
	
	// --- construct, maximize and sample the posterior ---
	cbl::statistics::Posterior posterior(prior_distributions, likelihood, seedlikelihood);
	
	// here is starting from the maximum of the posterior
	posterior.initialize_chains(chain_size, nwalkers, 1.e-5, {valA, valB, valC, valD, valE, valF},10000,1e-6,epsilon);
	posterior.sample_stretch_move(2);
	
	// show the median MCMC values of the four parameters on screen
	cout << endl;
	for (size_t i=0; i<posterior.parameters()->nparameters(); ++i)
	  cout << setprecision(4) << "Posterior median of " << posterior.parameters()->name(i) << " = " << posterior.parameters()->bestfit_value(i) << endl;
	
	// posterior.show_results(burn_in, thin,nbins,show_mode,ns,nb);
	posterior.show_results(burn_in, thin,nbins,show_mode);
	
	// store the chain ouputs
	posterior.write_results(cbl::par::DirLoc+outputdir+"/", "chains_linear_relation", burn_in, thin);
	
	// store the best-fit model
	// posterior.write_model_from_chain(cbl::par::DirLoc+outputdir+"/", "model_from_chain.dat", {}, {}, burn_in, thin);
      }else{
	
	std:: cout << "  " << std:: endl;
	std:: cout << " using the Gaussian Covariance  "   << std:: endl;
	cbl::statistics::Likelihood likelihood(ptr_data, ptr_model, cbl::statistics::LikelihoodType::_Gaussian_Covariance_,
					       x_index,w_index,model_parameters,prec,10000);   // uses full covariance
	std:: cout << "  " << std:: endl;
	
	// maximize the likelihood and write the output
	likelihood.maximize(start, limits);
	// likelihood.maximize(start, limits, max_iter, tol);
	const string file_output_bestfit = "model_bestfit.dat";
	likelihood.write_model_at_bestfit(cbl::par::DirLoc+outputdir+"/",file_output_bestfit);
	
	std:: cout << "  " << std:: endl;
	
	// --- construct, maximize and sample the posterior ---
	cbl::statistics::Posterior posterior(prior_distributions, likelihood, seedlikelihood);
	
	// here is starting from the maximum of the posterior
	posterior.initialize_chains(chain_size, nwalkers, 1.e-5, {valA, valB, valC, valD, valE, valF},10000,1e-6,epsilon);
	posterior.sample_stretch_move(2);

	// show the median MCMC values of the four parameters on screen
	cout << endl;
	for (size_t i=0; i<posterior.parameters()->nparameters(); ++i)
	  cout << setprecision(4) << "Posterior median of " << posterior.parameters()->name(i) << " = " << posterior.parameters()->bestfit_value(i) << endl;

	posterior.show_results(burn_in, thin,nbins,show_mode,10000,nn);

	// store the chain ouputs
	posterior.write_results(cbl::par::DirLoc+outputdir+"/", "chains_linear_relation", burn_in, thin);

	// store the best-fit model
	// posterior.write_model_from_chain(cbl::par::DirLoc+outputdir+"/", "model_from_chain.dat", {}, {}, burn_in, thin);
      }
    }
  }
  catch(cbl::glob::Exception &exc) { cerr << exc.what() << endl; exit(1); }
  std:: cout << " ... end of work ... ;-) " << std:: endl;
  return 0;
}
