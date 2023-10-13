#include <time.h>

#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phaseEVP.h"
#include "log-conform.h"
#include "tension.h"
#include "fluidlab_pack.h"

int mesh_level = 11;

// Parameters of the navier-stokes and saramito physics model
double Re = 0.01; // Reynolds number
double Ca = 0.0; // Weber number
double Bi = 0.0; // Bingham number
double Wi = 0.01; // Weissenberg number
double beta = 0.11111; // Viscosity ratio
double ratio_density = 1.0; // Ratio between density of two fluids
double ratio_viscosity = 1.0; // Ratio between density of two fluids

// Parameters of the numerical method
double regularization = 1e-7; // Yield-stress regularization parameter
double domain_size = 5.0; // Total size of the square domain

// Fields used by basilisk and for visualization
scalar lambdav[], mupv[]; // polymer relaxation time and polymeric viscosity (nondimensional: (Wi) and (1-beta))
scalar yielded[]; // Tracking yielded/plugged regions. Only used for visualization
scalar norm_tau_dev_field[]; // Storing the norm of tau_dev. Only used for visualization
scalar D2[]; // Storing the norm of strain rate. Only used for visualization

// Solid boundary on top and bottom
u.t[top] = dirichlet(y);
u.t[bottom] = dirichlet(y);

int main (int argc, char * argv[]) {

  if( argc<11 ) {
    printf("Not enough input arguments..\n");
    exit(0);
  }

  // ./droplet_shear 0.05 0.3 0.8 2.5 0.5 1e-07 16.0 1e-04 0.01 8
  Re = atof(argv[1]);
  Ca = atof(argv[2]);
  Bi = atof(argv[3]);
  Wi = atof(argv[4]);
  beta = atof(argv[5]);
  regularization = atof(argv[6]);
  domain_size = atof(argv[7]);
  TOLERANCE = atof(argv[8]);
  DT = atof(argv[9]);
  mesh_level = atof(argv[10]);
  ratio_density = 1.0;
  ratio_viscosity = 1.0;


  printf("Re = %g\n", Re);
  printf("Ca = %g\n", Ca);
  printf("Bi = %g\n", Bi);
  printf("Wi = %g\n", Wi);
  printf("beta = %g\n", beta);
  printf("regularization = %g\n", regularization);
  printf("domain_size = %g\n", domain_size);
  printf("tolerance = %g\n", TOLERANCE);
  printf("dt_max = %g\n", DT);

  OpenSimulationFolder("new_evp_droplet_mesh%d_Bi%g_Wi%g_beta%g_Re%g_Ca%g", mesh_level, Bi, Wi, beta, Re, Ca);

  // If Wi==0, we just make beta=1, so all the viscosity comes from the generalized newtonian formulation
  // that is, we just ignore the polymeric stress completely and do a traditional Bingham (generalized newtonian)
  if( Wi==0.0 )
    beta = 1.0;

  // Creating the domain (rectangle [-10, 10] x [-1, 1])
  size(domain_size);
  origin(-0.5*domain_size, -0.5*domain_size);
  init_grid(1 << 8);

  // Periodic boundaries left/right
  periodic(left);

  // Setting density and viscosities coefficients for both fluids
  // Fluid 1: EVP droplet
  // Fluid 2: surrounding air, newtonian
  rho2 = Re;
  mu2 = 1.0; 
  rho1 = Re*ratio_density;
  mu1 = beta*ratio_viscosity;
  

  // Setting the surface tension coefficient
  f.sigma = 1.0/(Ca);

  /**
  The viscoelastic fields (polymeric viscosity and relaxation time) will be set in the properties event. */
  mup = mupv;
  lambda = lambdav;

  run();
  CloseSimulationFolder();
}

double initial_shape(double x, double y) {
  return - sq(x) - sq(y) + sq(1.0);
}

// === This function was implemented by Tom Appleford
void CalculateDropletDeformation(double *deformation1, double *deformation2, double *deformation3)
{
  // We begin by calculating the center of mass (xcm, ycm)
  double wt = 0., xcm = 0., ycm = 0.;
  foreach(reduction(+:wt), reduction(+:xcm), reduction(+:ycm)) {
    wt += f[]*sq(Delta);
    xcm += f[]*x*sq(Delta);
    ycm += f[]*y*sq(Delta);
  }
  xcm /= wt;
  ycm /= wt;


  // Now we calculate the deformation with one method (by Vatsal and Tom)
  double r_max2 = 0., r_min2 = HUGE; // rmin is intialized with a large number
  face vector s; // need s for normal construction (see below)
  s.x.i = -1; // just book keeping -- ignore
  foreach(reduction(max:r_max2), reduction(min:r_min2)) {
    if (y > ycm) {// assuming rotational symmetry
      // check whether or not this is an interfacial cell
      if (f[] > 1e-6 && f[] < 1. - 1e-6) {
        coord n = facet_normal (point, f, s); // construct normal to the interface
        double alpha = plane_alpha (f[], n); // construct the intercept. remember y = mx+ alpha
        coord segment[2]; // where the interfacial coordinates are saved

        // we should only move forward if there both points of the interfacial line segment is in this cell. We want to skip cells which merely touches the interface
        if (facets (n, alpha, segment) == 2) {
          /* 
          the two points on the interface are essentially segment[0] and segment[1].
          So, we can find the centroid of the interface in this cell.
          */
          double xc = x + (segment[0].x+segment[1].x)*Delta/2.;
          double yc = y + (segment[0].y+segment[1].y)*Delta/2.;

          // calculate the distance of this interfacial location from the ceneter of mass of the droplet
          double rTemp = sqrt(sq(xc-xcm) + sq(yc-ycm));

          // checl for rmax2 and rmin2
          if (rTemp > r_max2) {
            r_max2 = rTemp;
            // xImax = xc; yImax = yc;
          }
          if (rTemp < r_min2) {
            r_min2 = rTemp;
            // xImin = xc; yImin = yc;
          }
        }
      }
          
    }
  }
  *deformation1 = (r_max2/r_min2 - 1)/(r_max2/r_min2 + 1);


  // Now we calculate again with another method, just to double check
  double fmin = 1e-3; // do not reconstruct fragments smaller than this
  double r_max = -1e+10, r_min = 1e+10;
  foreach(serial) {
    if( cfilter (point, f, fmin) ) {
      coord n = interface_normal (point, f);
      double alpha = plane_alpha (f[], n);
      coord v[2];
      int m = facets (n, alpha, v);
      if( m==2 ) {
        double p1x = x + v[0].x*Delta;
        double p1y = y + v[0].y*Delta;
        double p2x = x + v[0].x*Delta;
        double p2y = y + v[0].y*Delta;
        double centerx = 0.5*(p1x + p2x);
        double centery = 0.5*(p1y + p2y);

        double radius = sqrt( sq(centerx - xcm) + sq(centery - ycm) );
        r_max = (radius>r_max) ? radius : r_max;
        r_min = (radius<r_min) ? radius : r_min;
      }
    }
  }
  *deformation2 = (r_max - r_min)/(r_max + r_min);

  // Yet another method because why not, let's triple check
  double rmax = -HUGE, rmin = HUGE;
  foreach (reduction(max:rmax) reduction(min:rmin)) {
    if (f[] > 0 && f[] < 1) {
      coord p;
      coord n = mycs (point, f);
      double alpha = plane_alpha (f[], n);
      plane_area_center (n, alpha, &p);
      double rad  = sqrt(sq(x + Delta*p.x - xcm) + sq(y + Delta*p.y - ycm)); 
      if (rad > rmax)
	      rmax = rad;
      if (rad < rmin)
	      rmin = rad;
    }
  }
  *deformation3 = (rmax - rmin)/(rmax + rmin);

  delete((scalar *){s});
  return;
}

// Measures how long it took (real time) to run the latest timestep
clock_t start_time = -1.0, end_time = -1.0;
double cpu_time_used = 0.0; // Im actually measuring real time
event step_cputime(i++) {
  if( pid() )
    return 0;

  end_time = clock();
  cpu_time_used  = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
  start_time = end_time;
}

event logfile (t+=0.01; t<=50.0) {

  // === Calculating the kinetic energy inside the droplet and D2
  double kinetic = 0.0;
  foreach (reduction(+:kinetic)){
    kinetic += f[]*0.5*(sq(u.x[]) + sq(u.y[]))*sq(Delta);

    double D11 = (u.x[1, 0] - u.x[-1, 0])/(2.0*Delta);
    double D22 = (u.y[0, 1] - u.y[0, -1])/(2.0*Delta);
    double D12 = 0.5*( (u.x[0, 1] - u.x[0, -1])/(2.0*Delta) + (u.y[1, 0] - u.y[-1, 0])/(2.0*Delta) );
    D2[] = sqrt( sq(D11) + sq(D22) + 2.0*sq(D12) );
  }

  // === Calculating the droplet deformation
  double deformation1, deformation2, deformation3;
  CalculateDropletDeformation(&deformation1, &deformation2, &deformation3);

  PrintLog("%d %lf %e %lf %e %lf %lf %lf\n", i, t, dt, cpu_time_used, kinetic, deformation1, deformation2, deformation3);
}

event print_solution_mesh(t+=0.1) {
  // scalar *list = {f, u.x, u.y, p, tau_p.x.x, tau_p.x.y, tau_p.y.y, yielded, norm_tau_dev_field, D2};
  // PrintMeshVTK_Binary_Float(i, t, list, 
  //     (const char *[]){"fractions", "vel-u", "vel-v", "pressure", "taup_xx", "taup_xy", "taup_yy", "yielded", "norm_tau_dev", "D2"});

  PrintInterfaceVTK(i, t);
}

event print_solution_mesh2(t+=10.0) {
  scalar *list = {f, u.x, u.y, p, tau_p.x.x, tau_p.x.y, tau_p.y.y, yielded, norm_tau_dev_field, D2};
  PrintMeshVTK_Binary_Float(i, t, list, 
      (const char *[]){"fractions", "vel-u", "vel-v", "pressure", "taup_xx", "taup_xy", "taup_yy", "yielded", "norm_tau_dev", "D2"});
}

event init (t = 0) {

  // === Masking the domain into a rectangle
  mask (y >= 4 ? top : none);
  mask (y <= -4 ? bottom : none);

  // === Refining mesh at the center of the domain where the droplet is placed
  refine( (sq(x) + sq(y) < sq(1.0)) && (level<mesh_level) );

  // === Filling in the volume fraction field
  fraction (f, initial_shape(x, y));

  // === Initializing the velocity field
  foreach()
    u.x[] = y;

  boundary(all);
}

#if TREE
event adapt (i++) {
  // Mesh adaptation
  adapt_wavelet ({f, u.x, u.y}, (double[]){1e-2, 1e-03, 1e-03}, maxlevel = mesh_level, minlevel = mesh_level - 2);
}
#endif

event properties (i++) {

  // If Wi = 0 exactly, then I will just use a pure Bingham model (generalized newtonian) with regularization. (is this cheating?)
  // This section was partially taken from Vatsal's sandbox (I changed the regularization method).
  if( Wi==0.0 ) {
    foreach_face(x) {
      double ff = (sf[] + sf[-1])/2.;
      alphav.x[] = fm.x[]/rho(ff);
      face vector muv = mu;
      double D11 = 0.5*(u.y[0,1] - u.y[0,-1] + u.y[-1,1] - u.y[-1,-1])/(2.0*Delta);
      double D22 = (u.x[] - u.x[-1,0])/Delta;
      double D12 = 0.5*( (u.y[] - u.y[-1, 0])/Delta + 0.5*(u.x[0,1] - u.x[0,-1] + u.x[-1,1] - u.x[-1,-1])/(2.0*Delta) );
      double norm_D = sqrt( sq(D11) + sq(D22) + 2.0*sq(D12) );
      double apparent_visc = 1.0*ratio_viscosity + Bi/(2.0*norm_D + 1e-10)*(1.0 - exp(-norm_D/regularization));
      muv.x[] = fm.x[]*( clamp(ff, 0.0, 1.0)*(apparent_visc - mu2) + mu2 );
    }

    foreach_face(y) {
      double ff = (sf[0,0] + sf[0,-1])/2.;
      alphav.y[] = fm.y[]/rho(ff);
      face vector muv = mu;
      double D11 = (u.y[0,0] - u.y[0,-1])/Delta;
      double D22 = 0.5*( (u.x[1,0] - u.x[-1,0] + u.x[1,-1] - u.x[-1,-1])/(2.*Delta) );
      double D12 = 0.5*( (u.x[0,0] - u.x[0,-1])/Delta + 0.5*( (u.y[1,0] - u.y[-1,0] + u.y[1,-1] - u.y[-1,-1])/(2.*Delta) ) );
      double norm_D = sqrt( sq(D11) + sq(D22) + 2.0*sq(D12) );
      double apparent_visc = 1.0*ratio_viscosity + Bi/(2.0*norm_D + 1e-10)*(1.0 - exp(-norm_D/regularization));
      muv.y[] = fm.y[]*( clamp(ff, 0.0, 1.0)*(apparent_visc - mu2) + mu2 );
    }

    // Regarding the oldroyd-B solver
    // We just want to force tau_p = 0, so I set mupv=0 and lambdav=1 (lambda is irrelevant if mupv=0)
    foreach() {
      mupv[] = 0.0;
      lambdav[] = 1.0; // This value doesnt really matter because mupv=0
    }

  }
  else { // If Wi!=0, we really have to solve for the polymeric part of the stress tensor
    
    // This loop will set lambda and mup (relaxation time and polymeric viscosity)
    foreach() {
      /// === Norm of the polymeric stress tensor
      double trace_taup = tau_p.x.x[] + tau_p.y.y[];
      double tau_dev_xx = tau_p.x.x[] - 0.5*trace_taup; 
      double tau_dev_yy = tau_p.y.y[] - 0.5*trace_taup; 
      double tau_dev_xy = tau_p.x.y[];
      double norm_tau_dev = sqrt( 1.0*(tau_dev_xx*tau_dev_xx + 2.0*tau_dev_xy*tau_dev_xy + tau_dev_yy*tau_dev_yy) );
      double yield_term = 1.0;
      yielded[] = 1.0;
      if( Bi>0.0 ) {
        yield_term = (norm_tau_dev < (Bi + regularization)) ? regularization : (norm_tau_dev - Bi)/norm_tau_dev;
        yielded[] = (norm_tau_dev < (Bi + regularization)) ? 0.0 : 1.0;
      }
      // mupv[] = ratio_viscosity*(1. - beta)*clamp(f[],0,1)/(yield_term);
      mupv[] = ratio_viscosity*(1. - beta)*clamp(f[],0,1)/(yield_term);
      lambdav[] = (Wi/yield_term)*clamp(f[],0,1);
    }

    // This loop will set solvent viscosity (copy-paste from two-phase.h)
    foreach_face() {
      double ff = (sf[] + sf[-1])/2.;
      alphav.x[] = fm.x[]/rho(ff);
      if (mu1 || mu2) {
        face vector muv = mu;
        muv.x[] = fm.x[]*mu(ff);
      }
    }
  }

  foreach()
    rhov[] = cm[]*rho(sf[]);

#if TREE  
  sf.prolongation = fraction_refine;
  sf.dirty = true; // boundary conditions need to be updated
#endif

  boundary(all);
}


