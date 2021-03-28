#include <stdlib.h>
#include <stdio.h>
#include "FTAPI.h"
#include "constants.h"
#include <assert.h>

static double flash_level_func( void* func_params, double *coords) 
{
    double level;

    // Call the implementation from Simulation tree
    simulation_frontlevelfunc(coords, &level);
    return level;
}

static int fr_gmFlash2FT(int flashGM)
{
    int ftgm;
    switch(flashGM)
    {
	case CARTESIAN:
		return FT_GEOMETRY_CARTESIAN;
	case POLAR:
		return FT_GEOMETRY_POLAR;
	case CYLINDRICAL:
		return FT_GEOMETRY_CYLINDRICAL;
	case SPHERICAL:
		return FT_GEOMETRY_SPHERICAL;
	default:
		assert(0);
    }
}

static int fr_bcFlash2FT(int flashBC)
{
	int ftBC;
	if (flashBC==PERIODIC || flashBC==NOT_BOUNDARY)
	    ftBC = FT_BOUNDARY_PERIODIC;
	else if (flashBC==REFLECTING)
	    ftBC =  FT_BOUNDARY_REFLECTING;
	else if (flashBC == OUTFLOW)
	    ftBC = FT_BOUNDARY_FLOW;
	else
	{
	    printf("FrontTracking: Unsupported boundary type %d\n",flashBC);
	    assert(0);
	}
	//printf("Got BC: %s\n",f_wave_type_as_string(ftBC));
	return ftBC;
}

void FNAME(fr_init, FR_INIT)(int *dim,
	int * procGrid,
        double *xmin, double*xmax, double *ymin, double *ymax, 
        double *zmin, double *zmax,
        int *isize, int *jsize, 
        int *ksize,
	int *ibuf, int *jbuf, int *kbuf,
        int *xlbdry, int *xubdry,
        int *ylbdry, int *yubdry,
        int *zlbdry, int *zubdry,
        double *vx, double *vy,
        double *vz, int *geometry, int *restart,
	char *restart_filename,
	int restart_filename_length)
{
    // Look up FronTier enum types for boundary conditions
    int xlbdryft = fr_bcFlash2FT(*xlbdry);
    int xubdryft = fr_bcFlash2FT(*xubdry);
    int ylbdryft = fr_bcFlash2FT(*ylbdry);
    int yubdryft = fr_bcFlash2FT(*yubdry);
    int zlbdryft = 0; 
    int zubdryft = 0; 
    if(*dim==3)
    {
	zlbdryft = fr_bcFlash2FT(*zlbdry);
	zubdryft = fr_bcFlash2FT(*zubdry);
    }

    // Look up FronTier enum type for geometry
    int geometryft = fr_gmFlash2FT(*geometry);

    // Call FronTier Initialization function
    ftapi_init_cartesian(dim, procGrid,
         xmin, xmax,
         ymin,  ymax,
         zmin,  zmax,
         isize,  jsize,  ksize,
         ibuf,  jbuf,  kbuf,
         &xlbdryft,  &xubdryft,
         &ylbdryft,  &yubdryft,
         &zlbdryft,  &zubdryft,
         vx,  vy,  vz,
	 flash_level_func, NULL,
	 &geometryft, restart, restart_filename,
	 restart_filename_length);

}


/*void debug_point(POINT *p, HYPER_SURF_ELEMENT *hse, HYPER_SURF *hs)
{
        double tol = sqr(1.e-5);
        double coords1[3] = {-0.241666916780643892,   0.00998216360650780628, -0.25};
        double coords2[3] ={-0.241666916780643892 + 0.5,   0.00998216360650780628, -0.25};
        double dist1 = sqr(Coords(p)[0] - coords1[0]) + sqr(Coords(p)[1] - coords1[1]) + sqr(Coords(p)[2] - coords1[2]);
        double dist2 = sqr(Coords(p)[0] - coords2[0]) + sqr(Coords(p)[1] - coords2[1]) + sqr(Coords(p)[2] - coords2[2]);
        
        double vel[3];
        if(dist1 < tol)
        {
            flash_vel_intrp( front.vparams, &front, p, hse, hs, vel);
            printf("P1: coords = %f %f %f : vel = %g %g %g\n",
            Coords(p)[0], Coords(p)[1], Coords(p)[2],
            vel[0], vel[1],vel[2]);
        }
        
        if(dist2 < tol)
        {
            flash_vel_intrp( front.vparams, &front, p, hse, hs, vel);
            printf("P2: coords = %f %f %f : vel = %g %g %g\n",
            Coords(p)[0], Coords(p)[1], Coords(p)[2],
            vel[0], vel[1],vel[2]);
        }
}*/                   
