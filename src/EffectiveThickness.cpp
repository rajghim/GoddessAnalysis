#include "EffectiveThickness.h"
#include <iostream>
#include <fstream>
#include <math.h> 


//This function converts cartesion co-ordinates to sphercial polar co-ordinates
void carttosph(float x, float y, float z, float *r, float *theta, float *phi)
{
    float temp;
    if (x==0. && y==0. && z==0.){
        *r = 0;
        *theta = 0;
        *phi = 0;
    } else {
        *r = std::sqrt((std::pow(x,2)) + (std::pow(y,2)) + (std::pow(z,2)));
        float R = std::sqrt((std::pow(x,2)) + (std::pow(y,2)) + (std::pow(z,2)));
        *theta = acos(z/R);
        if (x==0. && y==0){
            *phi = 0.;
        }else{
            temp = std::sqrt((std::pow(x,2.)) + (std::pow(y,2.)));
            if (x>=0.) *phi = acos(y/temp);
            if (x<0.) *phi = 2*M_PI - acos(y/temp);
        }
    }
}


//This function converts spherical polar co-ordinates to cartesian co-ordinates
void sphtocart(float r, float theta, float phi, float *x, float *y, float *z){
    *z = r*cos(theta);
    float temp = r*sin(theta);
    *x = temp*sin(phi);
    *y = temp*cos(phi);
}


//This function outputs effective thickness
void targ_thick(float theta_in, float phi_in, float thickness_Z, float depth, float theta_targ, float *eff_thickness){
    //Note: theta_targ = 27 degrees for target, 0 degrees for QQQ5 deadlayer and 90 degrees for SX3 detectors
    //Converting the ion's direction to cartesian coordinates
    float dummy = 1.0;
    float x,y,z;
    sphtocart(dummy, theta_in, phi_in, &x, &y, &z);
    float length = (std::pow(x,2.)) + (std::pow(y,2.)) + (std::pow(z,2.));

    //Rotating the ion's vector so it is measured w.r.to the target
    float newy = y*cos(theta_targ) - z*sin(theta_targ);
    float newz = z*cos(theta_targ) + y*sin(theta_targ);
    float newx = x;

    //Convert the ion's vector back to spherical polars, now wrt the target
    float dummy_r, dummy_theta, dummy_phi;
    carttosph(newx, newy, newz, &dummy_r, &dummy_theta, &dummy_phi);

    //Calculating the thickness seen by the ion
    if(dummy_theta <= (0.5*M_PI)) *eff_thickness = (thickness_Z-depth)/cos(dummy_theta);
	if(dummy_theta > (0.5*M_PI)) *eff_thickness = (-depth)/cos(dummy_theta);
}
