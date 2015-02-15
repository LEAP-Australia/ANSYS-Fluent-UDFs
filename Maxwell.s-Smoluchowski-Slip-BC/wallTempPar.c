/*
Copyright (c) 2014, LEAP Pty Ltd
All rights reserved.

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following
conditions are met:

1. Redistributions of source code must retain the above
copyright notice, this list of conditions and the following
disclaimer.

2. Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided
with the distribution.

3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived
from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Discription:
Using DEFINE_PROFILE Macro to Implement Maxwell's Slip and Smoluchowski
Temperature Jump for Density Based Solver

This example a real values scheme variable named velocity-inlet is
defined in scheme. This value is used to do define a profile.
*/
#include "udf.h"
#include <math.h>

 DEFINE_PROFILE(wall_temp, thread, position)
 {
   
    #if !RP_HOST
		real FaceArea[ND_ND];
        real Cell_Center[ND_ND];
		real Face_Center[ND_ND];
		real n1,n2;
		real *T_Grad=NULL;
		real T_Grad_N;
		real FaceTemp,FaceTemp_new,FaceTemp_old;
		real rho;
		real mu;
		const real pi=M_PI;
		real P;
		const real u=0.4987445;
		real lambda;
		real fact;
		real Pr, k, Cp;
		real sigma=0.9;
		real gamma=1.4;		
		face_t f;
		cell_t c0;
		Thread *t0=NULL;
		real Rlx_f=0.1;

        if(Data_Valid_P())
		{
		    begin_f_loop(f,thread)
		    {
		        /*Cell Info corresponding to Boundary face*/
		        c0=F_C0(f,thread);
		        t0=THREAD_T0(thread);
		        
		        /*Cell center location*/
		        C_CENTROID(Cell_Center,c0,t0);
		        
		        /*Face center location*/
		        F_CENTROID(Face_Center,f,thread);

		        /*Vector of face area*/
		        F_AREA(FaceArea, f, thread);
		        
		        /*Normal vector of face in direction into the domain*/
		        n1=-FaceArea[0]/sqrt(pow(FaceArea[0],2.0)+pow(FaceArea[1],2.0));
		        n2=-FaceArea[1]/sqrt(pow(FaceArea[0],2.0)+pow(FaceArea[1],2.0));
	     
		        /* Vector of Temp Gradient @ Cell*/
		        T_Grad=C_T_RG(c0,t0);
		        
		        /*Density @ Cell */   
		        rho=C_R(c0,t0);
		        
		        /*Pressure @ Cell*/
		        P=C_P(c0,t0);
		        
		        /*Viscosity @ Cell*/
		        mu=C_MU_L(c0,t0);
		        
		        /*Thermal Conducitivity @ Cell*/
		        k=C_K_L(c0,t0);
		        
		        /*Specific Heat @ Cell*/
		        Cp=C_CP(c0,t0);
		        
		        /*Prandtl number*/
		        Pr=Cp*mu/k;	
		        
		        /*Mean Free Path @ Cell from Jennings(1987)*/
		        lambda=sqrt(pi/8.0)*mu/u*1.0/sqrt(rho*P);
		        
		        /*Current wall Temp value*/
		        FaceTemp_old=C_T(c0,t0)+(Face_Center[0]-Cell_Center[0])*T_Grad[0]+(Face_Center[1]-Cell_Center[1])*T_Grad[1];
		        		        		    
		        /*Temp gradient in the direction normal to face*/		    	
		        T_Grad_N=n1*T_Grad[0]+n2*T_Grad[1];
		        
		        FaceTemp_new=550.0+(2.0-sigma)/(sigma)*2.0*gamma/(gamma+1.0)*lambda/Pr*T_Grad_N;
		        
		        if(FaceTemp_new<=0.0)
		        	FaceTemp_new=550.0;
		        	
		       	if(FaceTemp_new>600)
		       		FaceTemp_new=550.0;
		        
		      	FaceTemp=(1.0-Rlx_f)*FaceTemp_old+Rlx_f*FaceTemp_new;
		      	
		        F_PROFILE(f,thread,position)= FaceTemp;
	        }
		    end_f_loop(f, thread)
	    }
	    else
	    {
	        Message("Data not Valid for T: Setting T=550.0\n");
		    begin_f_loop(f,thread)
		    {
		        FaceTemp=550.0;
		        F_PROFILE(f,thread,position)= FaceTemp;
		    }
		    end_f_loop(f, thread)
	    }
	#endif
 }
