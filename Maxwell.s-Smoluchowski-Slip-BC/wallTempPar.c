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
