 #include "udf.h"
 #include <math.h>

 DEFINE_PROFILE(U_Slip_X_Profile, thread, position)
 {
    #if !RP_HOST
		real FaceArea[ND_ND];
		real Cell_Center[ND_ND];
		real Face_Center[ND_ND];
		real n1,n2;
		real *T_Grad=NULL;
		real *U_Grad=NULL;
		real q1,q2;
		real Txx,Txy,Tyx,Tyy;
		real U_Slip_X;
		real k;
		real rho;
		real mu;
		real Cp;
		const real pi=M_PI;
		real P;
		const real u=0.4987445;
		real lambda;
		real Pr;
		face_t f;
		cell_t c0;
		Thread *t0=NULL;
		real gamma = 1.4;
		real sigma=0.9;
		real U_old,U_new;
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
		        
		        /*Density @ Cell */   
		        rho=C_R(c0,t0);
		        
		        /*Pressure @ Cell*/
		        P=C_P(c0,t0);
		        
		        /* Vector of Temp Gradient @ Cell*/
		        T_Grad=C_T_RG(c0,t0);    
		        
		        /* Vector of U velocity Gradient @ Cell*/
		        U_Grad=C_U_RG(c0,t0);                               
		        
		        /*Thermal Conductivity @ Cell*/
		        k=C_K_L(c0,t0);                                         
		        
		        /*Wall Heat Flux \vec{q}=q1 \hat{i} + q2 \hat{j}*/
		        q1=-k*T_Grad[0];                                      
		        q2=-k*T_Grad[1];
		        
		        /*Viscosity @ Cell*/
		        mu=C_MU_L(c0,t0);
		        
		        /*Stress Tensor near the wall */
		            /*Tau_XX*/
		            Txx=-mu *2.0*C_DUDX(c0,t0);
		            /*Tau_XY*/
		            Txy=-mu*(C_DVDX(c0,t0)+C_DUDY(c0,t0));
		            /*Tau_YX=Tau_XY*/
		            Tyx=Txy;
		            /*Tau_YY*/
		            Tyy=-mu *2.0*C_DVDY(c0,t0);
		        
		        /*Specific Heat @ Cell*/
		        Cp=C_CP(c0,t0);
		        	    
		        /*Prandtl number*/
		        Pr=Cp*mu/k;		    
       
		        /*Mean Free Path @ Cell from Jennings(1987)*/
		        lambda=sqrt(pi/8.0)*mu/u*1.0/sqrt(rho*P);
		        
		        /*Current U velocity*/
		        U_old=C_U(c0,t0)+(Face_Center[0]-Cell_Center[0])*U_Grad[0]+(Face_Center[1]-Cell_Center[1])*U_Grad[1];
		        	        
		        /*X Slip velocity at wal*/	    
		        U_new=3.0*Cp*mu*(pow(n1,2.0)*q1+n1*n2*q2-q1)*(gamma-1.0)/(4.0*k*P*gamma)
		                 -
		                 (-n1*Txx+pow(n1,3.0)*Txx+pow(n1,2.0)*n2*Txy-n2*Tyx+pow(n1,2.0)*n2*Tyx+n1*pow(n2,2.0)*Tyy)*lambda*(sigma-2.0)/(mu*sigma);
                 		    
		        /*Limiting Slip values*/
		        if(U_new<=-500.0)
	            {
	                Message("U_Slip_X= %g is limited to -500.0\n",U_new);
		        	U_new=0.0;
	        	}		    	
		       	if(U_new>500.0)
		       	{
		       	    Message("U_Slip_X= %g is limited to 500.0\n",U_new);
		       		U_new=500.0;
		        }
		      	
		      	/*Relaxing slip values*/
		      	U_Slip_X=(1.0-Rlx_f)*U_old+Rlx_f*U_new;
		      	
		        F_PROFILE(f,thread,position)=  U_Slip_X;
		    }
		    end_f_loop(f, thread)
	    }
	    else
	    {
	        Message("Data not Valid for U: Setting U=0.0\n");
		    begin_f_loop(f,thread)
		    {
		        U_Slip_X=0.0;
		        F_PROFILE(f,thread,position)=  U_Slip_X;
		    }
		    end_f_loop(f, thread)
	    } 
	#endif
 }
