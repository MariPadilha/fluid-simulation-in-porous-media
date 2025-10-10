#include "comum.h"

double maior_valor(double **array, int linha, int coluna, int i, int j){
    double max = fabs(array[i][j]);
    for(int i = 3; i < linha; i++){
        for(int j = 2; j < coluna; j++){
            if(array[i][j] > max)
                max = fabs(array[i][j]);
        }
    }
    return max;
}

//--- upwind_U ---
void upwind_Ui(double **um, double **vm, double **p, double **ru, int j){
    int i;

    for(i = 2; i <= imax; i++){
        //compute x-direction velocity component un
        fn[i][j] = 0.5 * (vm[i][j+1]+vm[i-1][j+1]) * areau_n[i] / epsilon1[i][j];
        fs[i][j] = 0.5 * (vm[i][j]+vm[i-1][j]) * areau_s[i] / epsilon1[i][j];
        fe[i][j] = 0.5 * (um[i+1][j]+um[i][j]) * areau_e[j] / epsilon1[i][j];
        fw[i][j] = 0.5 * (um[i][j]+um[i-1][j]) * areau_w[j] / epsilon1[i][j];

        df[i][j] = fe[i][j] - fw[i][j] + fn[i][j] - fs[i][j];

        dn[i][j] = (epsilon1[i][j]/re) * areau_n[i] / (y[j+1]-y[j]);
        ds[i][j] = (epsilon1[i][j]/re) * areau_s[i] / (y[j]-y[j-1]);
        de[i][j] = (epsilon1[i][j]/re) * areau_e[j] / (xm[i+1]-xm[i]);
        dw[i][j] = (epsilon1[i][j]/re) * areau_w[j] / (xm[i]-xm[i-1]);

        //upwind
        aw[i][j] = dw[i][j] + max(fw[i][j], 0.0);
        as[i][j] = ds[i][j] + max(fs[i][j], 0.0);
        ae[i][j] = de[i][j] + max(0.0, -fe[i][j]);
        an[i][j] = dn[i][j] + max(0.0, -fn[i][j]);

        ap[i][j] = aw[i][j] + ae[i][j] + as[i][j] + an[i][j] + df[i][j];

        u_w[i][j] = um[i-1][j];
        u_e[i][j] = um[i+1][j];
        u_s[i][j] = um[i][j-1];
        u_n[i][j] = um[i][j+1];
        u_p[i][j] = um[i][j];
        v_p[i][j] = vm[i][j];

        dudxdx[i][j] = areau_e[j] * (u_e[i][j]-u_p[i][j]) / (xm[i+1]-xm[i])
                     - areau_w[j] * (u_p[i][j]-u_w[i][j]) / (xm[i]-xm[i-1]); 
    
        dxdvdy[i][j] = areau_e[j] * (vm[i][j+1]-vm[i][j]) / (ym[j+1]-ym[j])
                     - areau_w[j] * (vm[i-1][j+1]-vm[i-1][j]) / (ym[j+1]-ym[j]);

        //bulk artificial viscosity term from Ramshaw(1990)
        q_art[i][j] = epsilon1[i][j] * (p[i][j]-p[i-1][j]) / (x[i]-x[i-1]) - iterations.b_art * (dudxdx[i][j]+dxdvdy[i][j]);

        ru[i][j] = 1.0 / (x[i]-x[i-1]) / (y[j]-y[j-1]) * (-ap[i][j]*u_p[i][j]
                +  aw[i][j] * u_w[i][j] + ae[i][j] * u_e[i][j]
                +  as[i][j] * u_s[i][j] + an[i][j] * u_n[i][j])
                -  q_art[i][j] - epsilon1[i][j] * (u_p[i][j]/(re*darcy_number) 
                +  cf/pow((epsilon1[i][j]*darcy_number),0.5) * u_p[i][j]
                *  (pow((pow(u_p[i][j],2.0) + pow(v_p[i][j],2.0)),0.5))) 
                *  liga_poros[i][j] - g * epsilon1[i][j];
    }
}

//--- upwind_U ---
void upwind_Uj(double **um, double **vm, double **p, double **ru, int i){
    int j;

    for(j = 2; j <= jmax-1; j++){
        //compute x-direction velocity component un
        fn[i][j] = 0.5 * (vm[i][j+1]+vm[i-1][j+1]) * areau_n[i] / epsilon1[i][j];
        fs[i][j] = 0.5 * (vm[i][j]+vm[i-1][j]) * areau_s[i] / epsilon1[i][j];
        fe[i][j] = 0.5 * (um[i+1][j]+um[i][j]) * areau_e[j] / epsilon1[i][j];
        fw[i][j] = 0.5 * (um[i][j]+um[i-1][j]) * areau_w[j] / epsilon1[i][j];

        df[i][j] = fe[i][j] - fw[i][j] + fn[i][j] - fs[i][j];

        dn[i][j] = (epsilon1[i][j]/re) * areau_n[i] / (y[j+1]-y[j]);
        ds[i][j] = (epsilon1[i][j]/re) * areau_s[i] / (y[j]-y[j-1]);
        de[i][j] = (epsilon1[i][j]/re) * areau_e[j] / (xm[i+1]-xm[i]);
        dw[i][j] = (epsilon1[i][j]/re) * areau_w[j] / (xm[i]-xm[i-1]);

        //upwind
        aw[i][j] = dw[i][j] + max(fw[i][j], 0.0);
        as[i][j] = ds[i][j] + max(fs[i][j], 0.0);
        ae[i][j] = de[i][j] + max(0.0, -fe[i][j]);
        an[i][j] = dn[i][j] + max(0.0, -fn[i][j]);

        ap[i][j] = aw[i][j] + ae[i][j] + as[i][j] + an[i][j] + df[i][j];

        u_w[i][j] = um[i-1][j];
        u_e[i][j] = um[i+1][j];
        u_s[i][j] = um[i][j-1];
        u_n[i][j] = um[i][j+1];
        u_p[i][j] = um[i][j];
        v_p[i][j] = vm[i][j];

        dudxdx[i][j] = areau_e[j] * (u_e[i][j]-u_p[i][j]) / (xm[i+1]-xm[i])
                     - areau_w[j] * (u_p[i][j]-u_w[i][j]) / (xm[i]-xm[i-1]); 
    
        dxdvdy[i][j] = areau_e[j] * (vm[i][j+1]-vm[i][j]) / (ym[j+1]-ym[j])
                     - areau_w[j] * (vm[i-1][j+1]-vm[i-1][j]) / (ym[j+1]-ym[j]);

        //bulk artificial viscosity term from Ramshaw(1990)
        q_art[i][j] = epsilon1[i][j] * (p[i][j]-p[i-1][j]) / (x[i]-x[i-1]) 
                    - (iterations.b_art) * (dudxdx[i][j]+dxdvdy[i][j]);

        ru[i][j] = 1.0 / (x[i]-x[i-1]) / (y[j]-y[j-1]) * (-ap[i][j]*u_p[i][j]
                +  aw[i][j] * u_w[i][j] + ae[i][j] * u_e[i][j]
                +  as[i][j] * u_s[i][j] + an[i][j] * u_n[i][j])
                -  q_art[i][j] - epsilon1[i][j] * (u_p[i][j]/(re*darcy_number) 
                +  cf/(pow((epsilon1[i][j]*darcy_number),0.5)) * u_p[i][j]
                *  (pow((pow(u_p[i][j],2.0) + pow(v_p[i][j],2.0)),0.5))) 
                *  liga_poros[i][j] - g * epsilon1[i][j];
    }
}

//--- upwind_V ---
void upwind_Vi(double **um, double **vm, double **p, double **rv, double **t, int j){
    int i, itc;

    for(i = 2; i <= imax-1; i++){
        fn[i][j] = 0.5 * (vm[i][j]+vm[i][j+1]) * areav_n[i] / epsilon1[i][j]; 
        fs[i][j] = 0.5 * (vm[i][j]+vm[i][j-1]) * areav_s[i] / epsilon1[i][j]; 
        fe[i][j] = 0.5 * (um[i+1][j]+um[i+1][j-1]) * areav_e[j] / epsilon1[i][j];  
        fw[i][j] = 0.5 * (um[i][j] + um[i][j-1]) * areav_w[j] / epsilon1[i][j];

        df[i][j] = fe[i][j] - fw[i][j] + fn[i][j] - fs[i][j];

        dn[i][j] = (epsilon1[i][j]/re) * areav_n[i] / (ym[j+1]-ym[j]);
        ds[i][j] = (epsilon1[i][j]/re) * areav_s[i] / (ym[j]-ym[j-1]);
        de[i][j] = (epsilon1[i][j]/re) * areav_e[j] / (x[i+1]-x[i]);
        dw[i][j] = (epsilon1[i][j]/re) * areav_w[j] / (x[i]-x[i-1]);

        //upwind
        aw[i][j] = dw[i][j] + max(fw[i][j], 0.0);
        as[i][j] = ds[i][j] + max(fs[i][j], 0.0);
        ae[i][j] = de[i][j] + max(0.0, -fe[i][j]);
        an[i][j] = dn[i][j] + max(0.0, -fn[i][j]);

        ap[i][j] = aw[i][j] + ae[i][j] + as[i][j] + an[i][j] + df[i][j];

        v_w[i][j] = vm[i-1][j];
        v_e[i][j] = vm[i+1][j];
        v_s[i][j] = vm[i][j-1];
        v_n[i][j] = vm[i][j+1];
        v_p[i][j] = vm[i][j];
        u_p[i][j] = um[i][j];

        dvdydy[i][j] = areav_n[i] * (v_n[i][j]-v_p[i][j]) / (ym[j+1]-ym[j])
                     - areav_s[i] * (v_p[i][j]-v_s[i][j]) / (ym[j]-ym[j-1]); 
            
        dydudx[i][j] = areav_n[i] * (um[i+1][j]-um[i][j]) / (xm[i+1]-xm[i])
                     - areav_s[i] * (um[i+1][j-1]-um[i][j-1]) / (xm[i+1]-xm[i]);
            
        //bulk artificial viscosity term from Ramshaw(1990)
        q_art[i][j] = epsilon1[i][j] * (p[i][j]-p[i][j-1]) / (y[j]-y[j-1]) - iterations.b_art * (dydudx[i][j]+dvdydy[i][j]);
        rv[i][j] = 1.0 / (x[i]-x[i-1]) / (y[j]-y[j-1]) * (-ap[i][j] * v_p[i][j]
                +  aw[i][j] * v_w[i][j] + ae[i][j] * v_e[i][j] 
                +  as[i][j] * v_s[i][j] + an[i][j] * v_n[i][j])  
                -  q_art[i][j] + invfr2 * (1.0 - 1.0 / ((t[i][j]+t[i][j-1]) * 0.5))
                -  epsilon1[i][j] * (v_p[i][j]/(re*darcy_number) 
                +  cf/(pow((epsilon1[i][j]*darcy_number), 0.5)) * v_p[i][j]
                *  (pow((pow(u_p[i][j], 2.0) + pow(v_p[i][j], 2.0)), 0.5))) * liga_poros[i][j]; 
    }
}

//--- upwind_V ---
void upwind_Vj(double **um, double **vm, double **p, double **rv, double **t, int i){
    int j;
    for(j=2; j <= jmax; j++){
        fn[i][j] = 0.5 * (vm[i][j] + vm[i][j+1]) * areav_n[i] / epsilon1[i][j]; 
        fs[i][j] = 0.5 * (vm[i][j] + vm[i][j-1]) * areav_s[i] / epsilon1[i][j]; 
        fe[i][j] = 0.5 * (um[i+1][j] + um[i+1][j-1]) * areav_e[j] / epsilon1[i][j];  
        fw[i][j] = 0.5 * (um[i][j] + um[i][j-1]) * areav_w[j] / epsilon1[i][j];  

        df[i][j] = fe[i][j] - fw[i][j] + fn[i][j] - fs[i][j];

        dn[i][j] = (epsilon1[i][j]/re) * areav_n[i] / (ym[j+1]-ym[j]);
        ds[i][j] = (epsilon1[i][j]/re) * areav_s[i] / (ym[j]-ym[j-1]);
        de[i][j] = (epsilon1[i][j]/re) * areav_e[j] / (x[i+1]-x[i]);
        dw[i][j] = (epsilon1[i][j]/re) * areav_w[j] / (x[i]-x[i-1]);

        //upwind
        aw[i][j] = dw[i][j] + max(fw[i][j], 0.0);
        as[i][j] = ds[i][j] + max(fs[i][j], 0.0);

        ae[i][j] = de[i][j] + max(0.0, -fe[i][j]);
        an[i][j] = dn[i][j] + max(0.0, -fn[i][j]);

        ap[i][j] = aw[i][j] + ae[i][j] + as[i][j] + an[i][j] + df[i][j];

        v_w[i][j] = vm[i-1][j];
        v_e[i][j] = vm[i+1][j];
        v_s[i][j] = vm[i][j-1];
        v_n[i][j] = vm[i][j+1];
        v_p[i][j] = vm[i][j];
        u_p[i][j] = um[i][j];

        dvdydy[i][j] = areav_n[i] * (v_n[i][j]-v_p[i][j]) / (ym[j+1]-ym[j])
                    -  areav_s[i] * (v_p[i][j]-v_s[i][j]) / (ym[j]-ym[j-1]); 
            
        dydudx[i][j] = areav_n[i] * (um[i+1][j]-um[i][j]) / (xm[i+1]-xm[i])
                    -  areav_s[i] * (um[i+1][j-1]-um[i][j-1]) / (xm[i+1]-xm[i]);
            
        //bulk artificial viscosity term from Ramshaw(1990)
        q_art[i][j] = epsilon1[i][j] * (p[i][j]-p[i][j-1]) / (y[j]-y[j-1]) - iterations.b_art * (dydudx[i][j]+dvdydy[i][j]);
        rv[i][j] = 1.0 / (x[i]-x[i-1]) / (y[j]-y[j-1]) * (-ap[i][j] * v_p[i][j]
                +  aw[i][j] * v_w[i][j] + ae[i][j] * v_e[i][j]
                +  as[i][j] * v_s[i][j] + an[i][j] * v_n[i][j])
                -  q_art[i][j] + invfr2 * (1.0 - 1.0 / ((t[i][j]+t[i][j-1]) * 0.5))
                -  epsilon1[i][j] * (v_p[i][j]/(re*darcy_number) 
                +  cf/(pow((epsilon1[i][j]*darcy_number),0.5)) * v_p[i][j] 
                *  (pow((pow(u_p[i][j], 2.0) + pow(v_p[i][j], 2.0)), 0.5)))*liga_poros[i][j];
    }
}

//resu////////////////////
void RESU(double **um, double **vm, double **p, double **ru){
    int i, j;   
    for(j = 3; j <= jmax-2; j++){
        for(i = 3; i <= imax-1; i++){
            fn[i][j] = 0.5 * (vm[i][j+1] + vm[i-1][j+1]) * areau_n[i] / epsilon1[i][j];
            fs[i][j] = 0.5 * (vm[i][j] + vm[i-1][j]) * areau_s[i] / epsilon1[i][j];
            fe[i][j] = 0.5 * (um[i+1][j] + um[i][j]) * areau_e[j] / epsilon1[i][j];
            fw[i][j] = 0.5 * (um[i][j] + um[i-1][j]) * areau_w[j] / epsilon1[i][j];
            df[i][j] = fe[i][j] - fw[i][j] + fn[i][j] - fs[i][j];

            dn[i][j] = (epsilon1[i][j]/re) * areau_n[i] / (y[j+1] - y[j]);
            ds[i][j] = (epsilon1[i][j]/re) * areau_s[i] / (y[j] - y[j-1]);
            de[i][j] = (epsilon1[i][j]/re) * areau_e[j] / (xm[i+1] - xm[i]);
            dw[i][j] = (epsilon1[i][j]/re) * areau_w[j] / (xm[i] - xm[i-1]);

        //quick
            if(fw[i][j] > 0.0){
            	afw[i][j] = 1.0;
			}else if(fw[i][j] < 0.0){
				afw[i][j] = 0.0;
			}

			if(fe[i][j] > 0.0){
				afe[i][j] = 1.0;
			}else if(fe[i][j] < 0.0){
				afe[i][j] = 0.0;
			}

			if(fn[i][j] > 0.0){
				afn[i][j] = 1.0;
			}else if(fn[i][j] < 0.0){
				afn[i][j] = 0.0;
			}

			if(fs[i][j] > 0.0){
				afs[i][j] = 1.0;
			}else if(fs[i][j] < 0.0){
				afs[i][j] = 0.0;
			}

			aw[i][j] = dw[i][j] + 0.75  * afw[i][j] * fw[i][j]
					+  0.125 * afe[i][j] * fe[i][j]
					+  0.375 * (1.0 - afw[i][j]) * fw[i][j];

			ae[i][j] = de[i][j] - 0.375 * afe[i][j] * fe[i][j]
					-  0.75  * (1.0 - afe[i][j]) * fe[i][j]
					-  0.125 * (1.0 - afw[i][j]) * fw[i][j];

			as[i][j] = ds[i][j] + 0.75 * afs[i][j] * fs[i][j]
					+  0.125 * afn[i][j] * fn[i][j]
					+  0.375 * (1.0 - afs[i][j]) * fs[i][j];

			an[i][j] = dn[i][j] - 0.375 * afn[i][j] * fn[i][j]
					-  0.75 * (1.0 - afn[i][j]) * fn[i][j]
					-  0.125 * (1.0 - afs[i][j]) * fs[i][j];

			aww[i][j] = -0.125 * afw[i][j] * fw[i][j];
			aee[i][j] =  0.125 * (1.0 - afe[i][j]) * fe[i][j];
			ass[i][j] = -0.125 * afs[i][j] * fs[i][j];
			ann[i][j] =  0.125 * (1.0 - afn[i][j]) * fn[i][j];
			ap[i][j] = aw[i][j] + ae[i][j] + as[i][j] + an[i][j] 
					+  aww[i][j] + aee[i][j] + ass[i][j] + ann[i][j] + df[i][j];
			//end Quick//////////////////////////////////////////////////////////////

			u_w[i][j]  = um[i-1][j];
			u_ww[i][j] = um[i-2][j];
			u_e[i][j]  = um[i+1][j];
			u_ee[i][j] = um[i+2][j];
			u_s[i][j]  = um[i][j-1];
			u_ss[i][j] = um[i][j-2];
			u_n[i][j]  = um[i][j+1];
			u_nn[i][j] = um[i][j+2];        
			u_p[i][j]  = um[i][j];
			v_p[i][j]  = vm[i][j];

			dudxdx[i][j] = areau_e[j] * (u_e[i][j] - u_p[i][j]) / (xm[i+1] - xm[i])
					    -  areau_w[j] * (u_p[i][j] - u_w[i][j]) / (xm[i] - xm[i-1]);
     
			dxdvdy[i][j] = areau_e[j] * (vm[i][j+1] - vm[i][j]) / (ym[j+1] - ym[j])
                        -  areau_w[j] * (vm[i-1][j+1] - vm[i-1][j]) / (ym[j+1] - ym[j]);
			
			artdivu[i][j] = -iterations.b_art * (dudxdx[i][j] + dxdvdy[i][j]);   
			
			//bulk artificial viscosity term from Ramshaw(1990)
			q_art[i][j] = epsilon1[i][j] * (p[i][j] - p[i-1][j]) / (x[i] - x[i-1]) + artdivu[i][j];

			ru[i][j] = 1.0 / (x[i]-x[i-1]) / (y[j] - y[j-1]) * (-ap[i][j] * u_p[i][j]
					+  aww[i][j] * u_ww[i][j] + aw[i][j] * u_w[i][j] 
					+  aee[i][j] * u_ee[i][j] + ae[i][j] * u_e[i][j]  
					+  ass[i][j] * u_ss[i][j] + as[i][j] * u_s[i][j]  
					+  ann[i][j] * u_nn[i][j] + an[i][j] * u_n[i][j]) 
					-  q_art[i][j] - epsilon1[i][j] * (u_p[i][j] / (re * darcy_number) 
                    +  cf / (pow((epsilon1[i][j] * darcy_number),0.5)) * u_p[i][j] 
                    *  (pow((pow(u_p[i][j],2.0) + pow(v_p[i][j],2.0)), 0.5))) * liga_poros[i][j] 
                    -  g * epsilon1[i][j];
		}
    }

    upwind_Ui(um,vm,p,ru,2);
    upwind_Ui(um,vm,p,ru,jmax-1);
    upwind_Uj(um,vm,p,ru,2);
    upwind_Uj(um,vm,p,ru,imax);
}

//--- ResV ---
void RESV(double **um, double **vm, double **p, double **t, double **rv){
    int i, j;

    for(j = 3; j <= jmax-1; j++){
        for(i = 3; i <= imax-2; i++){
            fn[i][j] = 0.5 * (vm[i][j]+vm[i][j+1]) * areav_n[i] / epsilon1[i][j];
            fs[i][j] = 0.5 * (vm[i][j]+vm[i][j-1]) * areav_s[i] / epsilon1[i][j];
            fe[i][j] = 0.5 * (um[i+1][j]+um[i+1][j-1]) * areav_e[j] / epsilon1[i][j];
            fw[i][j] = 0.5 * (um[i][j]+um[i][j-1]) * areav_w[j] / epsilon1[i][j];

            df[i][j] = fe[i][j] - fw[i][j] + fn[i][j] - fs[i][j];

            dn[i][j] = (epsilon1[i][j]/re) * areav_n[i] / (ym[j+1]-ym[j]);
            ds[i][j] = (epsilon1[i][j]/re) * areav_s[i] / (ym[j]-ym[j-1]);
            de[i][j] = (epsilon1[i][j]/re) * areav_e[j] / (x[i+1]-x[i]);
            dw[i][j] = (epsilon1[i][j]/re) * areav_w[j] / (x[i]-x[i-1]);

            //quick
            if(fw[i][j] > 0.0){
                afw[i][j] = 1.0;
            }else if(fw[i][j] < 0.0){
                afw[i][j] = 0.0;
            }

            if(fe[i][j] > 0.0){
                afe[i][j] = 1.0;
            }else if(fe[i][j] < 0.0){
                afe[i][j] = 0.0;
            }

            if(fn[i][j] > 0.0){
                afn[i][j] = 1.0;
            }else if(fn[i][j] < 0.0){
                afn[i][j] = 0.0;
            }

            if(fs[i][j] > 0.0){
                afs[i][j] = 1.0;
            }else if(fs[i][j] < 0.0){
                afs[i][j] = 0.0;
            }

            aw[i][j] = dw[i][j] + 0.75 * afw[i][j] * fw[i][j]
                    +  0.125 * afe[i][j] * fe[i][j]
                    +  0.375 * (1.0-afw[i][j]) * fw[i][j];

            ae[i][j] = de[i][j] - 0.375 * afe[i][j] * fe[i][j]
                    - 0.75 * (1.0-afe[i][j]) * fe[i][j]
                    - 0.125 * (1.0-afw[i][j]) * fw[i][j];

            as[i][j] = ds[i][j] + 0.75  * afs[i][j] * fs[i][j] 
                    +  0.125 * afn[i][j] * fn[i][j] 
                    +  0.375 * (1.0-afs[i][j]) * fs[i][j];

            an[i][j] = dn[i][j] - 0.375 * afn[i][j] * fn[i][j]
                    - 0.75 * (1.0-afn[i][j]) * fn[i][j]
                    - 0.125 * (1.0-afs[i][j]) * fs[i][j];

            aww[i][j] = -0.125 * afw[i][j] * fw[i][j];
            aee[i][j] =  0.125 * (1.0-afe[i][j]) * fe[i][j];
            ass[i][j] = -0.125 * afs[i][j] * fs[i][j];
            ann[i][j] =  0.125 * (1.0-afn[i][j]) * fn[i][j];

            ap[i][j] = aw[i][j] + ae[i][j] + as[i][j] + an[i][j] 
                    +  aww[i][j] + aee[i][j] + ass[i][j] + ann[i][j] + df[i][j];
            //end Quick///////////////////////////////////////////////

            v_w[i][j]  = vm[i-1][j];
            v_ww[i][j] = vm[i-2][j];
            v_e[i][j]  = vm[i+1][j];
            v_ee[i][j] = vm[i+2][j];
            v_s[i][j]  = vm[i][j-1];
            v_ss[i][j] = vm[i][j-2];
            v_n[i][j]  = vm[i][j+1];
            v_nn[i][j] = vm[i][j+2];         
            v_p[i][j]  = vm[i][j];
            u_p[i][j]  = um[i][j];

            dvdydy[i][j] = areav_n[i] * (v_n[i][j]-v_p[i][j]) / (ym[j+1]-ym[j])
                        -  areav_s[i] * (v_p[i][j]-v_s[i][j]) / (ym[j]-ym[j-1]);

            dydudx[i][j] = areav_n[i] * (um[i+1][j]-um[i][j]) / (xm[i+1]-xm[i])
                        -  areav_s[i] * (um[i+1][j-1]-um[i][j-1]) / (xm[i+1]-xm[i]);

            artdivv[i][j] = -(iterations.b_art) * (dydudx[i][j]+dvdydy[i][j]);            

            //bulk artificial viscosity term from Ramshaw(1990)
            q_art[i][j] = epsilon1[i][j] * (p[i][j]-p[i][j-1]) / (y[j]-y[j-1]) + artdivv[i][j];

            rv[i][j] = 1.0 / (x[i]-x[i-1]) / (y[j]-y[j-1]) 
                    * (-ap[i][j] * v_p[i][j] + aww[i][j] * v_ww[i][j] + aw[i][j] 
                    * v_w[i][j] + aee[i][j] * v_ee[i][j] + ae[i][j] * v_e[i][j] 
                    + ass[i][j] * v_ss[i][j] + as[i][j] * v_s[i][j] 
                    + ann[i][j] * v_nn[i][j] + an[i][j] * v_n[i][j]) 
                    - q_art[i][j] + invfr2 * (1.0 - 1.0 / ((t[i][j]+t[i][j-1]) * 0.5))
                    - epsilon1[i][j]*(v_p[i][j]/(re*darcy_number) 
                    + cf/(pow((epsilon1[i][j]*darcy_number), 0.5)) * v_p[i][j] 
                    * (pow((pow(u_p[i][j], 2.0) + pow(v_p[i][j], 2.0)), 0.5))) * liga_poros[i][j]; 
        }
    }

    upwind_Vi(um,vm,p,rv,t,2);
    upwind_Vi(um,vm,p,rv,t,jmax);    
    upwind_Vj(um,vm,p,rv,t,2);
    upwind_Vj(um,vm,p,rv,t,imax-1);
}

//--- ResZ ---
void RESZ(double **um_n, double **vm_n, double **z, double **rz){
    int i, j;

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            dzudx[i][j] = 0.5 * (z[i+1][j]+z[i][j]) * um_n[i+1][j] * areau_e[j]
                        - 0.5 * (z[i-1][j]+z[i][j]) * um_n[i][j] * areau_w[j];
            dzvdy[i][j] = 0.5 * (z[i][j+1]+z[i][j]) * vm_n[i][j+1] * areav_n[i]
                        - 0.5 * (z[i][j-1]+z[i][j]) * vm_n[i][j] * areav_s[i];

            de[i][j] = (ym[j+1]-ym[j]) * (1.0/pe) / (x[i+1]-x[i]);  
            dw[i][j] = (ym[j+1]-ym[j]) * (1.0/pe) / (x[i]-x[i-1]); 
            dn[i][j] = (xm[i+1]-xm[i]) * (1.0/pe) / (y[j+1]-y[j]);  
            ds[i][j] = (xm[i+1]-xm[i]) * (1.0/pe) / (y[j]-y[j-1]);  
            dp[i][j] = de[i][j] + dw[i][j] + dn[i][j] + ds[i][j];

            rz[i][j] = 1.0 / (xm[i+1]-xm[i]) / (ym[j+1]-ym[j]) 
                    *  (-dp[i][j]*z[i][j] + de[i][j]*z[i+1][j] 
                    +  dw[i][j]*z[i-1][j] + dn[i][j]*z[i][j+1] 
                    +  ds[i][j]*z[i][j-1] - (1.0-liga_poros[i][j])
                    *  (dzudx[i][j]+dzvdy[i][j])) / (liga_poros[i][j]
                    *  (epsilon1[i][j]-1.0)+1.0);
        }
    }
}

//--- ResC ---
void RESC(double **um_n, double **vm_n, double **c, double **rc){
    int i, j;

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            dcudx[i][j] = 0.5 * (c[i+1][j]+c[i][j]) * um_n[i+1][j] * areau_e[j]
                        - 0.5 * (c[i-1][j]+c[i][j]) * um_n[i][j] * areau_w[j];
            dcvdy[i][j] = 0.5 * (c[i][j+1]+c[i][j]) * vm_n[i][j+1] * areav_n[i]
                        - 0.5 * (c[i][j-1]+c[i][j]) * vm_n[i][j] * areav_s[i];

            de[i][j] = (ym[j+1]-ym[j]) * (1.0/re/sc) / (x[i+1]-x[i]);
            dw[i][j] = (ym[j+1]-ym[j]) * (1.0/re/sc) / (x[i]-x[i-1]);
            dn[i][j] = (xm[i+1]-xm[i]) * (1.0/re/sc) / (y[j+1]-y[j]);
            ds[i][j] = (xm[i+1]-xm[i]) * (1.0/re/sc) / (y[j]-y[j-1]);
            dp[i][j] = de[i][j] + dw[i][j] + dn[i][j] + ds[i][j];

            rc[i][j] = 1.0 / (xm[i+1]-xm[i]) / (ym[j+1]-ym[j]) 
                    *  (-dp[i][j]*c[i][j] + de[i][j]*c[i+1][j] 
                    +  dw[i][j]*c[i-1][j] + dn[i][j]*c[i][j+1] 
                    +  ds[i][j]*c[i][j-1] - (dcudx[i][j] + dcvdy[i][j]) 
                    /  (liga_poros[i][j]*(epsilon1[i][j]-1.0)+1.0));
        }
    }
}

//--- solve_U ---
void solve_U(double **um, double **vm, double **um_n, double **um_tau, double **vm_tau, double **um_n_tau, double **p, double *residual_u){
    int i, j;
    double **ui, **ru, **res_u;

    ui = (double**)malloc(sizeof(double*)*(imax+2));
    ru = (double**)malloc(sizeof(double*)*(imax+2));
    for(i = 0; i < (imax+2); i++){
        ui[i] = (double*)malloc(sizeof(double)*(jmax+1));
        ru[i] = (double*)malloc(sizeof(double)*(jmax+1));
    }

    res_u = (double**)malloc(sizeof(double*)*(imax+2));
    for(i = 0; i < (imax+2); i++){
        res_u[i] = (double*)malloc(sizeof(double)*(jmax+1));
    }

    RESU(um_tau,vm_tau,p,ru);
    
    for(j = 2; j <= jmax-1; j++){
        for(i = 3; i <= imax-1; i++){
            res_u[i][j] = ((um[i][j]-um_tau[i][j]) + ru[i][j]*dt) * dtau;
            ui[i][j] = um_tau[i][j] + res_u[i][j];
        }
    }
    
    bcUV(ui,vm_tau);

    RESU(ui,vm_tau,p,ru);
    
    for(j = 2; j <= jmax-1; j++){
        for(i = 3; i <= imax-1; i++){
            res_u[i][j] = ((um[i][j]-um_tau[i][j]) + ru[i][j]*dt) * dtau;
            ui[i][j] = 0.75 * um_tau[i][j] + 0.25 * (ui[i][j]+res_u[i][j]);            
        }
    }
    
    bcUV(ui,vm_tau);
    RESU(ui,vm_tau,p,ru);
    
    for(j = 2; j <= jmax-1; j++){
        for(i = 3; i <= imax-1; i++){
            res_u[i][j] = ((um[i][j]-um_tau[i][j]) + ru[i][j]*dt) * dtau;
            um_n_tau[i][j] = 1.0 / 3.0 * um_tau[i][j] + 2.0 / 3.0 * (ui[i][j]+res_u[i][j]); 
        }
    }
    
    bcUV(um_n_tau, vm_tau);

    (*residual_u) =  maior_valor(res_u, imax+2, jmax+1, 3, 2); 

    //desalocando
    for(i = 0; i < (imax+2); i++){
        free(ui[i]);
        free(ru[i]);
    }
    free(ui);
    free(ru);

    for(i = 0; i < (imax+2); i++){
        free(res_u[i]);
    }
    free(res_u);
}

//--- solve_V ---
void solve_V(double **um, double **vm, double **vm_n, double **um_tau, double **vm_tau, double **vm_n_tau, double **p, double **T, double *residual_v){
    int i, j;

    double **vi, **RV, **res_v;
    vi = (double**)malloc(sizeof(double*)*(imax+1));
    RV = (double**)malloc(sizeof(double*)*(imax+1));
    res_v = (double**)malloc(sizeof(double*)*(imax+1));
    for(i = 0; i < imax+1; i++){
        vi[i] = (double*)malloc(sizeof(double)*(jmax+2));
        RV[i] = (double*)malloc(sizeof(double)*(jmax+2));
        res_v[i] = (double*)malloc(sizeof(double)*(jmax+2));
    }
    
    /*for(i = 1; i <= imax; i++){
        for(j = 1; j <= jmax; j++){                
            printf("vm_n_tau [%d][%d] = %.14lf\n", i, j, vm_n_tau[i][j]);
            }
            }*/
           
    RESV(um_tau,vm_tau,p,T,RV);
    
    for(j = 3; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_v[i][j] = ((vm[i][j]-vm_tau[i][j]) + RV[i][j]*dt) * dtau;
            vi[i][j] = (vm_tau[i][j] + res_v[i][j]);
        }
    }
    
    for(int i = 1; i <= imax; i++){
        for(int j = 1; j <= jmax+1; j++){
            printf("[%i][%i] res_v = %lf\n", i, j, res_v[i][j]);
        }
    }

    bcUV(um_tau,vi);
    RESV(um_tau,vi,p,T,RV);

    for(j = 3; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_v[i][j] = ((vm[i][j]-vm_tau[i][j]) + RV[i][j]*dt) * dtau;
            vi[i][j] = ((double)(3.0/4.0) * vm_tau[i][j] + (double)(1.0/4.0) * (vi[i][j] + res_v[i][j]));
        }   
    }

    bcUV(um_tau,vi);
    RESV(um_tau,vi,p,T,RV);

    for(j = 3; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_v[i][j] =((vm[i][j]-vm_tau[i][j]) +  RV[i][j]*dt) * dtau;
            vm_n_tau[i][j] = (double)(1.0 / 3.0) * vm_tau[i][j] + (double)(2.0 / 3.0) * (vi[i][j] + res_v[i][j]); 
        }
    }


    bcUV(um_tau,vm_n_tau);
    (*residual_v) = maior_valor(res_v, imax+1, jmax+2, 2, 3); 

    //desalocando
    for(i = 0; i < imax+1; i++){
        free(vi[i]);
        free(RV[i]);
        free(res_v[i]);
    }
    free(vi);
    free(RV);
    free(res_v);
}

//--- Solve Continuity Equation ---
void solve_P(double **p, double **um_n, double **vm_n, double **pn, double *residual_p){
    int i, j;

    // RALSTON'S METHOD (Second Order Runge-Kutta) mass conservation
    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            dudx[i][j] = um_n[i+1][j] * areau_e[j] - um_n[i][j] * areau_w[j];
            dvdy[i][j] = vm_n[i][j+1] * areav_n[i] - vm_n[i][j] * areav_s[i];
            rp[i][j] = - (dudx[i][j]+dvdy[i][j]); 
            pi[i][j] = p[i][j] + dtau * rp[i][j] * iterations.beta;
        }
    }

    bcP(pi);

    for(j = 3; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){ 
            pi[i][j] = (double)(3.0/4.0) * p[i][j] + (double)(1.0/4.0) * (pi[i][j] + dtau * rp[i][j] * iterations.beta);
        }
    }

    bcP(pi);

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_p[i][j] = dtau * rp[i][j] * iterations.beta;
            pn[i][j] = (double)(1.0 / 3.0) * p[i][j] + (double)(2.0 / 3.0) * (pi[i][j] + res_p[i][j]);
        }
    }


    bcP(pn);
            //for(i = 1; i <= imax; i++){
            //    for(j = 1; j <= jmax; j++){
            //        printf("pn [%d][%d] = %lf\n", i, j, pn[i][j]);
            //    }
            //}
    (*residual_p) = maior_valor(res_p, imax+1, jmax+1, 1, 1);
}

//--- Solve Mixture Fraction ---
void solve_Z(double **um_n, double **vm_n, double **z, double **z_n_tau, double **z_tau){
    int i, j;

    // RALSTON'S METHOD (Second Order Runge-Kutta)
    RESZ(um_n,vm_n, z_tau, rz);

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){    
            res_z[i][j] = ((z[i][j]-z_tau[i][j]) + rz[i][j]*dt) * dtau; 
            zi[i][j] = z_tau[i][j] + res_z[i][j];             
        }
    }

    bcZ(zi);
    RESZ(um_n,vm_n,zi,rz);

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){    
            res_z[i][j] = ((z[i][j]-z_tau[i][j]) + rz[i][j]*dt) * dtau;
            zi[i][j] = 0.75 * z_tau[i][j] + 0.25 * (zi[i][j]+res_z[i][j]);
        }
    }
    
    bcZ(zi);
    RESZ(um_n,vm_n,zi,rz);

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_z[i][j] = ((z[i][j]-z_tau[i][j]) + rz[i][j]*dt) * dtau; 
            z_n_tau[i][j] = 1.0 / 3.0 * z_tau[i][j] + 2.0 / 3.0 * (zi[i][j]+res_z[i][j]);
        }
    }

    bcZ(z_n_tau);
}

//--- solve_C - concentration ---
void solve_C(double **um_n, double **vm_n, double **c, double **c_n_tau, double **c_tau){
    int i, j;

    // RALSTON'S METHOD (Second Order Runge-Kutta)
    RESC(um_n, vm_n, c_tau, rc);

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_c[i][j] = ((c[i][j]-c_tau[i][j]) + rc[i][j]*dt) * dtau; 
            ci[i][j] = c_tau[i][j] + res_c[i][j];
        }
    }

    bcC(ci);
    RESC(um_n,vm_n,ci,rc);
    
    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_c[i][j] = ((c[i][j]-c_tau[i][j]) + rc[i][j]*dt) * dtau; 
            ci[i][j] = 0.75 * c_tau[i][j] + 0.25 * (ci[i][j] + res_c[i][j]);   
        }
    }

    bcC(ci);
    RESC(um_n, vm_n, ci, rc);

    for(j = 2; j <= jmax-1; j++){
        for(i = 2; i <= imax-1; i++){
            res_c[i][j] = ((c[i][j]-c_tau[i][j]) + rc[i][j]*dt) * dtau; 
            c_n_tau[i][j] = 1.0 / 3.0 * c_tau[i][j] + 2.0 / 3.0 * (ci[i][j] + res_c[i][j]);            
        }
    }

    bcC(c_n_tau);
}