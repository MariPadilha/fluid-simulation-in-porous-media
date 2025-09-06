#include "comum.h"
#define idx 1

__global__ void calc_resu_parte1(){
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
        }
    }
}

__global__ void calc_resu_parte2(){

    for(j = 3; j <= jmax-2; j++){
        for(i = 3; i <= imax-1; i++){

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
}

//resu////////////
void RESU(double **um, double **vm, double **p, double **ru){
    int i, j;   

    calc_resu_parte1<<<thrads, blocks>>>>

    for(j = 3; j <= jmax-2; j++){
        for(i = 3; i <= imax-1; i++){

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

        }
    }

    calc_resu_parte2<<<thread, blocks>>>;

    upwind_Ui(um,vm,p,ru,2);
    upwind_Ui(um,vm,p,ru,jmax-1);
    upwind_Uj(um,vm,p,ru,2);
    upwind_Uj(um,vm,p,ru,imax);
}