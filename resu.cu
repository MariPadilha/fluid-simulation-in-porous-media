#include "comum.h"
#define idx i*(jmax+1)+j

__global__ void calc_resu(double *dev_fn, double *dev_fs, double *dev_fe, double *dev_fw, double *dev_df, double *dev_dn, double *dev_ds, double *dev_de, 
	double *dev_dw, double *dev_areau_e, double *dev_areau_n, double *dev_areau_s, double *dev_areau_w, double *dev_epsilon1, double *dev_y, double *dev_x,
	double *dev_xm, double *dev_afw, double *dev_afe, double *dev_afn, double *dev_afs, double *dev_aw, double *dev_ae, double *dev_as, double *dev_an,
    double *dev_aww, double *dev_aee, double *dev_ass, double *dev_ann, double *dev_ap, double *dev_u_e, double *dev_u_ee, double *dev_u_n, double *dev_u_nn, 
	double *dev_u_p, double *dev_u_s, double *dev_u_ss, double *dev_u_w, double *dev_u_ww, double *dev_v_p, double *dev_ym, double *dev_dudxdx, double *dev_dxdvdy,
	double *dev_q_art, double *dev_artdivu, double *dev_liga_poros, double b_art, int imax, int jmax, double re, double darcy_number, double g, double cf,
	double *dev_um, double *dev_vm, double *dev_p, double *dev_ru){

    int i = blockIdx.x * blockDim.x + threadIdx.x + 3;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
	if(i <= imax-1 && j <= jmax-2){
		dev_fn[idx] = 0.5 * (dev_vm[i*(jmax+2)+(j+1)] + dev_vm[(i-1)*(jmax+2)+(j+1)]) * dev_areau_n[i] / dev_epsilon1[idx];
		dev_fs[idx] = 0.5 * (dev_vm[i*(jmax+2)+j] + dev_vm[(i-1)*(jmax+2)+j]) * dev_areau_s[i] / dev_epsilon1[idx];
		dev_fe[idx] = 0.5 * (dev_um[(i+1)*(jmax+1)+j] + dev_um[idx]) * dev_areau_e[j] / dev_epsilon1[idx];
		dev_fw[idx] = 0.5 * (dev_um[idx] + dev_um[(i-1)*(jmax+1)+j]) * dev_areau_w[j] / dev_epsilon1[idx];
		dev_df[idx] = dev_fe[idx] - dev_fw[idx] + dev_fn[idx] - dev_fs[idx];
		dev_dn[idx] = (dev_epsilon1[idx]/re) * dev_areau_n[i] / (dev_y[j+1] - dev_y[j]);
		dev_ds[idx] = (dev_epsilon1[idx]/re) * dev_areau_s[i] / (dev_y[j] - dev_y[j-1]);
		dev_de[idx] = (dev_epsilon1[idx]/re) * dev_areau_e[j] / (dev_xm[i+1] - dev_xm[i]);
		dev_dw[idx] = (dev_epsilon1[idx]/re) * dev_areau_w[j] / (dev_xm[i] - dev_xm[i-1]);

        //quick
		dev_afw[idx] = (double)(dev_fw[idx] > 0.0);
		dev_afe[idx] = (double)(dev_fe[idx] > 0.0);
		dev_afn[idx] = (double)(dev_fn[idx] > 0.0);
		dev_afs[idx] = (double)(dev_fs[idx] > 0.0);
        dev_aw[idx] = dev_dw[idx] + 0.75  * dev_afw[idx] * dev_fw[idx]
					+  0.125 * dev_afe[idx] * dev_fe[idx]
					+  0.375 * (1.0 - dev_afw[idx]) * dev_fw[idx];

		dev_ae[idx] = dev_de[idx] - 0.375 * dev_afe[idx] * dev_fe[idx]
				-  0.75  * (1.0 - dev_afe[idx]) * dev_fe[idx]
				-  0.125 * (1.0 - dev_afw[idx]) * dev_fw[idx];

		dev_as[idx] = dev_ds[idx] + 0.75 * dev_afs[idx] * dev_fs[idx]
				+  0.125 * dev_afn[idx] * dev_fn[idx]
				+  0.375 * (1.0 - dev_afs[idx]) * dev_fs[idx];

		dev_an[idx] = dev_dn[idx] - 0.375 * dev_afn[idx] * dev_fn[idx]
				-  0.75 * (1.0 - dev_afn[idx]) * dev_fn[idx]
				-  0.125 * (1.0 - dev_afs[idx]) * dev_fs[idx];

		dev_aww[idx] = -0.125 * dev_afw[idx] * dev_fw[idx];
		dev_aee[idx] =  0.125 * (1.0 - dev_afe[idx]) * dev_fe[idx];
		dev_ass[idx] = -0.125 * dev_afs[idx] * dev_fs[idx];
		dev_ann[idx] =  0.125 * (1.0 - dev_afn[idx]) * dev_fn[idx];
		dev_ap[idx] = dev_aw[idx] + dev_ae[idx] + dev_as[idx] + dev_an[idx] 
				+  dev_aww[idx] + dev_aee[idx] + dev_ass[idx] + dev_ann[idx] + dev_df[idx];
		//end Quick//////////////////////////////////////////////////////////////

		dev_u_w[idx]  = dev_um[(i-1)*(jmax+1)+j];
		dev_u_ww[idx] = dev_um[(i-2)*(jmax+1)+j];
		dev_u_e[idx]  = dev_um[(i+1)*(jmax+1)+j];
		dev_u_ee[idx] = dev_um[(i+2)*(jmax+1)+j];
		dev_u_s[idx]  = dev_um[i*(jmax+1)+(j-1)];
		dev_u_ss[idx] = dev_um[i*(jmax+1)+(j-2)];
		dev_u_n[idx]  = dev_um[i*(jmax+1)+(j+1)];
		dev_u_nn[idx] = dev_um[i*(jmax+1)+(j+2)];        
		dev_u_p[idx]  = dev_um[idx];
		dev_v_p[idx]  = dev_vm[i*(jmax+2)+j];

		dev_dudxdx[idx] = dev_areau_e[j] * (dev_u_e[idx] - dev_u_p[idx]) / (dev_xm[i+1] - dev_xm[i])
					-  dev_areau_w[j] * (dev_u_p[idx] - dev_u_w[idx]) / (dev_xm[i] - dev_xm[i-1]);
	
		dev_dxdvdy[idx] = dev_areau_e[j] * (dev_vm[i*(jmax+2)+(j+1)] - dev_vm[i*(jmax+2)+j]) / (dev_ym[j+1] - dev_ym[j])
					-  dev_areau_w[j] * (dev_vm[(i-1)*(jmax+2)+(j+1)] - dev_vm[(i-1)*(jmax+2)+j]) / (dev_ym[j+1] - dev_ym[j]);
		
		dev_artdivu[idx] = -b_art * (dev_dudxdx[idx] + dev_dxdvdy[idx]);   
		
		//bulk artificial viscosity term from Ramshaw(1990)
		dev_q_art[idx] = dev_epsilon1[idx] * (dev_p[idx] - dev_p[(i-1)*(jmax+1)+j]) / (dev_x[i] - dev_x[i-1]) + dev_artdivu[idx];

		dev_ru[idx] = 1.0 / (dev_x[i]-dev_x[i-1]) / (dev_y[j] - dev_y[j-1]) * (-dev_ap[idx] * dev_u_p[idx]
				+  dev_aww[idx] * dev_u_ww[idx] + dev_aw[idx] * dev_u_w[idx] 
				+  dev_aee[idx] * dev_u_ee[idx] + dev_ae[idx] * dev_u_e[idx]  
				+  dev_ass[idx] * dev_u_ss[idx] + dev_as[idx] * dev_u_s[idx]  
				+  dev_ann[idx] * dev_u_nn[idx] + dev_an[idx] * dev_u_n[idx]) 
				-  dev_q_art[idx] - dev_epsilon1[idx] * (dev_u_p[idx] / (re * darcy_number) 
				+  cf / (pow((dev_epsilon1[idx] * darcy_number),0.5)) * dev_u_p[idx] 
				*  (pow((pow(dev_u_p[idx],2.0) + pow(dev_v_p[idx],2.0)), 0.5))) * dev_liga_poros[idx] 
				-  g * dev_epsilon1[idx];
    }
}

//resu////////////
void RESU(double *dev_um, double *dev_vm, double *dev_p, double *dev_ru){
    dim3 blockDim(16,16);
    dim3 gridDim((imax-4 + blockDim.x - 1)/blockDim.x, (jmax-5 + blockDim.y - 1)/blockDim.y);

    calc_resu<<<gridDim, blockDim>>>(dev_fn, dev_fs, dev_fe, dev_fw, dev_df, dev_dn, dev_ds, dev_de, 
	dev_dw, dev_areau_e, dev_areau_n, dev_areau_s, dev_areau_w, dev_epsilon1, dev_y, dev_x,
	dev_xm, dev_afw, dev_afe, dev_afn, dev_afs, dev_aw, dev_ae, dev_as, dev_an,
    dev_aww, dev_aee, dev_ass, dev_ann, dev_ap, dev_u_e, dev_u_ee, dev_u_n, dev_u_nn, 
	dev_u_p, dev_u_s, dev_u_ss, dev_u_w, dev_u_ww, dev_v_p, dev_ym, dev_dudxdx, dev_dxdvdy,
	dev_q_art, dev_artdivu, dev_liga_poros, iterations.b_art, imax, jmax, re, darcy_number, g, cf,
	dev_um, dev_vm, dev_p, dev_ru);

    upwind_Ui(dev_um, dev_vm, dev_p, dev_ru, 2);
    upwind_Ui(dev_um, dev_vm, dev_p, dev_ru, jmax-1);
    upwind_Uj(dev_um, dev_vm, dev_p, dev_ru, 2);
    upwind_Uj(dev_um, dev_vm, dev_p, dev_ru, imax);
}