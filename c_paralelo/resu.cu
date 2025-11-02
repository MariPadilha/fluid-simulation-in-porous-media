#include "comum.h"

__global__ void calc_resu(
	double *dev_areau_e, double *dev_areau_n, double *dev_areau_s, double *dev_areau_w, double *dev_epsilon1, double *dev_y, double *dev_x,
	double *dev_xm, double *dev_ym, double *dev_liga_poros, double b_art, int imax, int jmax, double re, double darcy_number, double g, 
	double cf, double *dev_um, double *dev_vm, double *dev_p, double *dev_ru){

    int i = blockIdx.x * blockDim.x + threadIdx.x + 3;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
	int idx = i*(jmax+1)+j;
	double df, dn, ds, de, dw;
	double fn, fs, fe, fw;
	double afw, afe, afn, afs;
	double aw, ae, as, an, ap;
	double aww, aee, ass, ann;
	double u_w, u_e, u_s, u_n, u_p, v_p;
	double u_ww, u_ee, u_ss, u_nn;
	double dudxdx, dxdvdy;
	double artdivu, q_art;
	double aux = dev_epsilon1[idx]/re;
	
	if(i <= imax-1 && j <= jmax-2){
		fn = 0.5 * (dev_vm[i*(jmax+2)+(j+1)] + dev_vm[(i-1)*(jmax+2)+(j+1)]) * dev_areau_n[i] / dev_epsilon1[idx];
		fs = 0.5 * (dev_vm[i*(jmax+2)+j] + dev_vm[(i-1)*(jmax+2)+j]) * dev_areau_s[i] / dev_epsilon1[idx];
		fe = 0.5 * (dev_um[(i+1)*(jmax+1)+j] + dev_um[idx]) * dev_areau_e[j] / dev_epsilon1[idx];
		fw = 0.5 * (dev_um[idx] + dev_um[(i-1)*(jmax+1)+j]) * dev_areau_w[j] / dev_epsilon1[idx];
		df = fe - fw + fn - fs;
		dn = aux * dev_areau_n[i] / (dev_y[j+1] - dev_y[j]);
		ds = aux * dev_areau_s[i] / (dev_y[j] - dev_y[j-1]);
		de = aux * dev_areau_e[j] / (dev_xm[i+1] - dev_xm[i]);
		dw = aux * dev_areau_w[j] / (dev_xm[i] - dev_xm[i-1]);

        //quick
		afw = (double)(fw > 0.0);
		afe = (double)(fe > 0.0);
		afn = (double)(fn > 0.0);
		afs = (double)(fs > 0.0);

        aw = dw + 0.75  * afw * fw
					+  0.125 * afe * fe
					+  0.375 * (1.0 - afw) * fw;

		ae = de - 0.375 * afe * fe
				-  0.75  * (1.0 - afe) * fe
				-  0.125 * (1.0 - afw) * fw;

		as = ds + 0.75 * afs * fs
				+  0.125 * afn * fn
				+  0.375 * (1.0 - afs) * fs;

		an = dn - 0.375 * afn * fn
				-  0.75 * (1.0 - afn) * fn
				-  0.125 * (1.0 - afs) * fs;

		aww = -0.125 * afw * fw;
		aee =  0.125 * (1.0 - afe) * fe;
		ass = -0.125 * afs * fs;
		ann =  0.125 * (1.0 - afn) * fn;
		ap = aw + ae + as + an + aww + aee + ass + ann + df;
		//end Quick//////////////////////////////////////////////////////////////

		u_w  = dev_um[(i-1)*(jmax+1)+j];
		u_ww = dev_um[(i-2)*(jmax+1)+j];
		u_e  = dev_um[(i+1)*(jmax+1)+j];
		u_ee = dev_um[(i+2)*(jmax+1)+j];
		u_s  = dev_um[i*(jmax+1)+(j-1)];
		u_ss = dev_um[i*(jmax+1)+(j-2)];
		u_n  = dev_um[i*(jmax+1)+(j+1)];
		u_nn = dev_um[i*(jmax+1)+(j+2)];        
		u_p  = dev_um[idx];
		v_p  = dev_vm[i*(jmax+2)+j];

		dudxdx = dev_areau_e[j] * (u_e - u_p) / (dev_xm[i+1] - dev_xm[i])
					-  dev_areau_w[j] * (u_p - u_w) / (dev_xm[i] - dev_xm[i-1]);
	
		dxdvdy = dev_areau_e[j] * (dev_vm[i*(jmax+2)+(j+1)] - dev_vm[i*(jmax+2)+j]) / (dev_ym[j+1] - dev_ym[j])
					-  dev_areau_w[j] * (dev_vm[(i-1)*(jmax+2)+(j+1)] - dev_vm[(i-1)*(jmax+2)+j]) / (dev_ym[j+1] - dev_ym[j]);

		artdivu = -b_art * (dudxdx + dxdvdy);

		//bulk artificial viscosity term from Ramshaw(1990)
		q_art = dev_epsilon1[idx] * (dev_p[idx] - dev_p[(i-1)*(jmax+1)+j]) / (dev_x[i] - dev_x[i-1]) + artdivu;

		dev_ru[idx] = 1.0 / (dev_x[i]-dev_x[i-1]) / (dev_y[j] - dev_y[j-1]) * (-ap * u_p
				+  aww * u_ww + aw * u_w
				+  aee * u_ee + ae * u_e
				+  ass * u_ss + as * u_s
				+  ann * u_nn + an * u_n)
				-  q_art - dev_epsilon1[idx] * (u_p / (re * darcy_number)
				+  cf / (pow((dev_epsilon1[idx] * darcy_number),0.5)) * u_p
				*  (pow((pow(u_p,2.0) + pow(v_p,2.0)), 0.5))) * dev_liga_poros[idx]
				-  g * dev_epsilon1[idx];
    }
}

//resu////////////
void RESU(double *dev_um, double *dev_vm, double *dev_p, double *dev_ru){
    dim3 blockDim(16,16);
    dim3 gridDim((imax-4 + blockDim.x - 1)/blockDim.x, (jmax-5 + blockDim.y - 1)/blockDim.y);

    calc_resu<<<gridDim, blockDim>>>(
	dev_areau_e, dev_areau_n, dev_areau_s, dev_areau_w, dev_epsilon1, dev_y, dev_x,
	dev_xm, dev_ym, dev_liga_poros, iterations.b_art, imax, jmax, re, darcy_number, g, cf,
	dev_um, dev_vm, dev_p, dev_ru);

    upwind_Ui(dev_um, dev_vm, dev_p, dev_ru, 2);
    upwind_Ui(dev_um, dev_vm, dev_p, dev_ru, jmax-1);
    upwind_Uj(dev_um, dev_vm, dev_p, dev_ru, 2);
    upwind_Uj(dev_um, dev_vm, dev_p, dev_ru, imax);
}