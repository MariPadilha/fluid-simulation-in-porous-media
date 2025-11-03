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
	
	double epsilon_idx = dev_epsilon1[idx];
	double aux = epsilon_idx/re;
	double areau_e_j = dev_areau_e[j];
	double areau_w_j = dev_areau_w[j]; 
	double areau_n_i = dev_areau_n[i];
	double areau_s_i = dev_areau_s[i];
	double xm_ip1 = dev_xm[i+1];
	double xm_i = dev_xm[i];
	double xm_im1 = dev_xm[i-1];
	double y_jp1 = dev_y[j+1];
	double y_j = dev_y[j];
	double y_jm1 = dev_y[j-1];
	double x_i = dev_x[i];
	double x_im1 = dev_x[i-1];
	double ym_jp1 = dev_ym[j+1];
	double ym_j = dev_ym[j];
	
	if(i <= imax-1 && j <= jmax-2){
		double inv_epsilon = 1.0 / epsilon_idx;
		fn = 0.5 * (dev_vm[i*(jmax+2)+(j+1)] + dev_vm[(i-1)*(jmax+2)+(j+1)]) * areau_n_i * inv_epsilon;
		fs = 0.5 * (dev_vm[i*(jmax+2)+j] + dev_vm[(i-1)*(jmax+2)+j]) * areau_s_i * inv_epsilon;
		fe = 0.5 * (dev_um[(i+1)*(jmax+1)+j] + dev_um[idx]) * areau_e_j * inv_epsilon;
		fw = 0.5 * (dev_um[idx] + dev_um[(i-1)*(jmax+1)+j]) * areau_w_j * inv_epsilon;
		df = fe - fw + fn - fs;
		
		double dy_n = y_jp1 - y_j;
		double dy_s = y_j - y_jm1;
		double dx_e = xm_ip1 - xm_i;
		double dx_w = xm_i - xm_im1;
		
		dn = aux * areau_n_i / dy_n;
		ds = aux * areau_s_i / dy_s;
		de = aux * areau_e_j / dx_e;
		dw = aux * areau_w_j / dx_w;

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

		dudxdx = areau_e_j * (u_e - u_p) / dx_e - areau_w_j * (u_p - u_w) / dx_w;
	
		double dy_ym = ym_jp1 - ym_j;
		double vm_n = dev_vm[i*(jmax+2)+(j+1)];
		double vm_s = dev_vm[i*(jmax+2)+j];
		double vm_nw = dev_vm[(i-1)*(jmax+2)+(j+1)];
		double vm_sw = dev_vm[(i-1)*(jmax+2)+j];
		
		dxdvdy = areau_e_j * (vm_n - vm_s) / dy_ym - areau_w_j * (vm_nw - vm_sw) / dy_ym;

		artdivu = -b_art * (dudxdx + dxdvdy);

		//bulk artificial viscosity term from Ramshaw(1990)
		double dx_main = x_i - x_im1;
		double dy_main = y_j - y_jm1;
		q_art = epsilon_idx * (dev_p[idx] - dev_p[(i-1)*(jmax+1)+j]) / dx_main + artdivu;

		double inv_vol = 1.0 / (dx_main * dy_main);
		double velocidade_mag = sqrt(u_p*u_p + v_p*v_p);
		double darcy_term = u_p / (re * darcy_number);
		double forchheimer_coef = cf / sqrt(epsilon_idx * darcy_number);
		double porous_drag = epsilon_idx * (darcy_term + forchheimer_coef * u_p * velocidade_mag) * dev_liga_poros[idx];
		
		dev_ru[idx] = inv_vol * (-ap * u_p
				+  aww * u_ww + aw * u_w
				+  aee * u_ee + ae * u_e
				+  ass * u_ss + as * u_s
				+  ann * u_nn + an * u_n)
				-  q_art - porous_drag - g * epsilon_idx;
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

    upwind_U_pair(dev_um, dev_vm, dev_p, dev_ru);
}