package Quickstart;

import java.io.*;
import java.util.Scanner;
import static java.lang.Math.pow;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;


public class Code extends Thread{

    private Scanner VScan, IScan;
    private String VString, IString;
    private double I, V; //input
    private double SOC_est, V_est, Vd_est, Vh_est, Cmax_est;                                    //output
    private double Rs_est, Rc_est, Cd_est, rho_est, Vh_max_est;                                //output
    private double[][] Dwork1 = new double[3][1];
    private double[][] Dwork2 = new double[1][36];
    private double Dwork3, Dwork4, Dwork5, Dwork6, Dwork7, Dwork8;       //dwork variable
    private double Dwork9, Dwork10, Dwork11, Dwork12, Dwork14, Dwork15; // dwork variable
    private double[][] Dwork13 = new double[1][18];
    private int length;
    StringBuffer s;
    double d[][];
  
    
    

    public void mul(int rwa, int cla, int rwb, int clb, double[][] a, double[][] b, double[][] c) {
        int i, j, k;
        for (i = 0; i < rwa; i++) {

            for (j = 0; j < clb; j++) {
                double sum = 0;
                for (k = 0; k < rwb; k++) {
                    sum += a[i][k] * b[k][j];
                }
                c[i][j] = sum;
            }
        }

    }

    public void add(int rw, int cl, double[][] a, double[][] b, double[][] c) {
        int i, j;
        for (i = 0; i < rw; i++) {

            for (j = 0; j < cl; j++) {

                c[i][j] = a[i][j] + b[i][j];

            }
        }

    }

    public void assign(int rw, int cl, double[][] a, double[][] b) {
        int i, j;
        for (i = 0; i < rw; i++) {

            for (j = 0; j < cl; j++) {

                b[i][j] = a[i][j];

            }
        }

    }

    public void const_mul(int rw, int cl, double val, double[][] a, double[][] b) {
        int i, j;
        for (i = 0; i < rw; i++) {

            for (j = 0; j < cl; j++) {

                b[i][j] = a[i][j] * val;

            }
        }

    }

    double Voc(double SOC) {
        return (-0.852 * Math.exp(-63.867 * SOC) + 3.6297 + 0.559 * SOC - 0.51 * Math.pow(SOC, 2) + 0.508 * Math.pow(SOC, 3));
    }

    double sat(double e, double bound) {
        if (e >= bound) {
            return 1;
        } else if (e <= -bound) {
            return -1;
        } else {
            return (e / bound);
        }
    }

    void initialize() throws FileNotFoundException, IOException {
        
              s=new StringBuffer();
         
        int i;
        length=0;
        Dwork1[0][0] = 0.8;
        Dwork1[1][0] = 0;
        Dwork1[2][0] = 0;//Xint [SOC_int, Vd, vh] SOC_int =0.8
        for (i = 0; i < 36; i++) {
            Dwork2[0][i] = 0;
        }
        Dwork2[0][0] = Math.pow(10, -6);
        Dwork2[0][7] = Math.pow(10, -9);
        Dwork2[0][14] = Math.pow(10, -13);
        Dwork2[0][21] = Math.pow(10, -5);
        Dwork2[0][28] = Math.pow(10, -9);
        Dwork2[0][35] = Math.pow(10, -5);
        //_theta_pri 6*6
        Dwork3 = 0; //error_pre SVSF Design
        Dwork4 = 0.008; //Rc_int
        Dwork5 = 4000; //Cd_int
        Dwork6 = 0.02; //Rs_int
        Dwork7 = 7 * Math.pow(10, -3); //rho_int
        Dwork8 = 0.03; //Vh_max_int
        Dwork9 = 0; //Vd_pre
        Dwork10 = 0; //I_pre
        Dwork11 = Math.exp(-1 / (0.008 * 4000)); //a1_pre=a1=exp(-Ts/(Rc_int*Cd_int))
        Dwork12 = 0; //Ha_pre
        //Dwork13= {0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}; //dxdTh_pre 3*6
        for (i = 0; i < 18; i++) {
            Dwork13[0][i] = 0;
        }
        Dwork14 = 16; //Cn_int initial capacaity [Ah]
        Dwork15 = 0; //vh_pre=0 V
    }

    void calculate() {
        int i, j, k, Ts, eta;

        double a1, b1, bound, gamma, Cn_int, n_C, R_theta;
        double[][] Xinit = new double[3][1];
        double[][] P_theta_pri_i = new double[1][36];
        double error_pre, Rc, Cd, Rs, rho, Vh_max; //dwork variable
        double Vd_pre, I_pre, a1_pre, Ha_pre, Cmax, vh_pre;//dwork variable
        double[][] dxdTh_pre_i = new double[1][18];
        double[][] P_theta_pri = new double[6][6];
        double[][] dxdTh_pre = new double[3][6];

        assign(3, 1, Dwork1, Xinit);   //Xinit = Dwork1;
        assign(1, 36, Dwork2, P_theta_pri_i);      //P_theta_pri_i = Dwork2;
        error_pre = Dwork3;
        Rc = Dwork4;
        Cd = Dwork5;   //This may be used later. So let's keep it.
        Rs = Dwork6;
        rho = Dwork7;
        Vh_max = Dwork8;
        Vd_pre = Dwork9;
        I_pre = Dwork10;
        a1_pre = Dwork11;
        Ha_pre = Dwork12;
        assign(1, 18, Dwork13, dxdTh_pre_i);   //dxdTh_pre_i = Dwork13;
        Cmax = Dwork14;
        vh_pre = Dwork15;

        Ts = 1;
        eta = 1;
        a1 = a1_pre;
        b1 = Rc * (1 - a1);
        Cn_int = Cmax * 3600;
        bound = 0.2;
        gamma = 0.1;
        n_C = 6;
        R_theta = pow(0.4, 2);
        double Q_theta[][] = {{pow(10, -7), 0, 0, 0, 0, 0}, {0, pow(10, -10), 0, 0, 0, 0}, {0, 0, 5 * pow(10, -15), 0, 0, 0}, {0, 0, 0, pow(10, -5), 0, 0}, {0, 0, 0, 0, pow(10, -10), 0}, {0, 0, 0, 0, 0, pow(10, -6)}};
        double Theta_est[][] = {{a1}, {b1}, {1 / Cn_int}, {Rs}, {rho}, {Vh_max}};

        for (i = 0; i < 6; i++) {
            for (j = 0; j < 6; j++) {
                P_theta_pri[i][j] = P_theta_pri_i[0][(j) + (i) * 6];
            }
        } // index change i-1 in C

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 6; j++) {
                dxdTh_pre[i][j] = dxdTh_pre_i[0][(j) + (i) * 6];
            }
        } // index change i-1 in C

        double Ha = Math.exp((-1) * Theta_est[5 - 1][0] * Math.floor(I) * Ts); // index change i-1 in C
        //System.out.println(Math.floor(I));
        double Ak_SVSF[][] = {{1, 0, 0}, {0, Theta_est[1 - 1][0], 0}, {0, 0, Ha}}; // index change i-1 in C
        double Bk_SVSF[][] = {{(-1) * Ts * Theta_est[3 - 1][0], 0}, {Theta_est[2 - 1][0], 0}, {0, (Ha - 1)}}; // index change i-1 in C

        //------------Calculate Xest_pri----------------matrix mul-----------------
        double[][] Xest_pri = new double[3][1];
        double[][] temp1 = new double[3][1];
        mul(3, 3, 3, 1, Ak_SVSF, Xinit, temp1); //temp1[3][1] = Ak_SVSF*Xinit
        double[][] temp2 = new double[3][1];
        double Im[][] = {{I}, {((I > 0) ? 1 : ((I < 0) ? -1 : 0))}};
        mul(3, 2, 2, 1, Bk_SVSF, Im, temp2); // temp2[3][1] = Bk_SVSF*[I;sign(I)]
        add(3, 1, temp1, temp2, Xest_pri);//Xest_pri = temp1 + temp2

        //------------------------------
        double temp = 54.4147 * Math.exp(-63.867 * Xest_pri[1 - 1][0]) + 0.559 - 1.02 * Xest_pri[1 - 1][0] + 1.524 * Math.pow(Xest_pri[1 - 1][0], 2);
        double C[][] = {{temp, -1, Theta_est[6 - 1][0]}}; //matlab to C index i-1

        double error = V - (Voc(Xest_pri[1 - 1][0]) - I * Theta_est[4 - 1][0] - Xest_pri[2 - 1][0] + Xest_pri[3 - 1][0] * Theta_est[6 - 1][0]);
        double Y = sat(error, bound);//saturation function
        //------------------------------------------------------------------
        double[][] Cprime = new double[3][1];
        for (i = 0; i < 3; i++) {
            Cprime[i][0] = C[0][i];
        }
        double CCprime = 0;
        CCprime = C[0][0] * Cprime[0][0] + C[0][1] * Cprime[1][0] + C[0][2] * Cprime[2][0];
        double[][] K_SVSF = new double[3][1];
        for (i = 0; i < 3; i++) {
            K_SVSF[i][0] = Cprime[i][0] / (CCprime + Math.pow(10, -12)) * (Math.floor(error) + gamma * Math.floor(error_pre)) * Y / (error + Math.pow(10, -12));
        }
        //----------------------------------------------------------------------------------
        double[][] K_SVSFerror = new double[3][1];
        double[][] Xest = new double[3][1];
        const_mul(3, 1, error, K_SVSF, K_SVSFerror);  //K_SVSFerror[i][0] = K_SVSF[i][0]*error;
        add(3, 1, Xest_pri, K_SVSFerror, Xest); //Xest[i][0] = Xest_pri[i][0]+K_SVSFerror[i][0];
        V_est = (Voc(Xest[1 - 1][0]) - I * Theta_est[4 - 1][0] - Xest[2 - 1][0] + Xest[3 - 1][0] * Theta_est[6 - 1][0]);
        error_pre = V - V_est;
        //---------------------------------------
        assign(3, 1, Xest, Xinit);  //Xinit[i][0]= Xest[i][0];
        //-----------------------------------------------------
        double[][] P_theta_pos = new double[6][6];       //P_theta_pos[i][j]=P_theta_pri[i][j]+Q_theta[i][j];
        add(6, 6, P_theta_pri, Q_theta, P_theta_pos);
        //----------------------------------------------------------
        double diag1[][] = {{1, 0, 0}, {0, a1_pre, 0}, {0, 0, Ha_pre}};
        double[][] diag1dxdTh_pre = new double[3][6];
        mul(3, 3, 3, 6, diag1, dxdTh_pre, diag1dxdTh_pre);
        double dxdTh_[][] = {{0, 0, -eta * Ts * I_pre, 0, 0, 0}, {Vd_pre, I_pre, 0, 0, 0, 0}, {0, 0, 0, 0, -Math.floor(I_pre) * Ts * Ha_pre * (vh_pre + ((I_pre > 0) ? 1 : ((I_pre < 0) ? -1 : 0))), 0}};
        add(3, 6, dxdTh_, diag1dxdTh_pre, dxdTh_);

        double dGdTh[][] = {{0, 0, 0, -I, 0, Xest[3 - 1][0]}};
        double[][] CdxdTh = new double[1][6];
        mul(1, 3, 3, 6, C, dxdTh_, CdxdTh);  // CdxdTh =C * dxdTh_
        add(1, 6, dGdTh, CdxdTh, dGdTh);   //dGdTh=CdxdTh + dGdTh
        double[][] H_theta = new double[1][6];
        assign(1, 6, dGdTh, H_theta);      //H_theta = dGdTh
        //---------------------------------------------------------------
        //K_theta = P_theta_pos*H_theta'*inv(H_theta*P_theta_pos*H_theta'+R_theta);
        double[][] H_thetaPrime = new double[6][1];
        double[][] K_theta = new double[6][1];
        double[][] temp3 = new double[6][1];
        double[][] temp4 = new double[1][1];
        double temp5;

        for (i = 0; i < 6; i++) {
            H_thetaPrime[i][0] = H_theta[0][i];
        }
        mul(6, 6, 6, 1, P_theta_pos, H_thetaPrime, temp3);
        mul(1, 6, 6, 1, H_theta, temp3, temp4);
        temp5 = temp4[0][0] + R_theta;
        temp5 = 1 / temp5;
        const_mul(6, 1, temp5, temp3, K_theta);
        //-----------------------------------------------------------------------
        //--Theta_new = Theta_est + K_theta*error;
        double[][] K_thetaerror = new double[6][1];
        double[][] Theta_new = new double[6][1];
        const_mul(6, 1, error, K_theta, K_thetaerror);
        add(6, 1, Theta_est, K_thetaerror, Theta_new);  // Theta_new = Theta_est + K_thetaerror;
        //--------------------------------------------------------
        //---P_theta = (eye(n_C)-K_theta*H_theta)*P_theta_pos;
        double[][] P_theta = new double[6][6];
        double eyenc[][] = {{1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}};

        double[][] KHtheta = new double[6][6];
        mul(6, 1, 1, 6, K_theta, H_theta, KHtheta);
        const_mul(6, 6, -1, KHtheta, KHtheta);
        add(6, 6, eyenc, KHtheta, KHtheta);
        mul(6, 6, 6, 6, KHtheta, P_theta_pos, P_theta);

        //-----P_theta_pri= P_theta
        assign(6, 6, P_theta, P_theta_pri);

        //--------dxdTh_pre=dxdTh_
        assign(3, 6, dxdTh_, dxdTh_pre);

        for (i = 0; i < 6; i++) {
            for (j = 0; j < 6; j++) {
                P_theta_pri_i[0][j + (i) * 6] = P_theta_pri[i][j];
            }
        }

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 6; j++) {
                dxdTh_pre_i[0][j + (i) * 6] = dxdTh_pre[i][j];
            }
        }
        Vh_max_est = Theta_new[6 - 1][0];         //index   i-1
        rho_est = Theta_new[5 - 1][0];            //index i-1
        Rc_est = Theta_new[2 - 1][0] / (1 - Theta_new[1 - 1][0]);  //index i-1
        Cd_est = -Ts / (Rc_est * Math.log(Theta_new[1 - 1][0]));     //index i-1
        Cmax_est = 1 / (Theta_new[3 - 1][0] * 3600);            //index i-1
        Rs_est = Theta_new[4 - 1][0];                       //index i-1
        SOC_est = Xinit[1 - 1][0];                          //index i-1
        Vd_est = Xinit[2 - 1][0];                       //index i-1
        Vh_est = Xinit[3 - 1][0] * Vh_max_est;            //index i-1
        // assigning temporary values to D work variable in c assign to eeprom memory
        assign(3, 1, Xest, Dwork1);         //   Dwork1= Xinit;
        assign(1, 36, P_theta_pri_i, Dwork2);  //Dwork2= P_theta_pri_i;
        Dwork3 = error_pre;
        Dwork4 = Rc_est;
        Dwork5 = Cd_est;
        Dwork6 = Rs_est;
        Dwork7 = rho_est;
        Dwork8 = Vh_max_est;
        Dwork9 = Vd_est; //Vd_pre
        Dwork10 = I; //I_pre
        Dwork11 = Theta_new[1 - 1][0]; //1_pre  index  1-1
        Dwork12 = Ha; // Ha_pre
        assign(1, 18, dxdTh_pre_i, Dwork13);         //)Dwork13= dxdTh_pre_i;
        Dwork14 = Cmax_est;
        Dwork15 = Xest[3 - 1][0];  //vh_pre  index i-1
        
        
      
        
        
      /*  System.out.println("\n\n Output Parameters are\n\n");
	System.out.println("\n SOC_est   = " + SOC_est);
	System.out.println("\n V_est   =  " +V_est);
	System.out.println("\n Vd_est   = " + Vd_est);
	System.out.println("\n Vh_est   = " + Vh_est);
	System.out.println("\n Cmax_est   = " + Cmax_est);
	System.out.println("\n Rs_est   = " + Rs_est);
	System.out.println("\n Rc_est   = " + Rc_est);
	System.out.println("\n Cd_est   = " + Cd_est);
	System.out.println("\n rho_est   = " + rho_est);
	System.out.println("\n Vh_max_est   = " + Vh_max_est);
        */
        s.append(SOC_est+"\n");
        s.append(V_est+"\n");
        s.append(Vd_est+"\n");
        s.append(Vh_est+"\n");
        s.append(Cmax_est+"\n");
        s.append(Rs_est+"\n");
        s.append(Rc_est+"\n");
        s.append(Cd_est+"\n");
        s.append(rho_est+"\n");
        s.append(Vh_max_est+"\n");
        
        

    }

   
    void readXlsx( List<List<Object>> values1) throws FileNotFoundException, IOException, InterruptedException {
        initialize();
       

         
          List<List<Object>> values=values1;
    
             for (List row : values) {
            // Print columns A and E, which correspond to indices 0 and 4.
            System.out.printf("%s, %s\n", row.get(0), row.get(1));
            
            I=Double.parseDouble(row.get(0).toString());
            V=Double.parseDouble(row.get(1).toString());
            calculate();
            
            Thread.sleep(100);
            
          }
            

        }
     
}
