#include<iostream>
#include<bits/stdc++.h>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string>
#include<omp.h>
#include<nvtx3/nvToolsExt.h>
using namespace std;

const double pi = 2.0*atan(1.0/0.0);
const double api = atan(1.0/0.0);
float Re = 300.0;
int GSite = 50, GSite2 = 350;
int G_Str = 1;
float alpha_SOR=1.5;
double dt=0.0002; // overall time step
double dt1=dt/G_Str; // time step to solve flapping filament
int o = 5000001, writeInterval=982;
double amp_x = 0, omega_x = 0, amp_y = 0.0f, omega_y = 0.00f; //body kinematics
double pitch, pitch_a = 0.0, omega_pitch = 0.0, pitch_dot, phase = 0;
double u_c_rigid = 0.0, v_c_rigid = 0.0 ;
double time1 = 0.0*dt, lift3=0.0;
//double hff=radius_b;

int main()
{  
    double theta; 
    int Nf1=201, Nf2=16, counter_f;
    double Lf1=1.0, Lf2=0.3;
    double ds1=Lf1/(Nf1-1), ds2=Lf2/(Nf2-1);
    double Lf=Lf1+Lf2;
    int Nf=Nf1+Nf2-1;
    int Nt=Nf+1;
    double radius_a = Lf1/2.0, radius_b = 0.12*radius_a;
    //double ds=Lf/(Nf-1);

    float x_lower_left = -12.50, x_lower_right = 50.0, y_lower_left = -15.0, y_upper_left = 15.0;
    float Hx2 = 2.7, Hx1 = -0.3, Hy2 = 1.5, Hy1 = -1.5;
    int m_div1/* = 15*/, m_div2 = 600, m_div3/* = 75*/, m, n_div1/* = 5*/, n_div2 = 600, n_div3/* = 15*/, n,m_u,n_u,m_v,n_v,c5,probe_i,probe_j;
    int m_div11, m_div22, n_div11, n_div22, c11, c22, c33, c44, c11_f, c22_f, c33_f, c44_f;
    int nn1=3, nn2=10;   // Number of pannels in between two nodes for interpolation
    int N_markerf = Nf+Nf, N_markerfc = Nf+Nf+1, N_markerf1=Nf1+Nf1-1, N_markerf2=Nf2+Nf2+1, N_marker=(N_markerf1-1)*nn1 + (N_markerf2-1)*nn2; //Please take a number divisible by 4
    int ctr1, ctr2;
    float rate_xu = 1.1, rate_xd = 1.02, rate_yu = 1.02, rate_yd = 1.02;
    float q_x1,q_x2,q_y1,q_y2,d_x1,d_x2,d_y1,d_y2,dx1, dx2, dy1, dy2, dx, dy;
    float probe_location = 5.50; /*with respect to origin along x-axis */
    cout<<N_markerf<< " "<< N_marker <<endl;

    // Filament parameters
    double hff_filament, hff;
    int temp_filament;
    int X_SOR=20000,T_SOR=20000;
    double alpha_SOR_f=1.9;
    double kappa=0.0*pi, gamma=0.1, Fr=0.0, g_xg=0.0, g_yg=0.0, beta=1.0;
    double xf[Nf], yf[Nf],xf_old1[Nf], yf_old1[Nf],xf_old2[Nf], yf_old2[Nf], Fb_x[Nf], Fb_y[Nf],x_star[Nf],y_star[Nf];
    double Zeta[Nt],A_f[Nt],B_f[Nt],C_f[Nt],A_t[Nt],B_t[Nt],C_t[Nt];
    double xf_0=0.0,yf_0=0.0;
    double u_c[N_marker], v_c[N_marker];
    double Zeta_oldf[Nt], xf_oldf[Nf], yf_oldf[Nf], tolf=0.000000000001, max_errorf1[T_SOR],max_errorf2[X_SOR],max_errorf3[X_SOR], hf[Nf], thetaf[Nf], thetaf_global[Nf];
    double yfc[Nt], xfc[Nt];
    int cf1, cf2, cf3, cf4;
    double temp_Df[Nf],temp_Dfc,temp_Lf[Nf],temp_Lfc,tempf[Nf],tempfc,temp1f[Nf],temp1fc,drag3f[Nf],lift3f[Nf];
    double r_f0,r_f1,r_f2,theta_intp;
    double x_sfc[N_markerfc],y_sfc[N_markerfc],thetafc[Nt], minfx1, minfx2, maxfx1, maxfx2, minfy1, minfy2, maxfy1, maxfy2;
    double Afc1,Afc2,Afc3,Afc4,Bfc1,Bfc2,Bfc3,Bfc4,Cfc1,Cfc2,Cfc3,Cfc4;
    double amplitude_f=0.375, theta_amplitude=0.0*pi, omega_f=4.0, phi_pitching=0.0, phi_heaving=0.0, theta_pitching;
    double x_ss[Nf], y_ss[Nf],xff[Nf];
    double yf_c;
    double x_s_max, x_s_min, y_s_max, y_s_min;
    //cout<<fixed <<setprecision(20)<<ds<<endl;
    theta = (pi)/(Nf1-1);
    for (int i=0; i<Nf1; i++)
    {
		xf[i]=radius_a*cos(pi-theta*i) + radius_a;
		yf[i]=0.0;
    }
    for (int i=Nf1; i<Nf; i++)
    {
		xf[i]=Lf1 + (i-Nf1+1)*ds2*cos(kappa);
		yf[i]=yf_0 + (i-Nf1+1)*ds2*sin(kappa);
    }
    for (int i=0; i<Nf; i++)
    {
		xff[i]=xf[i];
    }
    for (int i=0; i<Nf1; i++)
    {
        //hf[i]=hff/2.0;
		hf[i]=abs(radius_b*pow((1-(pow((xf[i]-radius_a),2)/pow(radius_a,2))),0.5));
		//cout<< "Height" << " " << i << " " << hf[i] << endl;
    }
    //Find the location for the filament to be attached at trailing edge of elliptical foil
    for (int i=(Nf1-1)/2; i<Nf1; i++)
    {
        if (hf[i]<=0.01)
        {
	        hff_filament=hf[i];
	        temp_filament=i;
	        //cout << "hff_filament" << hff_filament << " " << temp_filament << endl;
	        break;
        }
	
    }
   //temp_filament=Nf1-1;
   //hff_filament=0.01;
    hff = 2.0*hff_filament; // thickness of the tail filament
    for (int i=temp_filament; i<Nf; i++)
    {
		hf[i]=hff_filament;
		//cout << hf[i] << endl;
    }
    /*for (int i=0; i<Nf; i++)
    {
		cout << i << "count" << hf[i] << endl;
    }*/
// creating the chordlines at initial conditions    
    for (int i=0; i<Nf1; i++)
    {
		yf[i]=0.0;
		xf[i]=xff[i];
    }	
    yf_c = amplitude_f*sin(omega_f*(time1)-phi_heaving);
    theta_pitching = theta_amplitude*sin(omega_f*(time1)-phi_pitching);
    for (int i=0; i<Nf1; i++)
    {
		x_ss[i]=(xf[i]-radius_a)*cos(theta_pitching)+yf[i]*sin(theta_pitching)+radius_a;
		y_ss[i]=-(xf[i]-radius_a)*sin(theta_pitching)+yf[i]*cos(theta_pitching)+yf_c;
    }
    for (int i=0; i<Nf1; i++)
    {
		xf[i]=x_ss[i];
		yf[i]=y_ss[i];
    }    
    for (int i=Nf1; i<Nf; i++)
    {
		yf[i]=yf_c;
		xf[i]=xf[i];
    }	
// end    
    for (int i=0; i<Nf; i++)
    {
		xf_old1[i]=xf[i];
		yf_old1[i]=yf[i];
		xf_old2[i]=xf[i];
		yf_old2[i]=yf[i];
    }
    for (int i=0;i<Nt;i++)
    {
        if(i==0)
        {
            xfc[i]=xf_old1[i];
            yfc[i]=yf_old1[i];
        }
        else if(i==Nt-1)
        {
            xfc[i]=xf_old1[i-1];
            yfc[i]=yf_old1[i-1];
        }
        else
        {
            xfc[i]=(xf_old1[i-1]+xf_old1[i])/2.0;
            yfc[i]=(yf_old1[i-1]+yf_old1[i])/2.0;
        }
    }    
    for (int i=0;i<Nf;i++)
    {
        Fb_x[i]=0.0;
        Fb_y[i]=0.0;
        x_star[i]=0.0;
        y_star[i]=0.0;
        xf_oldf[i]=0.0;
        yf_oldf[i]=0.0;
        temp_Df[i]=0.0;
        temp_Lf[i]=0.0;
        tempf[i]=0.0;
        temp1f[i]=0.0;
        drag3f[i]=0.0;
        lift3f[i]=0.0;
        thetaf[i]=0.0;
        thetaf_global[i]=0.0;
    }
    for (int i=0;i<Nt;i++)
    {
        Zeta[i]=0.0;
        Zeta_oldf[i]=0.0;
        A_f[i]=0.0;
        B_f[i]=0.0;
        C_f[i]=0.0;
        A_t[i]=0.0;
        B_t[i]=0.0;
        C_t[i]=0.0;
	thetafc[i]=0.0;
    }
    for (int i=0; i<N_marker; i++)
    {
        u_c[i]=0.0;
        v_c[i]=0.0;
    }
// filament parameters end
    dx2 = (Hx2-Hx1)/m_div2;
    m_div1 = ceil((log((((Hx1-x_lower_left)/dx2)*(rate_xu-1))+1))/(log(rate_xu)));
    m_div3 = ceil((log((((x_lower_right-Hx2)/dx2)*(rate_xd-1))+1))/(log(rate_xd)));
    cout << m_div1 << " " << m_div3 << endl;
    d_x1 = (dx2*((pow(rate_xu,m_div1))-1))/(rate_xu-1);

    dy2 = (Hy2-Hy1)/n_div2;
    n_div1 = ceil((log((((Hy1-y_lower_left)/dy2)*(rate_yu-1))+1))/(log(rate_yu)));
    n_div3 = ceil((log((((y_upper_left-Hy2)/dy2)*(rate_yd-1))+1))/(log(rate_yd)));
    cout << n_div1 << " " << n_div3 << endl;
    d_y1 = (dy2*((pow(rate_yu,n_div1))-1))/(rate_yu-1);

    m=m_div1+m_div2+m_div3+2;
    n=n_div1+n_div2+n_div3+2;
    m_u = m-1;
    n_u = n;
    m_v = m;
    n_v = n-1;
    cout << m << "  " << n << endl;

    c5 = ceil((log((((probe_location-Hx2)/dx2)*(rate_xd-1))+1))/(log(rate_xd)));
    probe_i = n_div1+ceil((n_div2/2))+1;
    probe_j = m_div1+m_div2+c5;

    float x[n-1][m-1], y[n-1][m-1], x1[n][m], y1[n][m], x_u[n_u][m_u], y_u[n_u][m_u], x_v[n_v][m_v], y_v[n_v][m_v];
    int flag[n][m], flag_tu[n_u][m_u], flag_u[n_u][m_u], flag_tv[n_v][m_v], flag_v[n_v][m_v], flag1[n][m];
    double x_sf[N_markerf], y_sf[N_markerf], x_s_old[N_marker], y_s_old[N_marker], x_s[N_marker], y_s[N_marker], n1[N_marker], n2[N_marker], t1[N_marker], t2[N_marker],v1,v2;
    float distance, min_dis, phi, temp;
    int c1,c2,c3,c4, k1, sig, n1_i, n1_j;
    const float dist = sqrt(dx2*dx2+dy2*dy2);
    //Calculation of the coordinates of the grid points
    for(int i=0;i<n-1;i++)
    {
        dx1 = dx2*(pow(rate_xu,(m_div1-1)));
        for(int j=0;j<m-1;j++)
        {
            if (j<=m_div1)
            {
                if (j==0)
                    q_x1=Hx1-d_x1;
                else
                {
                    q_x1 = q_x1 + dx1;
                    dx1 = dx1/rate_xu;
                }
                x[i][j] = q_x1;
            }
            else if(j<=m_div1+m_div2)
            {
                x[i][j] = q_x1 + (j-m_div1)*dx2;
                q_x2 = x[i][j];
                d_x2 = dx2;
            }
            else if(j<=m_div1+m_div2+m_div3)
            {
                q_x2 = q_x2 + d_x2;
                x[i][j] = q_x2;
                d_x2 = d_x2 * rate_xd;
            }
        }
    }
    for(int j=0;j<m-1;j++)
    {
        dy1 = dy2*(pow(rate_yu,(n_div1-1)));
        for(int i=0;i<n-1;i++)
        {
            if (i<=n_div1)
            {
                if (i==0)
                    q_y1=Hy1-d_y1;
                else
                {
                    q_y1 = q_y1 + dy1;
                    dy1 = dy1/rate_yu;
                }
                y[i][j] = q_y1;
            }
            else if(i<=n_div1+n_div2)
            {
                y[i][j] = q_y1 + (i-n_div1)*dy2;
                q_y2 = y[i][j];
                d_y2 = dy2;
            }
            else if(i<=n_div1+n_div2+n_div3)
            {
                q_y2 = q_y2 + d_y2;
                y[i][j] = q_y2;
                d_y2 = d_y2 * rate_yd;
            }
        }
    }
    //Calculation of the coordinate of the center of each grid cell
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            if (j==0)
            {
                if(i<n-1)
                    x1[i][j] = x[i][j];
                else
                    x1[i][j] = x[i-1][j];
            }
            else if (j==m-1)
            {
                if(i<n-1)
                    x1[i][j] = x[i][j-1];
                else
                    x1[i][j] = x[i-1][j-1];
            }
            else
            {
                if(i<n-1)
                    x1[i][j] = (x[i][j]+x[i][j-1])/2.0;
                else
                    x1[i][j] = (x[i-1][j]+x[i-1][j-1])/2.0;
            }

            if (i==0)
            {
                if(j<m-1)
                    y1[i][j] = y[i][j];
                else
                    y1[i][j] = y[i][j-1];
            }
            else if (i==n-1)
            {
                if(j<m-1)
                    y1[i][j] = y[i-1][j];
                else
                    y1[i][j] = y[i-1][j-1];
            }
            else
            {
                if(j<m-1)
                    y1[i][j] = (y[i][j]+y[i-1][j])/2.0;
                else
                    y1[i][j] = (y[i][j-1]+y[i-1][j-1])/2.0;
            }

        }
    }
    //Calculation of the coordinate of the u and v velocity locations
    for (int i=0;i<n_u;i++)
    {
        for(int j=0;j<m_u;j++)
        {
            x_u[i][j] = x[0][j];
            y_u[i][j] = y1[i][0];
        }
    }
    for (int i=0;i<n_v;i++)
    {
        for(int j=0;j<m_v;j++)
        {
            x_v[i][j] = x1[0][j];
            y_v[i][j] = y[i][0];
        }
    }
    cout << probe_i << " " << probe_j << endl;
    cout << x[probe_i][probe_j] << " " << y[probe_i][probe_j] << endl;
    cout << x1[probe_i][probe_j] << " " << y1[probe_i][probe_j] << endl;
    cout << x_u[probe_i][probe_j] << " " << y_u[probe_i][probe_j] << endl;
    cout << x_v[probe_i][probe_j] << " " << y_v[probe_i][probe_j] << endl;
    //Defining the Solid body
    //for upper surface
    for (int i=0; i<Nf; i++)
    {
        if (i==0)
        {
			thetaf[i]=atan((yf_old1[i+1]-yf_old1[i])/(xf_old1[i+1]-xf_old1[i]));
        }
        else if (i==Nf1-1)
        {
			thetaf[i]=atan((yf_old1[i]-yf_old1[i-1])/(xf_old1[i]-xf_old1[i-1]));
        }
        else if (i==Nf-1)
        {
			thetaf[i]=atan((yf_old1[i]-yf_old1[i-1])/(xf_old1[i]-xf_old1[i-1]));
        }
        else
        {
			thetaf[i]=atan((yf_old1[i+1]-yf_old1[i-1])/(xf_old1[i+1]-xf_old1[i-1]));
        }
    }
    // to find global slope
    /*for (int i=0; i<Nf; i++)
    {
        if (i==0)
        {
           thetaf_global[i]=thetaf[i];
        }
        else
        {
           thetaf_global[i]=atan((yf_old1[i]-yf_old1[0])/(xf_old1[i]-xf_old1[0]));
        }
        thetaf_global[i]=atan((yf_old1[i]-yf_old1[0])/(xf_old1[i]-xf_old1[0]));
    }*/
    // to find the surface points
    // for upper surface
    for (int i=0; i<Nf; i++)
    {
        if (i==0)
        {
		x_sf[i]=xf_old1[i];
		y_sf[i]=yf_old1[i];
        }
		else if (i==Nf-1)
		{
			if ((xf_old1[i]-xf_old1[i-1])>=0.0 && thetaf[i]>=0.0)
			{
			x_sf[i]=xf_old1[i] - hf[i]*sin(thetaf[i]);
			y_sf[i]=yf_old1[i] + hf[i]*cos(thetaf[i]);
			}
			else if ((xf_old1[i]-xf_old1[i-1])<0.0 && thetaf[i]<=0.0)
			{
				x_sf[i]=xf_old1[i] + hf[i]*sin(thetaf[i]);
				y_sf[i]=yf_old1[i] - hf[i]*cos(thetaf[i]);
			}
			else if ((xf_old1[i]-xf_old1[i-1])>=0.0 && thetaf[i]<=0.0)
			{
				x_sf[i]=xf_old1[i] - hf[i]*sin(thetaf[i]);
				y_sf[i]=yf_old1[i] + hf[i]*cos(thetaf[i]);
			}
			else if ((xf_old1[i]-xf_old1[i-1])<0.0 && thetaf[i]>=0.0)
			{
				x_sf[i]=xf_old1[i] + hf[i]*sin(thetaf[i]);
				y_sf[i]=yf_old1[i] - hf[i]*cos(thetaf[i]);
			}
			cout<< " Dmension lal" << x_sf[i] << " " << y_sf[i] <<endl;
		}
		else
		{
			if ((xf_old1[i+1]-xf_old1[i-1])>=0.0 && thetaf[i]>=0.0)
			{
				x_sf[i]=xf_old1[i] - hf[i]*sin(thetaf[i]);
				y_sf[i]=yf_old1[i] + hf[i]*cos(thetaf[i]);
			}
			else if ((xf_old1[i+1]-xf_old1[i-1])<0.0 && thetaf[i]<=0.0)
			{
				x_sf[i]=xf_old1[i] + hf[i]*sin(thetaf[i]);
				y_sf[i]=yf_old1[i] - hf[i]*cos(thetaf[i]);
			}
			else if ((xf_old1[i+1]-xf_old1[i-1])>=0.0 && thetaf[i]<=0.0)
			{
				x_sf[i]=xf_old1[i] - hf[i]*sin(thetaf[i]);
				y_sf[i]=yf_old1[i] + hf[i]*cos(thetaf[i]);
			}
			else if ((xf_old1[i+1]-xf_old1[i-1])<0.0 && thetaf[i]>=0.0)
			{
				x_sf[i]=xf_old1[i] + hf[i]*sin(thetaf[1]);
				y_sf[i]=yf_old1[i] - hf[i]*cos(thetaf[i]);
			}
			else
			{
				cout<< " Top surface middle11 " <<endl;
			}
		}
    }
    // for lower surface

    for (int i=Nf; i<N_markerf; i++)
    {
        if(i==Nf)
        {
			x_sf[i]=xf_old1[N_markerf-i-1];
			y_sf[i]=yf_old1[N_markerf-i-1];
			cout<< " Dmension lal" << x_sf[i] << " " << y_sf[i] <<endl;
        }
        else if (i==Nf+1)
		{
			if ((xf_old1[N_markerf-i]-xf_old1[N_markerf-i-1])>=0.0 && thetaf[N_markerf-i]>=0.0)
			{
				x_sf[i]=xf_old1[N_markerf-i] + hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
				y_sf[i]=yf_old1[N_markerf-i] - hf[N_markerf-i]*cos(thetaf[N_markerf-i]);
			}
			else if ((xf[N_markerf-i]-xf[N_markerf-i-1])<0.0 && thetaf[N_markerf-i]<=0.0)
			{
				x_sf[i]=xf_old1[N_markerf-i] - hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
				y_sf[i]=yf_old1[N_markerf-i] + hf[N_markerf-i]*cos(thetaf[N_markerf-i]);
			}
			else if ((xf_old1[N_markerf-i]-xf_old1[N_markerf-i-1])>=0.0 && thetaf[N_markerf-i]<=0.0)
			{
				x_sf[i]=xf_old1[N_markerf-i] + hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
				y_sf[i]=yf_old1[N_markerf-i] - hf[N_markerf-i]*cos(thetaf[N_markerf-i]);
			}
			else if ((xf_old1[N_markerf-i]-xf_old1[N_markerf-i-1])<0.0 && thetaf[N_markerf-i]>=0.0)
			{
				x_sf[i]=xf_old1[N_markerf-i] - hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
				y_sf[i]=yf_old1[N_markerf-i] + hf[N_markerf-i]*cos(thetaf[N_markerf-i]);
			}
			cout<< " Dmension lal" << x_sf[i] << " " << y_sf[i] <<endl;
		}
		else
		{
			if ((xf_old1[N_markerf-i+1]-xf_old1[N_markerf-i-1])>=0.0 && thetaf[N_markerf-i]>=0.0)
			{
				x_sf[i]=xf_old1[N_markerf-i] + hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
				y_sf[i]=yf_old1[N_markerf-i] - hf[N_markerf-i]*cos(thetaf[N_markerf-i]);
			}
			else if ((xf_old1[N_markerf-i+1]-xf_old1[N_markerf-i-1])<0.0 && thetaf[N_markerf-i]<=0.0)
			{
				x_sf[i]=xf_old1[N_markerf-i] - hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
				y_sf[i]=yf_old1[N_markerf-i] + hf[N_markerf-i]*cos(thetaf[N_markerf-i]);
			}
			else if ((xf_old1[N_markerf-i+1]-xf_old1[N_markerf-i-1])>=0.0 && thetaf[N_markerf-i]<=0.0)
			{
				x_sf[i]=xf_old1[N_markerf-i] + hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
				y_sf[i]=yf_old1[N_markerf-i] - hf[N_markerf-i]*cos(thetaf[N_markerf-i]);
			}
			else if ((xf_old1[N_markerf-i+1]-xf_old1[N_markerf-i-1])<0.0 && thetaf[N_markerf-i]>=0.0)
			{
				x_sf[i]=xf_old1[N_markerf-i] - hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
				y_sf[i]=yf_old1[N_markerf-i] + hf[N_markerf-i]*cos(thetaf[N_markerf-i]);
			}
			else
			{
				cout<< " Bottom surface middle11 " <<endl;
			}
		}
    }

    x_s_max=x_sf[0];
    y_s_max=y_sf[0];
    x_s_min=x_sf[0];
    y_s_min=y_sf[0];
    for (int i=0; i<N_markerf; i++)
    {
        if(x_s_max<x_sf[i])
        {
                x_s_max=x_sf[i];
        }
        if(y_s_max<y_sf[i])
        {
                y_s_max=y_sf[i];
        }
    }
    for (int i=0; i<N_markerf; i++)
    {
        if(x_s_min>x_sf[i])
        {
                x_s_min=x_sf[i];
        }
        if(y_s_min>y_sf[i])
        {
                y_s_min=y_sf[i];
        }
    }
    
     // Boundary for flagging purpose
    m_div11 = m_div1+ceil((x_s_min-Hx1)/dx2);
    m_div22 = ceil((x_s_max-x_s_min)/dx2);
    n_div11 = n_div1+ceil((y_s_min-Hy1)/dy2);
    n_div22 = ceil((y_s_max-y_s_min)/dy2);
    
    c11 = m_div11-20;
    c22 = m_div11+m_div22+20;
    c33 = n_div11-20;
    c44 = n_div11+n_div22+20;
    // co-ordinate to compute force only for the rigid foil
    x_s_max=x_sf[0];
    y_s_max=y_sf[0];
    x_s_min=x_sf[0];
    y_s_min=y_sf[0];
    for (int i=0; i<Nf1; i++)
    {
        if(x_s_max<x_sf[i])
        {
                x_s_max=x_sf[i];
        }
        if(y_s_max<y_sf[i])
        {
                y_s_max=y_sf[i];
        }
    }
    for (int i=Nf1+N_markerf2-2; i<N_markerf; i++)
    {
        if(x_s_max<x_sf[i])
        {
                x_s_max=x_sf[i];
        }
        if(y_s_max<y_sf[i])
        {
                y_s_max=y_sf[i];
        }
    }
    for (int i=0; i<Nf1; i++)
    {
        if(x_s_min>x_sf[i])
        {
                x_s_min=x_sf[i];
        }
        if(y_s_min>y_sf[i])
        {
                y_s_min=y_sf[i];
        }
    }
    for (int i=Nf1+N_markerf2-2; i<N_markerf; i++)
    {
        if(x_s_min>x_sf[i])
        {
                x_s_min=x_sf[i];
        }
        if(y_s_min>y_sf[i])
        {
                y_s_min=y_sf[i];
        }
    }
    
     // Boundary for flagging purpose
    m_div11 = m_div1+ceil((x_s_min-Hx1)/dx2);
    m_div22 = ceil((x_s_max-x_s_min)/dx2);
    n_div11 = n_div1+ceil((y_s_min-Hy1)/dy2);
    n_div22 = ceil((y_s_max-y_s_min)/dy2);
    
    c11_f = m_div11-20;
    c22_f = m_div11+m_div22-1;
    c33_f = n_div11-4;
    c44_f = n_div11+n_div22+4;
//Generating surface points for flagging using linear interpolation
    for (int i=0; i<(Nf1*nn1)-(nn1-1); i++)
    {
		ctr1=trunc(i/nn1);
		if(i%nn1==0)
		{
			x_s[i]=x_sf[ctr1];
			y_s[i]=y_sf[ctr1];
		}
		else
		{
			dx = (x_sf[ctr1+1]-x_sf[ctr1])/(nn1);
			dy = (y_sf[ctr1+1]-y_sf[ctr1])/(nn1);
			if (atan((y_sf[ctr1+1]-y_sf[ctr1])/(x_sf[ctr1+1]-x_sf[ctr1]))==api || atan((y_sf[ctr1+1]-y_sf[ctr1])/(x_sf[ctr1+1]-x_sf[ctr1]))==-api || dx==0.0 || dx==-0.0)
			{
				x_s[i] = x_sf[ctr1];
				y_s[i] = y_sf[ctr1] + (i-ctr1*nn1)*dy;
			}
			else
			{
				x_s[i] = x_sf[ctr1] + (i-ctr1*nn1)*dx;
				y_s[i] = y_sf[ctr1] + (y_sf[ctr1+1]-y_sf[ctr1])/(x_sf[ctr1+1]-x_sf[ctr1])*(x_s[i]-x_sf[ctr1]);
			}
		}
    }
    counter_f=1;
    for (int i=(Nf1*nn1)-(nn1-1); i<((Nf2*2)*nn2+(Nf1*nn1)-(nn1-1)); i++)
    {
		ctr1=trunc((i+1-((Nf1*nn1)-(nn1-1)))/nn2);
		if((i+1-((Nf1*nn1)-(nn1-1)))%nn2==0)
		{
			x_s[i]=x_sf[(Nf1-1) + ctr1];
			y_s[i]=y_sf[(Nf1-1) + ctr1];
		}
		else
		{
			
			dx = (x_sf[(Nf1-1) + ctr1+1]-x_sf[(Nf1-1) + ctr1])/(nn2);
			dy = (y_sf[(Nf1-1) + ctr1+1]-y_sf[(Nf1-1) + ctr1])/(nn2);
			if (atan((y_sf[(Nf1-1) + ctr1+1]-y_sf[(Nf1-1) + ctr1])/(x_sf[(Nf1-1) + ctr1+1]-x_sf[(Nf1-1) + ctr1]))==api || atan((y_sf[(Nf1-1) + ctr1+1]-y_sf[(Nf1-1) + ctr1])/(x_sf[(Nf1-1) + ctr1+1]-x_sf[(Nf1-1) + ctr1]))==-api || dx==0.0 || dx==-0.0)
			{
				x_s[i] = x_sf[(Nf1-1) + ctr1];
				y_s[i] = y_sf[(Nf1-1) + ctr1] + (counter_f-ctr1*nn2)*dy;
			}
			else
			{
				x_s[i] = x_sf[(Nf1-1) + ctr1] + (counter_f-ctr1*nn2)*dx;
				y_s[i] = y_sf[(Nf1-1) + ctr1] + (y_sf[(Nf1-1) + ctr1+1]-y_sf[(Nf1-1) + ctr1])/(x_sf[(Nf1-1) + ctr1+1]-x_sf[(Nf1-1) + ctr1])*(x_s[i]-x_sf[(Nf1-1) + ctr1]);
			}
		}
		counter_f=counter_f+1;
    }
    counter_f=1;
    for (int i=((Nf2*2)*nn2+(Nf1*nn1)-(nn1-1)); i<N_marker; i++)
    {
		ctr1=trunc((i+1-((Nf2*2)*nn2+(Nf1*nn1)-(nn1-1)))/nn1);
		if((i+1-((Nf2*2)*nn2+(Nf1*nn1)-(nn1-1)))%nn1==0)
		{
			x_s[i]=x_sf[(Nf1+Nf2*2-1) + ctr1];
			y_s[i]=y_sf[(Nf1+Nf2*2-1) + ctr1];
		}
		else
		{
			//if (ctr1==N_markerf-1)
			if (((Nf1+Nf2*2-1) + ctr1)==N_markerf-1)
			{
				dx = (x_sf[0]-x_sf[(Nf1+Nf2*2-1) + ctr1])/(nn1);
				dy = (y_sf[0]-y_sf[(Nf1+Nf2*2-1) + ctr1])/(nn1);
			}
			else
			{
				dx = (x_sf[(Nf1+Nf2*2-1) + ctr1+1]-x_sf[(Nf1+Nf2*2-1) + ctr1])/(nn1);
				dy = (y_sf[(Nf1+Nf2*2-1) + ctr1+1]-y_sf[(Nf1+Nf2*2-1) + ctr1])/(nn1);
			}
			if (atan((y_sf[(Nf1+Nf2*2-1) + ctr1+1]-y_sf[(Nf1+Nf2*2-1) + ctr1])/(x_sf[(Nf1+Nf2*2-1) + ctr1+1]-x_sf[(Nf1+Nf2*2-1) + ctr1]))==api || atan((y_sf[(Nf1+Nf2*2-1) + ctr1+1]-y_sf[(Nf1+Nf2*2-1) + ctr1])/(x_sf[(Nf1+Nf2*2-1) + ctr1+1]-x_sf[(Nf1+Nf2*2-1) + ctr1]))==-api || dx==0.0 || dx==-0.0)
			{
				x_s[i] = x_sf[(Nf1+Nf2*2-1) + ctr1];
				y_s[i] = y_sf[(Nf1+Nf2*2-1) + ctr1] + (counter_f-ctr1*nn1)*dy;
			}
			else
			{
				x_s[i] = x_sf[(Nf1+Nf2*2-1) + ctr1] + (counter_f-ctr1*nn1)*dx;
				if (((Nf1+Nf2*2-1) + ctr1)==N_markerf-1)
				{
					y_s[i] = y_sf[(Nf1+Nf2*2-1) + ctr1] + (y_sf[0]-y_sf[(Nf1+Nf2*2-1) + ctr1])/(x_sf[0]-x_sf[(Nf1+Nf2*2-1) + ctr1])*(x_s[i]-x_sf[(Nf1+Nf2*2-1) + ctr1]);
				}
				else
				{
					y_s[i] = y_sf[(Nf1+Nf2*2-1) + ctr1] + (y_sf[(Nf1+Nf2*2-1) + ctr1+1]-y_sf[(Nf1+Nf2*2-1) + ctr1])/(x_sf[(Nf1+Nf2*2-1) + ctr1+1]-x_sf[(Nf1+Nf2*2-1) + ctr1])*(x_s[i]-x_sf[(Nf1+Nf2*2-1) + ctr1]);
				}
			}
		}
		counter_f=counter_f+1;
    }
    for (int i=0; i<N_marker; i++)
    {
        x_s_old[i]=x_s[i];
        y_s_old[i]=y_s[i];
    }
    for (int i = 0;i<N_marker;i++)
    {
        if(i==0)
        {
            v1 = x_s[i+1] - x_s[N_marker-1];
            v2 = y_s[i+1] - y_s[N_marker-1];
        }
        else if (i==N_marker-1)
        {
            v1 = x_s[0] - x_s[i-1];
            v2 = y_s[0] - y_s[i-1];
        }
        else
        {
            v1 = x_s[i+1] - x_s[i-1];
            v2 = y_s[i+1] - y_s[i-1];
        }
        n1[i] = -v2;
        n2[i] = v1;
        n1[i] = (-v2)/sqrt((-v2)*(-v2)+v1*v1);
        n2[i] = v1/sqrt((-v2)*(-v2)+v1*v1);
        //t1[i] = v1/sqrt(v1*v1+v2*v2);
        //t2[i] = v2/sqrt(v1*v1+v2*v2);
    }
    ofstream file("solidBoundary.dat");
    if (file.is_open())
    {
        for(int i=0; i<N_marker; i++)
        {
            file << x_s[i] << " " << y_s[i] << " " << n1[i] << " " << n2[i] << endl;
        }
        file << x_s[0] << " " << y_s[0] << " " << n1[0] << " " << n2[0];
        file.close();
    }
    //Identifying the fluid points, solid points and the forcing points
    /*c1 = m_div1 + (((xc-Hx1)-(radius_a+10*dx2))/dx2);
     c2 = m_div1+m_div2-(((Hx2-xc)-(radius_a+10*dx2))/dx2);
     c3 = n_div1 + (((yc-Hy1)-(radius_b+10*dy2))/dy2);
     c4 = n_div1+n_div2-(((Hy2-yc)-(radius_b+10*dy2))/dy2);*/
    c1 = m_div1;
    c2 = m_div1+m_div2;
    c3 = n_div1;
    c4 = n_div1+n_div2;
    cout<<"c1 = " <<c1<<" c11 = " <<c11<<endl;
    cout<<"c2 = " <<c2<<" c22 = " <<c22<<endl;
    cout<<"c3 = " <<c3<<" c33 = " <<c33<<endl;
    cout<<"c4 = " <<c4<<" c44 = " <<c44<<endl;
    for(int i=0;i<n;i++)
    {
        for (int j=0;j<m;j++)
        {
            if(j>c11 && j<c22 && i>c33 && i<c44)
            {
                min_dis = 10.0;
                for(int k=0;k<N_marker;k++)
                {
                    distance = sqrt((x_s[k]-x1[i][j])*(x_s[k]-x1[i][j]) + (y_s[k]-y1[i][j])*(y_s[k]-y1[i][j]));
                    if (distance < min_dis)
                    {
                        min_dis = distance;
                        k1 = k;
                    }
                }
                temp = (x1[i][j]-x_s[k1])*n1[k1] + (y1[i][j]-y_s[k1])*n2[k1];
                if (temp>0)
                    sig=1;
                else
                    sig=-1;
                phi = sig*min_dis;
                if (phi>0)
                    flag[i][j]=1;
                else
                    flag[i][j]=0;
            }
            else
                flag[i][j]=1;
        }
    }
    for(int i=0;i<n-1;i++)
    {
        for (int j=0;j<m-1;j++)
        {
            if(j>c11 && j<c22 && i>c33 && i<c44)
            {
                min_dis = 10.0;
                for(int k=0;k<N_marker;k++)
                {
                    distance = sqrt((x_s[k]-x[i][j])*(x_s[k]-x[i][j]) + (y_s[k]-y[i][j])*(y_s[k]-y[i][j]));
                    if (distance < min_dis)
                    {
                        min_dis = distance;
                        k1 = k;
                    }
                }
                temp = (x[i][j]-x_s[k1])*n1[k1] + (y[i][j]-y_s[k1])*n2[k1];
                if (temp>0)
                    sig=1;
                else
                    sig=-1;
                phi = sig*min_dis;
                if (phi>0)
                    flag1[i][j]=1;
                else
                    flag1[i][j]=0;
            }
            else
                flag1[i][j]=1;
        }
    }
    for(int i=0;i<n_u;i++)
    {
        for (int j=0;j<m_u;j++)
        {
            if(j>c11 && j<c22 && i>c33 && i<c44)
            {
                min_dis = 10.0;
                for(int k=0;k<N_marker;k++)
                {
                    distance = sqrt((x_s[k]-x_u[i][j])*(x_s[k]-x_u[i][j]) + (y_s[k]-y_u[i][j])*(y_s[k]-y_u[i][j]));
                    if (distance < min_dis)
                    {
                        min_dis = distance;
                        k1 = k;
                    }
                }
                temp = (x_u[i][j]-x_s[k1])*n1[k1] + (y_u[i][j]-y_s[k1])*n2[k1];
                if (temp>0)
                    sig=1;
                else
                    sig=-1;
                phi = sig*min_dis;
                if (phi>0)
                    flag_tu[i][j]=1;
                else
                    flag_tu[i][j]=0;
            }
            else
                flag_tu[i][j]=1;
        }
    }
    for(int i=0;i<n_u;i++)
    {
        for (int j=0;j<m_u;j++)
        {
            if (flag_tu[i][j]==0)
            {
                /*if(flag_tu[i+1][j]==0 && flag_tu[i+1][j+1]==0 && flag_tu[i][j+1]==0 && flag_tu[i-1][j+1]==0 && flag_tu[i-1][j]==0 && flag_tu[i-1][j-1]==0 && flag_tu[i][j-1]==0 && flag_tu[i+1][j-1]==0)
                 {
			flag_u[i][j]=0;
                 }
                 else*/
                flag_u[i][j]=3;
            }
            else
                flag_u[i][j]=1;
        }
    }
    for(int i=0;i<n_v;i++)
    {
        for (int j=0;j<m_v;j++)
        {
            if(j>c11 && j<c22 && i>c33 && i<c44)
            {
                min_dis = 10.0;
                for(int k=0;k<N_marker;k++)
                {
                    distance = sqrt((x_s[k]-x_v[i][j])*(x_s[k]-x_v[i][j]) + (y_s[k]-y_v[i][j])*(y_s[k]-y_v[i][j]));
                    if (distance < min_dis)
                    {
                        min_dis = distance;
                        k1 = k;
                    }
                }
                temp = (x_v[i][j]-x_s[k1])*n1[k1] + (y_v[i][j]-y_s[k1])*n2[k1];
                if (temp>0)
                    sig=1;
                else
                    sig=-1;
                phi = sig*min_dis;
                if (phi>0)
                    flag_tv[i][j]=1;
                else
                    flag_tv[i][j]=0;
            }
            else
                flag_tv[i][j]=1;
        }
    }
    for(int i=0;i<n_v;i++)
    {
        for (int j=0;j<m_v;j++)
        {
            if (flag_tv[i][j]==0)
            {
                /*if(flag_tv[i+1][j]==0 && flag_tv[i+1][j+1]==0 && flag_tv[i][j+1]==0 && flag_tv[i-1][j+1]==0 && flag_tv[i-1][j]==0 && flag_tv[i-1][j-1]==0 && flag_tv[i][j-1]==0 && flag_tv[i+1][j-1]==0)
                 {
                 flag_v[i][j]=0;
                 }
                 else*/
                flag_v[i][j]=3;
            }
            else
                flag_v[i][j]=1;
        }
    }
    ofstream file1("grid.dat");
    if (file1.is_open())
    {
        file1 << "ZONE T=DATA I=" << m-1 << " " << "J=" << n-1 << endl;
        for(int i=0;i<n-1;i++)
        {
            for(int j=0;j<m-1;j++)
            {
                file1 << x[i][j] << " " << y[i][j] << " " << flag1[i][j] << endl;
            }
        }
        file1.close();
    }
    ofstream file2("grid1.dat");
    if (file2.is_open())
    {
        file2 << "ZONE T=DATA I=" << m << " " << "J=" << n << endl;
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<m;j++)
            {
                file2 << x1[i][j] << " " << y1[i][j] << " " << flag[i][j] << endl;
            }
        }
        file2.close();
    }
    ofstream file3("grid_u.dat");
    if (file3.is_open())
    {
        file3 << "ZONE T=DATA I=" << m_u << " " << "J=" << n_u << endl;
        for(int i=0;i<n_u;i++)
        {
            for(int j=0;j<m_u;j++)
            {
                file3 << x_u[i][j] << " " << y_u[i][j] << " " << flag_u[i][j] << endl;
            }
        }
        file3.close();
    }
    ofstream file4("grid_v.dat");
    if (file4.is_open())
    {
        file4 << "ZONE T=DATA I=" << m_v << " " << "J=" << n_v << endl;
        for(int i=0;i<n_v;i++)
        {
            for(int j=0;j<m_v;j++)
            {
                file4 << x_v[i][j] << " " << y_v[i][j] << " " << flag_v[i][j] << endl;
            }
        }
        file4.close();
    }

    float u[n_u][m_u], v[n_v][m_v], u_old[n_u][m_u], v_old[n_v][m_v], pressure[n][m], uf1[n][m], vf1[n][m];
    float u_1[n_u][m_u], v_1[n_v][m_v];
    float ue,uw,un,us,ve,vw,vn,vs,ue_old,uw_old,un_old,us_old,ve_old,vw_old,vn_old,vs_old;
    float a_u[n_u][m_u], b_u[n_u][m_u], c_u[n_u][m_u], d_u[n_u][m_u], e_u[n_u][m_u], f_u[n_u][m_u], Bf_u[n_u][m_u];
    float a_v[n_v][m_v], b_v[n_v][m_v], c_v[n_v][m_v], d_v[n_v][m_v], e_v[n_v][m_v], f_v[n_v][m_v], Bf_v[n_v][m_v];
    float p_c[n][m], a_p[n][m], b_p[n][m], c_p[n][m], d_p[n][m], e_p[n][m], f_p[n][m];
    float big,residual_flux,error,A1,B1,C1, c_sasv = 1.0f,sum,sum1,sum2;
    double x_intersect, y_intersect, x_mirror, y_mirror, x_mirror1, y_mirror1, slope, dist_a, dist_b;
    float U_intp, V_intp, u_mirror, v_mirror, u_body, v_body;
    float uA_old[n_u][m_u], vA_old[n_v][m_v],pA_old[n][m],Bf_uA_old[n_u][m_u], Bf_vA_old[n_v][m_v],max_error1[GSite],max_error2[GSite],max_error_Bf_u,max_error_Bf_v,max_error3[GSite2], tol=0.000001;
    float surface_pr[N_marker],lift1,drag1,drag3;
    float uf, uf_n, uf_s, uf_old, vf, vf_e, vf_w, vf_old, temp1, temp_D, temp_L;
    float u_red[n_u][m_u],u_black[n_u][m_u],v_red[n_v][m_v],v_black[n_v][m_v],p_c_red[n][m],p_c_black[n][m];
    float u_red_old[n_u][m_u],u_black_old[n_u][m_u],v_red_old[n_v][m_v],v_black_old[n_v][m_v],p_c_red_old[n][m],p_c_black_old[n][m];
    //    float x_1,y_1,x_2,y_2,x_3,y_3,u_2,u_3,v_2,v_3,p_2,p_3,b2_u_x,b3_u_y,b2_v_x,b3_v_y,b1,b2,b3,Det,dpdn;

    sum2 = 0.0;
    for(int i=1;i<n_u-1;i++)
    {
        sum2 = sum2 + (y[i][m_u-1]-y[i-1][m_u-1]);
    }
     /*ifstream fileread1("Cavity_u.txt");
     if (fileread1.is_open())
     {
     for(int i=0; i<n_u; i++)
     {
     for(int j=0; j<m_u; j++)
     {
     fileread1 >> u[i][j];
     }
     }
     fileread1.close();
     }
     ifstream fileread1a("Cavity_u_old.txt");
     if (fileread1a.is_open())
     {
     for(int i=0; i<n_u; i++)
     {
     for(int j=0; j<m_u; j++)
     {
     fileread1a >> u_old[i][j];
     }
     }
     fileread1a.close();
     }
     ifstream fileread2("Cavity_v.txt");
     if (fileread2.is_open())
     {
     for(int i=0; i<n_v; i++)
     {
     for(int j=0; j<m_v; j++)
     {
     fileread2 >> v[i][j];
     }
     }
     fileread2.close();
     }
     ifstream fileread2a("Cavity_v_old.txt");
     if (fileread2a.is_open())
     {
     for(int i=0; i<n_v; i++)
     {
     for(int j=0; j<m_v; j++)
     {
     fileread2a >> v_old[i][j];
     }
     }
     fileread2a.close();
     }
     ifstream fileread3("pressure.txt");
     if (fileread3.is_open())
     {
     for(int i=0; i<n; i++)
     {
     for(int j=0; j<m; j++)
     {
     fileread3 >> pressure[i][j];
     }
     }
     fileread3.close();
     }
     ifstream fileread3a("p_c.txt");
     if (fileread3a.is_open())
     {
     for(int i=0; i<n; i++)
     {
     for(int j=0; j<m; j++)
     {
     fileread3a >> p_c[i][j];
     }
     }
     fileread3a.close();
     }
     ifstream fileread1f("Cavity_u_red.txt");
     if (fileread1f.is_open())
     {
     for(int i=0; i<n_u; i++)
     {
     for(int j=0; j<m_u; j++)
     {
     fileread1f >> u_red[i][j];
     }
     }
     fileread1f.close();
     }
     ifstream fileread1fa("Cavity_u_black.txt");
     if (fileread1fa.is_open())
     {
     for(int i=0; i<n_u; i++)
     {
     for(int j=0; j<m_u; j++)
     {
     fileread1fa >> u_black[i][j];
     }
     }
     fileread1fa.close();
     }
     ifstream fileread2f("Cavity_v_red.txt");
     if (fileread2f.is_open())
     {
     for(int i=0; i<n_v; i++)
     {
     for(int j=0; j<m_v; j++)
     {
     fileread2f >> v_red[i][j];
     }
     }
     fileread2f.close();
     }
     ifstream fileread2fa("Cavity_v_black.txt");
     if (fileread2fa.is_open())
     {
     for(int i=0; i<n_v; i++)
     {
     for(int j=0; j<m_v; j++)
     {
     fileread2fa >> v_black[i][j];
     }
     }
     fileread2fa.close();
     }
     ifstream fileread3f("p_c_red.txt");
     if (fileread3f.is_open())
     {
     for(int i=0; i<n; i++)
     {
     for(int j=0; j<m; j++)
     {
     fileread3f >> p_c_red[i][j];
     }
     }
     fileread3f.close();
     }
     ifstream fileread3fa("p_c_black.txt");
     if (fileread3fa.is_open())
     {
     for(int i=0; i<n; i++)
     {
     for(int j=0; j<m; j++)
     {
     fileread3fa >> p_c_black[i][j];
     }
     }
     fileread3fa.close();
     }*/
    //declaring the initial condition
    for(int i=0;i<n_u;i++)
    {
        for(int j=0;j<m_u;j++)
        {
            if(flag_u[i][j]==1)
                u[i][j]=1.0;
            else
                u[i][j]=0.0;
        }
    }
    for(int i=0;i<n_u;i++)
    {
        for(int j=0;j<m_u;j++)
        {
            u_red[i][j]=u[i][j];
            u_black[i][j]=u[i][j];
        }
    }
    for(int i=0;i<n_v;i++)
    {
        for(int j=0;j<m_v;j++)
        {
		    v[i][j]=0.0f;
		    v_red[i][j]=0.0f;
		    v_black[i][j]=0.0f;
        }
    }
    for (int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++)
        {
            pressure[i][j] = 0.0f;
            p_c[i][j]= 0.0f;
	    p_c_red[i][j]=0.0f;
            p_c_black[i][j]=0.0f;
        }
    }
    for(int i=0;i<n_u;i++)
    {
        for(int j=0;j<m_u;j++)
        {
            Bf_u[i][j]=0.0f;
        }
    }
    for(int i=0;i<n_v;i++)
    {
        for(int j=0;j<m_v;j++)
        {
            Bf_v[i][j]=0.0f;
        }
    }
    /*
     // Guess pressure field
     ifstream fileread("pressure1.txt");
     if (fileread.is_open())
     {
     for(i=0; i<n; i++)
     {
     for(j=0; j<m; j++)
     {
     fileread >> pressure[i][j];
     }
     }
     fileread.close();
     }*/

    // boundary conditions for u
    for(int i=0; i<n_u; i++)
    {
        for (int j=0; j<m_u; j++)
        {
            /*if ( i==0 ) // bottom
             {
             u[i][j]= 0.0;
             }
             else if ( i== n_u-1 ) // top
             {
             u[i][j]= 0.0;
             }*/
            if ( j==0 && i!=0 && i != n_u-1) //left
            {
                u[i][j]= 1.0f;
            }
            /*else if ( j==m_u-1 && i!=0 && i != n_u-1 ) //right
             {
             u[i][j]= 0.5;
             }*/
        }
    }
    //Boundary condition for v
    for(int i=0; i<n_v; i++)
    {
        for (int j=0; j<m_v; j++)
        {
            if ( i==0 ) // bottom
            {
                v[i][j]= 0.0f;
            }
            else if ( i== n_v-1 ) // top
            {
                v[i][j]= 0.0f;
            }
            else if ( j==0 && i!=0 && i != n_v-1) //left
            {
                v[i][j]= 0.0f;
            }
            /*else if ( j==m_v-1 && i!=0 && i != n_v-1 ) //right
             {
             v[i][j]= 0.0;
             }*/
        }
    }
    //storing the values of u and v
    for (int i=0; i<n_u; i++)
    {
        for(int j=0; j<m_u; j++)
        {
            u_old[i][j]= u[i][j];
            u_1[i][j]= u[i][j];
        }
    }
    for (int i=0; i<n_v; i++)
    {
        for(int j=0; j<m_v; j++)
        {
            v_old[i][j]= v[i][j];
            v_1[i][j]= v[i][j];
        }
    }
    // Estimation of u from x-momentum
    for(int i=0; i<n_u; i++)
    {
        for (int j=1; j<m_u; j++)
        {
            if (j==m_u-1 && i!=0 && i!=n_u-1) //right boundary
            {
                /*u[i][j] = u_old[i][j] - c_sasv*dt*((u_old[i][j]-u_old[i][j-1])/(x[0][j]-x[0][j-1]));*/
                u[i][j]=u_old[i][j-1];
            }
            else if (i==0)  //bottom boundary
            {
                u[i][j] = u_old[i+1][j];
            }
            else if (i==1 && j!=m_u-1)//bottom row
            {
                ue = ((x[i][j+1]-x1[i][j+1])/(x[i][j+1]-x[i][j]))*u_old[i][j]+((x1[i][j+1]-x[i][j])/(x[i][j+1]-x[i][j]))*u_old[i][j+1];
                uw = ((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u_old[i][j-1]+((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u_old[i][j];
                un = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i][j];
                us = u_old[i-1][j];
                vn = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j+1];
                vs = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j+1];
                u[i][j] = u_old[i][j] - (dt/(x1[i][j+1]-x1[i][j]))*(ue*ue-uw*uw) - (dt/(y[i][j]-y[i-1][j]))*(vn*un-vs*us) + (dt/((x1[i][j+1]-x1[i][j])*Re))*(((u_old[i][j+1]-u_old[i][j])/(x[i][j+1]-x[i][j]))-((u_old[i][j]-u_old[i][j-1])/(x[i][j]-x[i][j-1]))) + (dt/((y[i][j]-y[i-1][j])*Re))*(((u_old[i+1][j]-u_old[i][j])/(y1[i+1][j]-y1[i][j]))-((u_old[i][j]-u_old[i-1][j])/(y1[i][j]-y1[i-1][j]))) + (dt/(x1[i][j+1]-x1[i][j]))*(pressure[i][j]-pressure[i][j+1]);
            }
            else if (i==n_u-2 && j!=m_u-1)//top row
            {
                ue = ((x[i][j+1]-x1[i][j+1])/(x[i][j+1]-x[i][j]))*u_old[i][j]+((x1[i][j+1]-x[i][j])/(x[i][j+1]-x[i][j]))*u_old[i][j+1];
                uw = ((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u_old[i][j-1]+((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u_old[i][j];
                un = u_old[i+1][j];
                us = ((y[i-1][j]-y1[i-1][j])/(y1[i][j]-y1[i-1][j]))*u_old[i][j]+((y1[i][j]-y[i-1][j])/(y1[i][j]-y1[i-1][j]))*u_old[i-1][j];
                vn = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j+1];
                vs = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j+1];
                u[i][j] = u_old[i][j] - (dt/(x1[i][j+1]-x1[i][j]))*(ue*ue-uw*uw) - (dt/(y[i][j]-y[i-1][j]))*(vn*un-vs*us) + (dt/((x1[i][j+1]-x1[i][j])*Re))*(((u_old[i][j+1]-u_old[i][j])/(x[i][j+1]-x[i][j]))-((u_old[i][j]-u_old[i][j-1])/(x[i][j]-x[i][j-1]))) + (dt/((y[i][j]-y[i-1][j])*Re))*(((u_old[i+1][j]-u_old[i][j])/(y1[i+1][j]-y1[i][j]))-((u_old[i][j]-u_old[i-1][j])/(y1[i][j]-y1[i-1][j]))) + (dt/(x1[i][j+1]-x1[i][j]))*(pressure[i][j]-pressure[i][j+1]);
            }
            else if (i==n_u-1) //top boundary
            {
                u[i][j] = u_old[i-1][j];
            }
            else if (j!=m_u-1)
            {
                ue = ((x[i][j+1]-x1[i][j+1])/(x[i][j+1]-x[i][j]))*u_old[i][j]+((x1[i][j+1]-x[i][j])/(x[i][j+1]-x[i][j]))*u_old[i][j+1];
                uw = ((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u_old[i][j-1]+((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u_old[i][j];
                un = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i][j];
                us = ((y[i-1][j]-y1[i-1][j])/(y1[i][j]-y1[i-1][j]))*u_old[i][j]+((y1[i][j]-y[i-1][j])/(y1[i][j]-y1[i-1][j]))*u_old[i-1][j];
                vn = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j+1];
                vs = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j+1];
                u[i][j] = u_old[i][j] - (dt/(x1[i][j+1]-x1[i][j]))*(ue*ue-uw*uw) - (dt/(y[i][j]-y[i-1][j]))*(vn*un-vs*us) + (dt/((x1[i][j+1]-x1[i][j])*Re))*(((u_old[i][j+1]-u_old[i][j])/(x[i][j+1]-x[i][j]))-((u_old[i][j]-u_old[i][j-1])/(x[i][j]-x[i][j-1]))) + (dt/((y[i][j]-y[i-1][j])*Re))*(((u_old[i+1][j]-u_old[i][j])/(y1[i+1][j]-y1[i][j]))-((u_old[i][j]-u_old[i-1][j])/(y1[i][j]-y1[i-1][j]))) + (dt/(x1[i][j+1]-x1[i][j]))*(pressure[i][j]-pressure[i][j+1]);
            }
        }
    }
    //Estimation of v from y-momentum
    for(int i=1; i<n_v-1; i++)
    {
        for (int j=1; j<=m_v-1; j++)
        {
            if (j==1)//left column
            {
                ue = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i][j];
                uw = ((y[i][j-1]-y1[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i+1][j-1]+((y1[i+1][j-1]-y[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i][j-1];
                ve = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j+1];
                vw = v_old[i][j-1];
                vn = ((y1[i+1][j]-y[i][j])/(y[i+1][j]-y[i][j]))*v_old[i+1][j]+((y[i+1][j]-y1[i+1][j])/(y[i+1][j]-y[i][j]))*v_old[i][j];
                vs = ((y1[i][j]-y[i-1][j])/(y[i][j]-y[i-1][j]))*v_old[i][j]+((y[i][j]-y1[i][j])/(y[i][j]-y[i-1][j]))*v_old[i-1][j];
                v[i][j] = v_old[i][j] - (dt/(x[i][j]-x[i][j-1]))*(ue*ve-uw*vw) - (dt/(y1[i+1][j]-y1[i][j]))*(vn*vn-vs*vs) + (dt/((x[i][j]-x[i][j-1])*Re))*(((v_old[i][j+1]-v_old[i][j])/(x1[i][j+1]-x1[i][j]))-((v_old[i][j]-v_old[i][j-1])/(x1[i][j]-x1[i][j-1]))) + (dt/((y1[i+1][j]-y1[i][j])*Re))*(((v_old[i+1][j]-v_old[i][j])/(y[i+1][j]-y[i][j]))-((v_old[i][j]-v_old[i-1][j])/(y[i][j]-y[i-1][j]))) + (dt/(y1[i+1][j]-y1[i][j]))*(pressure[i][j]-pressure[i+1][j]);
            }
            else if (j==m_v-2)//right column
            {
                ue = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i][j];
                uw = ((y[i][j-1]-y1[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i+1][j-1]+((y1[i+1][j-1]-y[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i][j-1];
                ve = v_old[i][j+1];
                vw = ((x1[i][j]-x[i][j-1])/(x1[i][j]-x1[i][j-1]))*v_old[i][j-1]+((x[i][j-1]-x1[i][j-1])/(x1[i][j]-x1[i][j-1]))*v_old[i][j];
                vn = ((y1[i+1][j]-y[i][j])/(y[i+1][j]-y[i][j]))*v_old[i+1][j]+((y[i+1][j]-y1[i+1][j])/(y[i+1][j]-y[i][j]))*v_old[i][j];
                vs = ((y1[i][j]-y[i-1][j])/(y[i][j]-y[i-1][j]))*v_old[i][j]+((y[i][j]-y1[i][j])/(y[i][j]-y[i-1][j]))*v_old[i-1][j];
                v[i][j] = v_old[i][j] - (dt/(x[i][j]-x[i][j-1]))*(ue*ve-uw*vw) - (dt/(y1[i+1][j]-y1[i][j]))*(vn*vn-vs*vs) + (dt/((x[i][j]-x[i][j-1])*Re))*(((v_old[i][j+1]-v_old[i][j])/(x1[i][j+1]-x1[i][j]))-((v_old[i][j]-v_old[i][j-1])/(x1[i][j]-x1[i][j-1]))) + (dt/((y1[i+1][j]-y1[i][j])*Re))*(((v_old[i+1][j]-v_old[i][j])/(y[i+1][j]-y[i][j]))-((v_old[i][j]-v_old[i-1][j])/(y[i][j]-y[i-1][j]))) + (dt/(y1[i+1][j]-y1[i][j]))*(pressure[i][j]-pressure[i+1][j]);
            }
            else if (j==m_v-1) //right boundary
            {
                /*v[i][j] = v_old[i][j] - c_sasv*dt*((v_old[i][j]-v_old[i][j-1])/(x[0][j]-x1[0][j-1]));*/
                v[i][j] = v_old[i][j-1];
            }
            else
            {
                ue = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i][j];
                uw = ((y[i][j-1]-y1[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i+1][j-1]+((y1[i+1][j-1]-y[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i][j-1];
                ve = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j+1];
                vw = ((x1[i][j]-x[i][j-1])/(x1[i][j]-x1[i][j-1]))*v_old[i][j-1]+((x[i][j-1]-x1[i][j-1])/(x1[i][j]-x1[i][j-1]))*v_old[i][j];
                vn = ((y1[i+1][j]-y[i][j])/(y[i+1][j]-y[i][j]))*v_old[i+1][j]+((y[i+1][j]-y1[i+1][j])/(y[i+1][j]-y[i][j]))*v_old[i][j];
                vs = ((y1[i][j]-y[i-1][j])/(y[i][j]-y[i-1][j]))*v_old[i][j]+((y[i][j]-y1[i][j])/(y[i][j]-y[i-1][j]))*v_old[i-1][j];
                v[i][j] = v_old[i][j] - (dt/(x[i][j]-x[i][j-1]))*(ue*ve-uw*vw) - (dt/(y1[i+1][j]-y1[i][j]))*(vn*vn-vs*vs) + (dt/((x[i][j]-x[i][j-1])*Re))*(((v_old[i][j+1]-v_old[i][j])/(x1[i][j+1]-x1[i][j]))-((v_old[i][j]-v_old[i][j-1])/(x1[i][j]-x1[i][j-1]))) + (dt/((y1[i+1][j]-y1[i][j])*Re))*(((v_old[i+1][j]-v_old[i][j])/(y[i+1][j]-y[i][j]))-((v_old[i][j]-v_old[i-1][j])/(y[i][j]-y[i-1][j]))) + (dt/(y1[i+1][j]-y1[i][j]))*(pressure[i][j]-pressure[i+1][j]);
            }
        }
    }

    // Calculation of the coefficients of x-momentum equation
    for(int i=1; i<n_u-1; i++)
    {
        for (int j=1; j<m_u-1; j++)
        {
            if (i==1)//bottom row
            {
                b_u[i][j] = (dt/(2.0*Re*(y[i][j]-y[i-1][j])))*(1.0/(y1[i+1][j]-y1[i][j]));
                c_u[i][j] = 0.0;
                d_u[i][j] = (dt/(2.0*Re*(x1[i][j+1]-x1[i][j])))*(1.0/(x[i][j+1]-x[i][j]));
                e_u[i][j] = (dt/(2.0*Re*(x1[i][j+1]-x1[i][j])))*(1.0/(x[i][j]-x[i][j-1]));
            }
            else if (i==n_u-2)//top row
            {
                b_u[i][j] = 0.0;
                c_u[i][j] = (dt/(2.0*Re*(y[i][j]-y[i-1][j])))*(1.0/(y1[i][j]-y1[i-1][j]));
                d_u[i][j] = (dt/(2.0*Re*(x1[i][j+1]-x1[i][j])))*(1.0/(x[i][j+1]-x[i][j]));
                e_u[i][j] = (dt/(2.0*Re*(x1[i][j+1]-x1[i][j])))*(1.0/(x[i][j]-x[i][j-1]));
            }
            else
            {
                b_u[i][j] = (dt/(2.0*Re*(y[i][j]-y[i-1][j])))*(1.0/(y1[i+1][j]-y1[i][j]));
                c_u[i][j] = (dt/(2.0*Re*(y[i][j]-y[i-1][j])))*(1.0/(y1[i][j]-y1[i-1][j]));
                d_u[i][j] = (dt/(2.0*Re*(x1[i][j+1]-x1[i][j])))*(1.0/(x[i][j+1]-x[i][j]));
                e_u[i][j] = (dt/(2.0*Re*(x1[i][j+1]-x1[i][j])))*(1.0/(x[i][j]-x[i][j-1]));
            }
            a_u[i][j] = 1.0+b_u[i][j]+c_u[i][j]+d_u[i][j]+e_u[i][j];
        }
    }
    // Calculation of the co-efficients of the y-momentum equation
    for(int i=1; i<n_v-1; i++)
    {
        for (int j=1; j<m_v-1; j++)
        {
            b_v[i][j] = (dt/(2.0*Re*(y1[i+1][j]-y1[i][j])))*(1.0/(y[i+1][j]-y[i][j]));
            c_v[i][j] = (dt/(2.0*Re*(y1[i+1][j]-y1[i][j])))*(1.0/(y[i][j]-y[i-1][j]));
            d_v[i][j] = (dt/(2.0*Re*(x[i][j]-x[i][j-1])))*(1.0/(x1[i][j+1]-x1[i][j]));
            e_v[i][j] = (dt/(2.0*Re*(x[i][j]-x[i][j-1])))*(1.0/(x1[i][j]-x1[i][j-1]));
            a_v[i][j] = 1.0+b_v[i][j]+c_v[i][j]+d_v[i][j]+e_v[i][j];
        }
    }
    //Calculation of the coefficients of the pressure correction equation
    for(int i=1; i<n-1; i++)
    {
        for (int j=1; j<m-1; j++)
        {
            if (i==1 && j==1)//bottom left corner
            {
                b_p[i][j] = 1.0/((y1[i+1][j]-y1[i][j])*(y[i][j]-y[i-1][j]));
                c_p[i][j] = 0.0;
                d_p[i][j] = 1.0/((x1[i][j+1]-x1[i][j])*(x[i][j]-x[i][j-1]));
                e_p[i][j] = 0.0;
            }
            else if (i==1 && j>1 && j<m-2)//bottom row
            {
                b_p[i][j] = 1.0/((y1[i+1][j]-y1[i][j])*(y[i][j]-y[i-1][j]));
                c_p[i][j] = 0.0;
                d_p[i][j] = 1.0/((x1[i][j+1]-x1[i][j])*(x[i][j]-x[i][j-1]));
                e_p[i][j] = 1.0/((x1[i][j]-x1[i][j-1])*(x[i][j]-x[i][j-1]));
            }
            else if (i==1 && j==m-2)//bottom right corner
            {
                b_p[i][j] = 1.0/((y1[i+1][j]-y1[i][j])*(y[i][j]-y[i-1][j]));
                c_p[i][j] = 0.0;
                d_p[i][j] = 0.0;
                e_p[i][j] = 1.0/((x1[i][j]-x1[i][j-1])*(x[i][j]-x[i][j-1]));
            }
            else if (i>1 && i<n-2 && j==1)//left column
            {
                b_p[i][j] = 1.0/((y1[i+1][j]-y1[i][j])*(y[i][j]-y[i-1][j]));
                c_p[i][j] = 1.0/((y1[i][j]-y1[i-1][j])*(y[i][j]-y[i-1][j]));
                d_p[i][j] = 1.0/((x1[i][j+1]-x1[i][j])*(x[i][j]-x[i][j-1]));
                e_p[i][j] = 0.0;
            }
            else if (i>1 && i<n-2 && j==m-2)//right column
            {
                b_p[i][j] = 1.0/((y1[i+1][j]-y1[i][j])*(y[i][j]-y[i-1][j]));
                c_p[i][j] = 1.0/((y1[i][j]-y1[i-1][j])*(y[i][j]-y[i-1][j]));
                d_p[i][j] = 0.0;
                e_p[i][j] = 1.0/((x1[i][j]-x1[i][j-1])*(x[i][j]-x[i][j-1]));
            }
            else if (i==n-2 && j==1)//top left corner
            {
                b_p[i][j] = 0.0;
                c_p[i][j] = 1.0/((y1[i][j]-y1[i-1][j])*(y[i][j]-y[i-1][j]));
                d_p[i][j] = 1.0/((x1[i][j+1]-x1[i][j])*(x[i][j]-x[i][j-1]));
                e_p[i][j] = 0.0;
            }
            else if (i==n-2 && j>1 && j<m-2)//top row
            {
                b_p[i][j] = 0.0;
                c_p[i][j] = 1.0/((y1[i][j]-y1[i-1][j])*(y[i][j]-y[i-1][j]));
                d_p[i][j] = 1.0/((x1[i][j+1]-x1[i][j])*(x[i][j]-x[i][j-1]));
                e_p[i][j] = 1.0/((x1[i][j]-x1[i][j-1])*(x[i][j]-x[i][j-1]));
            }
            else if (i ==n-2 && j==m-2)//top right corner
            {
                b_p[i][j] = 0.0;
                c_p[i][j] = 1.0/((y1[i][j]-y1[i-1][j])*(y[i][j]-y[i-1][j]));
                d_p[i][j] = 0.0;
                e_p[i][j] = 1.0/((x1[i][j]-x1[i][j-1])*(x[i][j]-x[i][j-1]));
            }
            else
            {
                b_p[i][j] = 1.0/((y1[i+1][j]-y1[i][j])*(y[i][j]-y[i-1][j]));
                c_p[i][j] = 1.0/((y1[i][j]-y1[i-1][j])*(y[i][j]-y[i-1][j]));
                d_p[i][j] = 1.0/((x1[i][j+1]-x1[i][j])*(x[i][j]-x[i][j-1]));
                e_p[i][j] = 1.0/((x1[i][j]-x1[i][j-1])*(x[i][j]-x[i][j-1]));
            }
            a_p[i][j] = b_p[i][j] + c_p[i][j] + d_p[i][j] + e_p[i][j];
        }
    }
    /*ifstream fileread1xf("xf.dat");
     if (fileread1xf.is_open())
     {
	  for(int i=0; i<Nf; i++)
          {
          fileread1xf >> xf[i];
          }
          fileread1xf.close();
     }
    ifstream fileread1yf("yf.dat");
      if (fileread1yf.is_open())
      {
	  for(int i=0; i<Nf; i++)
          {
          fileread1yf >> yf[i];
          }
          fileread1yf.close();
      }
    ifstream fileread2xf("xf_old.dat");
     if (fileread2xf.is_open())
     {
	  for(int i=0; i<Nf; i++)
          {
          fileread2xf >> xf_old1[i];
          }
          fileread2xf.close();
     }
    ifstream fileread2yf("yf_old.dat");
      if (fileread2yf.is_open())
      {
	  for(int i=0; i<Nf; i++)
          {
          fileread2yf >> yf_old1[i];
          }
          fileread2yf.close();
      }
    ifstream filereadx_s("x_s.dat");
     if (filereadx_s.is_open())
     {
	  for(int i=0; i<N_marker; i++)
          {
          filereadx_s >> x_s_old[i];
          }
          filereadx_s.close();
     }
    ifstream fileready_s("y_s.dat");
      if (fileready_s.is_open())
      {
	  for(int i=0; i<N_marker; i++)
          {
          fileready_s >> y_s_old[i];
          }
          fileready_s.close();
      }
     ifstream filereadZeta("Zeta.dat");
      if (filereadZeta.is_open())
      {
	  for(int i=0; i<Nt; i++)
          {
          filereadZeta >> Zeta[i];
          }
          filereadZeta.close();
      }
     ifstream fileread1lift3f("lift3f.txt");
      if (fileread1lift3f.is_open())
      {
          for(int i=0; i<Nf; i++)
          {
          fileread1lift3f >> lift3f[i];
          }
          fileread1lift3f.close();
      }
     ifstream fileread1drag3f("drag3f.txt");
      if (fileread1drag3f.is_open())
      {
         for(int i=0; i<Nf; i++)
         {
         fileread1drag3f >> drag3f[i];
         }
         fileread1drag3f.close();
      }*/
    //ofstream file5a("forcea.txt");
    ofstream file5c("forcec.txt");
    //ofstream file6("residual_flux.txt");
    ofstream file6a("velatprobe.txt");
    //filament writing
    ofstream fileb("xf_total.dat");
    ofstream filec("yf_total.dat");
    ofstream filed("uf_total.dat");
    ofstream filee("vf_total.dat");

    ofstream filethetaf("thetaf.txt");
    ofstream filethetaf_global("thetaf_global.txt");
    //ofstream file6b("Bf_error.txt");
    for(int k=1;k<=o;k++)   //time marching
    {
        //time1 = k*dt;
		for (int kk=1;kk<=G_Str;kk++)
		{
		        for (int i=Nf1-1;i<Nf;i++)
			{
				xf_old2[i]=xf_old1[i];
				yf_old2[i]=yf_old1[i];
				xf_old1[i]=xf[i];
				yf_old1[i]=yf[i];
			}
			for (int i=Nf1-1;i<Nf;i++)
			{
				x_star[i]=2.0*xf_old1[i]-xf_old2[i];
				y_star[i]=2.0*yf_old1[i]-yf_old2[i];
			}
			for (int i=0; i<Nf1; i++)
			{
				yf[i]=0.0;
				xf[i]=xff[i];
			}	
			yf_c=amplitude_f*sin(omega_f*(time1+kk*dt1)-phi_heaving);
			v_c_rigid=amplitude_f*omega_f*cos(omega_f*(time1+kk*dt1)-phi_heaving);
			u_c_rigid=0.0;
			theta_pitching = theta_amplitude*sin(omega_f*(time1+kk*dt1)-phi_pitching);
			for (int i=0; i<Nf1; i++)
			{
				x_ss[i]=(xf[i]-radius_a)*cos(theta_pitching)+yf[i]*sin(theta_pitching)+radius_a;
				y_ss[i]=(xf[i]-radius_a)*sin(theta_pitching)+yf[i]*cos(theta_pitching)+yf_c;
			}
			for (int i=0; i<Nf1; i++)
			{
				xf[i]=x_ss[i];
				yf[i]=y_ss[i];
			}
			
			theta_pitching = -theta_pitching; // slope of the attached filament
			
			for (int i=Nf1;i<Nf;i++)
			{
				if (i==Nf1)
				{
					Fb_x[i]=(-gamma/pow(ds2,4))*(x_star[i+2]-4.0*x_star[i+1]+7.0*x_star[i]-4.0*x_star[i-1]-2.0*ds2*cos(theta_pitching));
					Fb_y[i]=(-gamma/pow(ds2,4))*(y_star[i+2]-4.0*y_star[i+1]+7.0*y_star[i]-4.0*y_star[i-1]-2.0*ds2*sin(theta_pitching));
				}
				else if (i==Nf-2)
				{
					Fb_x[i]=(-gamma/pow(ds2,4))*(-2.0*x_star[i+1]+5.0*x_star[i]-4.0*x_star[i-1]+x_star[i-2]);
					Fb_y[i]=(-gamma/pow(ds2,4))*(-2.0*y_star[i+1]+5.0*y_star[i]-4.0*y_star[i-1]+y_star[i-2]);
				}
				else if (i==Nf-1)
				{
					Fb_x[i]=(-2.0*gamma/pow(ds2,4))*(x_star[i]-2.0*x_star[i-1]+x_star[i-2]);
					Fb_y[i]=(-2.0*gamma/pow(ds2,4))*(y_star[i]-2.0*y_star[i-1]+y_star[i-2]);
				}
				else
				{
					Fb_x[i]=(-gamma/pow(ds2,4))*(x_star[i+2]-4.0*x_star[i+1]+6.0*x_star[i]-4.0*x_star[i-1]+x_star[i-2]);
					Fb_y[i]=(-gamma/pow(ds2,4))*(y_star[i+2]-4.0*y_star[i+1]+6.0*y_star[i]-4.0*y_star[i-1]+y_star[i-2]);
				}
			}
			for (int i=Nf1;i<Nt-1;i++)
			{
				A_f[i]=beta*(1.0/(2.0*pow(dt1,2)))*(1.0-(1.0/pow(ds2,2))*(2.0*pow((xf_old1[i]-xf_old1[i-1]),2)+2.0*pow((yf_old1[i]-yf_old1[i-1]),2)-pow((xf_old2[i]-xf_old2[i-1]),2)-pow((yf_old2[i]-yf_old2[i-1]),2)));
				B_f[i]=beta*(1.0/(pow(dt1,2)*pow(ds2,2)))*(pow((xf_old1[i]-xf_old2[i]-xf_old1[i-1]+xf_old2[i-1]),2) + pow((yf_old1[i]-yf_old2[i]-yf_old1[i-1]+yf_old2[i-1]),2));
				if (i==Nf1)
				{
					C_f[i]=((x_star[i]-x_star[i-1])*((Fb_x[i]+(drag3f[i]/(2.0*ds2)))+(beta*Fr*g_xg)) + (y_star[i]-y_star[i-1])*((Fb_y[i]+(lift3f[i]/(2.0*ds2)))+(beta*Fr*g_yg+beta*pow(omega_f,2)*amplitude_f*sin(omega_f*(time1+kk*dt1)-phi_heaving))))/(pow(ds2,2));
				}
				else
				{
					C_f[i]=((x_star[i]-x_star[i-1])*((Fb_x[i]+(drag3f[i]/(2.0*ds2)))-(Fb_x[i-1]+(drag3f[i-1]/(2.0*ds2)))) + (y_star[i]-y_star[i-1])*((Fb_y[i]+(lift3f[i]/(2.0*ds2)))-(Fb_y[i-1]+(lift3f[i-1]/(2.0*ds2)))))/(pow(ds2,2));
				}
				if (i==Nf1)
				{
					A_t[i]=-(1.0/(pow(ds2,4)))*(pow((x_star[i]-x_star[i-1]),2) + pow((y_star[i]-y_star[i-1]),2));
					B_t[i]=(1.0/(pow(ds2,4)))*((x_star[i]-x_star[i-1])*(x_star[i+1]-x_star[i]) + (y_star[i]-y_star[i-1])*(y_star[i+1]-y_star[i]));
					C_t[i]=0.0;
				}
				else if (i==Nt-2)
				{
					A_t[i]=-(3.0/(pow(ds2,4)))*(pow((x_star[i]-x_star[i-1]),2) + pow((y_star[i]-y_star[i-1]),2));
					B_t[i]=0.0;
					C_t[i]=(1.0/(pow(ds2,4)))*((x_star[i]-x_star[i-1])*(x_star[i-1]-x_star[i-2]) + (y_star[i]-y_star[i-1])*(y_star[i-1]-y_star[i-2]));
				}
				else
				{
					A_t[i]=-(2.0/(pow(ds2,4)))*(pow((x_star[i]-x_star[i-1]),2) + pow((y_star[i]-y_star[i-1]),2));
					B_t[i]=(1.0/(pow(ds2,4)))*((x_star[i]-x_star[i-1])*(x_star[i+1]-x_star[i]) + (y_star[i]-y_star[i-1])*(y_star[i+1]-y_star[i]));
					C_t[i]=(1.0/(pow(ds2,4)))*((x_star[i]-x_star[i-1])*(x_star[i-1]-x_star[i-2]) + (y_star[i]-y_star[i-1])*(y_star[i-1]-y_star[i-2]));
				}
			}
	// Solving Zeta (Tension coefficient) value
			for (int PT=1;PT<=T_SOR;PT++)
			{
				for (int i=Nf1-1;i<Nt;i++)
				{
				   Zeta_oldf[i]=Zeta[i];
				}
				for (int i=Nf1;i<Nt-1;i++)
				{
					Zeta[i]=alpha_SOR_f*((-B_t[i]*Zeta[i+1] - C_t[i]*Zeta[i-1] + A_f[i] - B_f[i] - C_f[i])/A_t[i])+(1.0-alpha_SOR_f)*Zeta[i];
				}
				big=0.0;
				for(int i=Nf1;i<Nt;i++)
				{
					error=abs(Zeta[i]-Zeta_oldf[i]);
					if (error>big)
					big = error;
				}
				max_errorf1[PT]=big;
				if (max_errorf1[PT] <= tolf)
				{
					break;
				}
			}
	// Solving filament positions (xf and yf)
			for (int PX=1;PX<=X_SOR;PX++)
			{
				for (int i=Nf1-1; i<Nf;i++)
				{
					xf_oldf[i]=xf[i];
					yf_oldf[i]=yf[i];
				}
				for (int i=Nf1-1;i<Nf;i++)
				{
					if (i==Nf1-1)
					{
						xf[i]=xf[i];
						yf[i]=yf[i];
					}
					else if (i==(Nf1))
					{
						xf[i]=alpha_SOR_f*(((-gamma*pow(dt1,2)/pow(ds2,4))*xf[i+2] + pow(dt1,2)*((Zeta[i+1]/pow(ds2,2))+(4.0*gamma/pow(ds2,4)))*xf[i+1] + pow(dt1,2)*((Zeta[i]/pow(ds2,2))+(4.0*gamma/pow(ds2,4)))*xf[i-1] + 2.0*cos(theta_pitching)*gamma*pow(dt1,2)/pow(ds2,3) + (beta*pow(dt1,2)*Fr*g_xg) + pow(dt1,2)*(drag3f[i]/(2.0*ds2)) + (2.0*beta*xf_old1[i]) - (beta*xf_old2[i]))/(beta+(pow((dt1/ds2),2)*(Zeta[i+1]+Zeta[i]))+(7.0*gamma*pow(dt1,2)/pow(ds2,4)))) + (1.0-alpha_SOR_f)*xf[i];

						yf[i]=alpha_SOR_f*(((-gamma*pow(dt1,2)/pow(ds2,4))*yf[i+2] + pow(dt1,2)*((Zeta[i+1]/pow(ds2,2))+(4.0*gamma/pow(ds2,4)))*yf[i+1] + pow(dt1,2)*((Zeta[i]/pow(ds2,2))+(4.0*gamma/pow(ds2,4)))*yf[i-1] + 2.0*sin(theta_pitching)*gamma*pow(dt1,2)/pow(ds2,3) + (beta*pow(dt1,2)*Fr*g_yg) + pow(dt1,2)*(lift3f[i]/(2.0*ds2)) + (2.0*beta*yf_old1[i]) - (beta*yf_old2[i]))/(beta+(pow((dt1/ds2),2)*(Zeta[i+1]+Zeta[i]))+(7.0*gamma*pow(dt1,2)/pow(ds2,4)))) + (1.0-alpha_SOR_f)*yf[i];
					}
					else if (i==Nf-2)
					{
						xf[i]=alpha_SOR_f*((pow(dt1,2)*((Zeta[i+1]/pow(ds2,2))+(2.0*gamma/pow(ds2,4)))*xf[i+1] + pow(dt1,2)*((Zeta[i]/pow(ds2,2))+(4.0*gamma/pow(ds2,4)))*xf[i-1] - (gamma*pow(dt1,2)/pow(ds2,4))*xf[i-2] + (beta*pow(dt1,2)*Fr*g_xg) + pow(dt1,2)*(drag3f[i]/(2.0*ds2)) + (2.0*beta*xf_old1[i]) - (beta*xf_old2[i]))/(beta+(pow((dt1/ds2),2)*(Zeta[i+1]+Zeta[i]))+(5.0*gamma*pow(dt1,2)/pow(ds2,4)))) + (1.0-alpha_SOR_f)*xf[i];

						yf[i]=alpha_SOR_f*((pow(dt1,2)*((Zeta[i+1]/pow(ds2,2))+(2.0*gamma/pow(ds2,4)))*yf[i+1] + pow(dt1,2)*((Zeta[i]/pow(ds2,2))+(4.0*gamma/pow(ds2,4)))*yf[i-1] - (gamma*pow(dt1,2)/pow(ds2,4))*yf[i-2] + (beta*pow(dt1,2)*Fr*g_yg) + pow(dt1,2)*(lift3f[i]/(2.0*ds2)) + (2.0*beta*yf_old1[i]) - (beta*yf_old2[i]))/(beta+(pow((dt1/ds2),2)*(Zeta[i+1]+Zeta[i]))+(5.0*gamma*pow(dt1,2)/pow(ds2,4)))) + (1.0-alpha_SOR_f)*yf[i];
					}
					else if (i==Nf-1)
					{
						xf[i]=alpha_SOR_f*((pow(dt1,2)*((2.0*Zeta[i]/pow(ds2,2))+(4.0*gamma/pow(ds2,4)))*xf[i-1] - (2.0*gamma*pow(dt1,2)/pow(ds2,4))*xf[i-2] + (beta*pow(dt1,2)*Fr*g_xg) + pow(dt1,2)*(drag3f[i]/(2.0*ds2)) + (2.0*beta*xf_old1[i]) - (beta*xf_old2[i]))/(beta+(2.0*pow((dt1/ds2),2)*Zeta[i])+(2.0*gamma*pow(dt1,2)/pow(ds2,4)))) + (1.0-alpha_SOR_f)*xf[i];

						yf[i]=alpha_SOR_f*((pow(dt1,2)*((2.0*Zeta[i]/pow(ds2,2))+(4.0*gamma/pow(ds2,4)))*yf[i-1] - (2.0*gamma*pow(dt1,2)/pow(ds2,4))*yf[i-2] + (beta*pow(dt1,2)*Fr*g_yg) + pow(dt1,2)*(lift3f[i]/(2.0*ds2)) + (2.0*beta*yf_old1[i]) - (beta*yf_old2[i]))/(beta+(2.0*pow((dt1/ds2),2)*Zeta[i])+(2.0*gamma*pow(dt1,2)/pow(ds2,4)))) + (1.0-alpha_SOR_f)*yf[i];
					}
					else
					{
						xf[i]=alpha_SOR_f*(((-gamma*pow(dt1,2)/pow(ds2,4))*xf[i+2] + pow(dt1,2)*((Zeta[i+1]/pow(ds2,2))+(4.0*gamma/pow(ds2,4)))*xf[i+1] + pow(dt1,2)*((Zeta[i]/pow(ds2,2))+(4.0*gamma/pow(ds2,4)))*xf[i-1] - (gamma*pow(dt1,2)/pow(ds2,4))*xf[i-2] + (beta*pow(dt1,2)*Fr*g_xg) + pow(dt1,2)*(drag3f[i]/(2.0*ds2)) + (2.0*beta*xf_old1[i]) - (beta*xf_old2[i]))/(beta+(pow((dt1/ds2),2)*(Zeta[i+1]+Zeta[i]))+(6.0*gamma*pow(dt1,2)/pow(ds2,4)))) + (1.0-alpha_SOR_f)*xf[i];

						yf[i]=alpha_SOR_f*(((-gamma*pow(dt1,2)/pow(ds2,4))*yf[i+2] + pow(dt1,2)*((Zeta[i+1]/pow(ds2,2))+(4.0*gamma/pow(ds2,4)))*yf[i+1] + pow(dt1,2)*((Zeta[i]/pow(ds2,2))+(4.0*gamma/pow(ds2,4)))*yf[i-1] - (gamma*pow(dt1,2)/pow(ds2,4))*yf[i-2] + (beta*pow(dt1,2)*Fr*g_yg) + pow(dt1,2)*(lift3f[i]/(2.0*ds2)) + (2.0*beta*yf_old1[i]) - (beta*yf_old2[i]))/(beta+(pow((dt1/ds2),2)*(Zeta[i+1]+Zeta[i]))+(6.0*gamma*pow(dt1,2)/pow(ds2,4)))) + (1.0-alpha_SOR_f)*yf[i];
					}
				}
				big=0.0;
				for(int i=Nf-1;i<Nf;i++)
				{
					error=abs(xf[i]-xf_oldf[i]);
					if (error>big)
					big = error;
				}
				max_errorf2[PX]=big;
				big=0.0;
				for(int i=Nf1-1;i<Nf;i++)
				{
				   error=abs(yf[i]-yf_oldf[i]);
				   if (error>big)
				   big = error;
				}
				max_errorf3[PX]=big;
				if (max_errorf2[PX] <= tolf && max_errorf3[PX] <= tolf)
				{
					break;
				}
			}
		}
		time1 = k*dt;
        // To find the center nodes along the chord length
        for (int i=0;i<Nt;i++)
        {
            if(i==0)
            {
                xfc[i]=xf[i];
                yfc[i]=yf[i];
            }
            else if(i==Nt-1)
            {
                xfc[i]=xf[i-1];
                yfc[i]=yf[i-1];
            }
            else
            {
                xfc[i]=(xf[i-1]+xf[i])/2.0;
                yfc[i]=(yf[i-1]+yf[i])/2.0;
            }
        }
        // local slope
/*        for (int i=0; i<(Nf-1)/2+1; i++)
        {
        	thetaf[i]=theta_amplitude*sin(omega_f*k*dt);
        }
*/        
        for (int i=0; i<Nf; i++)
        {
            if (i==0)
            {
               thetaf[i]=atan((yf[i+1]-yf[i])/(xf[i+1]-xf[i]));
            }
            else if (i==Nf1-1)
            {
                thetaf[i]=atan((yf[i]-yf[i-1])/(xf[i]-xf[i-1]));
            }
            else if (i==Nf-1)
            {
                thetaf[i]=atan((yf[i]-yf[i-1])/(xf[i]-xf[i-1]));
            }
            else
            {
                thetaf[i]=atan((yf[i+1]-yf[i-1])/(xf[i+1]-xf[i-1]));
            }
        }
        // to find the surface points
        // upper surface
        for (int i=0; i<Nf; i++)
        {
	    	if(i==0)
	    	{
				x_sf[i]=xf[i];
				y_sf[i]=yf[i];
	    	}
			else if (i==Nf-1 || i==Nf1-1)
			{
				if ((xf[i]-xf[i-1])>=0.0 && thetaf[i]>=0.0)
				{
					x_sf[i]=xf[i] - hf[i]*sin(thetaf[i]);
					y_sf[i]=yf[i] + hf[i]*cos(thetaf[i]);
				}
				else if ((xf[i]-xf[i-1])<0.0 && thetaf[i]<=0.0)
				{
					x_sf[i]=xf[i] + hf[i]*sin(thetaf[i]);
					y_sf[i]=yf[i] - hf[i]*cos(thetaf[i]);
				}
				else if ((xf[i]-xf[i-1])>=0.0 && thetaf[i]<=0.0)
				{
					x_sf[i]=xf[i] - hf[i]*sin(thetaf[i]);
					y_sf[i]=yf[i] + hf[i]*cos(thetaf[i]);
				}
				else if ((xf[i]-xf[i-1])<0.0 && thetaf[i]>=0.0)
				{
					x_sf[i]=xf[i] + hf[i]*sin(thetaf[i]);
					y_sf[i]=yf[i] - hf[i]*cos(thetaf[i]);
				}
				else
				{
					cout<< " Top surface end " <<endl;
				}
			}
			else
			{
				if ((xf[i+1]-xf[i-1])>=0.0 && thetaf[i]>=0.0)
				{
					x_sf[i]=xf[i] - hf[i]*sin(thetaf[i]);
					y_sf[i]=yf[i] + hf[i]*cos(thetaf[i]);
				}
				else if ((xf[i+1]-xf[i-1])<0.0 && thetaf[i]<=0.0)
				{
					x_sf[i]=xf[i] + hf[i]*sin(thetaf[i]);
					y_sf[i]=yf[i] - hf[i]*cos(thetaf[i]);
				}
				else if ((xf[i+1]-xf[i-1])>=0.0 && thetaf[i]<=0.0)
				{
					x_sf[i]=xf[i] - hf[i]*sin(thetaf[i]);
					y_sf[i]=yf[i] + hf[i]*cos(thetaf[i]);
				}
				else if ((xf[i+1]-xf[i-1])<0.0 && thetaf[i]>=0.0)
				{
					x_sf[i]=xf[i] + hf[i]*sin(thetaf[i]);
					y_sf[i]=yf[i] - hf[i]*cos(thetaf[i]);
				}
				else
				{
					cout<< " Top surface middle " << " " << i <<endl;
				}
			}
        }
        // for lower surface
        for (int i=Nf; i<N_markerf; i++)
        {
			if(i==Nf)
			{
				x_sf[i]=xf[N_markerf-i-1];
				y_sf[i]=yf[N_markerf-i-1];
			}
			else if (i==Nf+1 || i==2*Nf2+Nf1-1)
			{
				if ((xf[N_markerf-i]-xf[N_markerf-i-1])>=0.0 && thetaf[N_markerf-i]>=0.0)
				{
					x_sf[i]=xf[N_markerf-i] + hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
					y_sf[i]=yf[N_markerf-i] - hf[N_markerf-i]*cos(thetaf[N_markerf-i]);
				}
				else if ((xf[N_markerf-i]-xf[N_markerf-i-1])<0.0 && thetaf[N_markerf-i]<=0.0)
				{
					x_sf[i]=xf[N_markerf-i] - hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
					y_sf[i]=yf[N_markerf-i] + hf[N_markerf-i]*cos(thetaf[N_markerf-i-1]);
				}
				else if ((xf[N_markerf-i]-xf[N_markerf-i-1])>=0.0 && thetaf[N_markerf-i]<=0.0)
				{
					x_sf[i]=xf[N_markerf-i] + hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
					y_sf[i]=yf[N_markerf-i] - hf[N_markerf-i]*cos(thetaf[N_markerf-i]);
				}
				else if ((xf[N_markerf-i]-xf[N_markerf-i-1])<0.0 && thetaf[N_markerf-i]>=0.0)
				{
					x_sf[i]=xf[N_markerf-i] - hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
					y_sf[i]=yf[N_markerf-i] + hf[N_markerf-i]*cos(thetaf[N_markerf-i]);
				}
				else
				{
					cout<< " Bottom surface end " <<endl;
				}
			}
			else
			{
				if ((xf[N_markerf-i+1]-xf[N_markerf-i-1])>=0.0 && thetaf[N_markerf-i]>=0.0)
				{
					x_sf[i]=xf[N_markerf-i] + hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
					y_sf[i]=yf[N_markerf-i] - hf[N_markerf-i]*cos(thetaf[N_markerf-i]);
				}
				else if ((xf[N_markerf-i+1]-xf[N_markerf-i-1])<0.0 && thetaf[N_markerf-i]<=0.0)
				{
					x_sf[i]=xf[N_markerf-i] - hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
					y_sf[i]=yf[N_markerf-i] + hf[N_markerf-i]*cos(thetaf[N_markerf-i]);
				}
				else if ((xf[N_markerf-i+1]-xf[N_markerf-i-1])>=0.0 && thetaf[N_markerf-i]<=0.0)
				{
					x_sf[i]=xf[N_markerf-i] + hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
					y_sf[i]=yf[N_markerf-i] - hf[N_markerf-i]*cos(thetaf[N_markerf-i]);
				}
				else if ((xf[N_markerf-i+1]-xf[N_markerf-i-1])<0.0 && thetaf[N_markerf-i]>=0.0)
				{
					x_sf[i]=xf[N_markerf-i] - hf[N_markerf-i]*sin(thetaf[N_markerf-i]);
					y_sf[i]=yf[N_markerf-i] + hf[N_markerf-i]*cos(thetaf[N_markerf-i]);
				}
				else
				{
					cout<< " Bottom surface middle " << " " << i <<endl;
				}
			}
        }

		// Generating surface points along the centre points
		// local slope along center points of filamment
	/*	for (int i=0; i<(Nf-1)/2+1; i++)
		{
			thetafc[i]=theta_amplitude*sin(omega_f*k*dt);
		}
	*/
        for (int i=0; i<Nt; i++)
        {
            if (i==0)
            {
               thetafc[i]=atan((yfc[i+1]-yfc[i])/(xfc[i+1]-xfc[i]));
            }
            else if (i==Nt-1)
            {
                thetafc[i]=atan((yfc[i]-yfc[i-1])/(xfc[i]-xfc[i-1]));
            }
            else
            {
                thetafc[i]=atan((yfc[i+1]-yfc[i-1])/(xfc[i+1]-xfc[i-1]));
            }
        }
		// upper surface
        for (int i=0; i<Nf+1; i++)
        {
            if (i==0)
            {
                x_sfc[i]=xfc[i];
		y_sfc[i]=yfc[i];
            }
			else if (i==Nf)
			{
				if ((xfc[i]-xfc[i-1])>=0.0 && thetafc[i]>=0.0)
				{
					x_sfc[i]=xfc[i] - hf[i]*sin(thetafc[i]);
					y_sfc[i]=yfc[i] + hf[i]*cos(thetafc[i]);
				}
				else if ((xfc[i]-xfc[i-1])<0.0 && thetafc[i]<=0.0)
				{
					x_sfc[i]=xfc[i] + hf[i]*sin(thetafc[i]);
					y_sfc[i]=yfc[i] - hf[i]*cos(thetafc[i]);
				}
				else if ((xfc[i]-xfc[i-1])>=0.0 && thetafc[i]<=0.0)
				{
					x_sfc[i]=xfc[i] - hf[i]*sin(thetafc[i]);
					y_sfc[i]=yfc[i] + hf[i]*cos(thetafc[i]);
				}
				else if ((xfc[i]-xfc[i-1])<0.0 && thetafc[i]>=0.0)
				{
					x_sfc[i]=xfc[i] + hf[i]*sin(thetafc[i]);
					y_sfc[i]=yfc[i] - hf[i]*cos(thetafc[i]);
				}
			}
			else
			{
				if ((xfc[i+1]-xfc[i-1])>=0.0 && thetafc[i]>=0.0)
				{
					x_sfc[i]=xfc[i] - hf[i]*sin(thetafc[i]);
					y_sfc[i]=yfc[i] + hf[i]*cos(thetafc[i]);
				}
				else if ((xfc[i+1]-xfc[i-1])<0.0 && thetafc[i]<=0.0)
				{
					x_sfc[i]=xfc[i] + hf[i]*sin(thetafc[i]);
					y_sfc[i]=yfc[i] - hf[i]*cos(thetafc[i]);
				}
				else if ((xfc[i+1]-xfc[i-1])>=0.0 && thetafc[i]<=0.0)
				{
					x_sfc[i]=xfc[i] - hf[i]*sin(thetafc[i]);
					y_sfc[i]=yfc[i] + hf[i]*cos(thetafc[i]);
				}
				else if ((xfc[i+1]-xfc[i-1])<0.0 && thetafc[i]>=0.0)
				{
					x_sfc[i]=xfc[i] + hf[i]*sin(thetafc[i]);
					y_sfc[i]=yfc[i] - hf[i]*cos(thetafc[i]);
				}
			}
        }
        // for lower surface
        for (int i=Nf+1; i<N_markerfc; i++)
        {
            if (i==Nf+1)
			{
				if ((xfc[N_markerf-i]-xfc[N_markerf-i-1])>=0.0 && thetafc[N_markerf-i]>=0.0)
				{
					x_sfc[i]=xfc[N_markerf-i] + hf[N_markerf-i]*sin(thetafc[N_markerf-i]);
					y_sfc[i]=yfc[N_markerf-i] - hf[N_markerf-i]*cos(thetafc[N_markerf-i]);
				}
				else if ((xfc[N_markerf-i]-xfc[N_markerf-i-1])<0.0 && thetafc[N_markerf-i]<=0.0)
				{
					x_sfc[i]=xfc[N_markerf-i] - hf[N_markerf-i]*sin(thetafc[N_markerf-i]);
					y_sfc[i]=yfc[N_markerf-i] + hf[N_markerf-i]*cos(thetafc[N_markerf-i]);
				}
				else if ((xfc[N_markerf-i]-xfc[N_markerf-i-1])>=0.0 && thetafc[N_markerf-i]<=0.0)
				{
					x_sfc[i]=xfc[N_markerf-i] + hf[N_markerf-i]*sin(thetafc[N_markerf-i]);
					y_sfc[i]=yfc[N_markerf-i] - hf[N_markerf-i]*cos(thetafc[N_markerf-i]);
				}
				else if ((xfc[N_markerf-i]-xfc[N_markerf-i-1])<0.0 && thetafc[N_markerf-i]>=0.0)
				{
					x_sfc[i]=xfc[N_markerf-i] - hf[N_markerf-i]*sin(thetafc[N_markerf-i]);
					y_sfc[i]=yfc[N_markerf-i] + hf[N_markerf-i]*cos(thetafc[N_markerf-i]);
				}
			}
			else
			{
				if ((xfc[N_markerfc-i+1]-xfc[N_markerfc-i-1])>=0.0 && thetafc[N_markerfc-i]>=0.0)
				{
					x_sfc[i]=xfc[N_markerfc-i] + hf[N_markerfc-i]*sin(thetafc[N_markerfc-i]);
					y_sfc[i]=yfc[N_markerfc-i] - hf[N_markerfc-i]*cos(thetafc[N_markerfc-i]);
				}
				else if ((xfc[N_markerfc-i+1]-xfc[N_markerfc-i-1])<0.0 && thetafc[N_markerfc-i]<=0.0)
				{
					x_sfc[i]=xfc[N_markerfc-i] - hf[N_markerfc-i]*sin(thetafc[N_markerfc-i]);
					y_sfc[i]=yfc[N_markerfc-i] + hf[N_markerfc-i]*cos(thetafc[N_markerfc-i]);
				}
				else if ((xfc[N_markerfc-i+1]-xfc[N_markerfc-i-1])>=0.0 && thetafc[N_markerfc-i]<=0.0)
				{
					x_sfc[i]=xfc[N_markerfc-i] + hf[N_markerfc-i]*sin(thetafc[N_markerfc-i]);
					y_sfc[i]=yfc[N_markerfc-i] - hf[N_markerfc-i]*cos(thetafc[N_markerfc-i]);
				}
				else if ((xfc[N_markerfc-i+1]-xfc[N_markerfc-i-1])<0.0 && thetafc[N_markerfc-i]>=0.0)
				{
					x_sfc[i]=xfc[N_markerfc-i] - hf[N_markerfc-i]*sin(thetafc[N_markerfc-i]);
					y_sfc[i]=yfc[N_markerfc-i] + hf[N_markerfc-i]*cos(thetafc[N_markerfc-i]);
				}
			}
        }
        x_s_max=x_sf[0];
        y_s_max=y_sf[0];
        x_s_min=x_sf[0];
        y_s_min=y_sf[0];
        for (int i=0; i<N_markerf; i++)
        {
	        if(x_s_max<x_sf[i])
	        {
		        x_s_max=x_sf[i];
	        }
	        if(y_s_max<y_sf[i])
	        {
		        y_s_max=y_sf[i];
	        }
        }
        for (int i=0; i<N_markerf; i++)
        {
	        if(x_s_min>x_sf[i])
	        {
		        x_s_min=x_sf[i];
	        }
	        if(y_s_min>y_sf[i])
	        {
		        y_s_min=y_sf[i];
	        }
        }

        // Boundary for flagging purpose
        m_div11 = m_div1+ceil((x_s_min-Hx1)/dx2);
        m_div22 = ceil((x_s_max-x_s_min)/dx2);
        n_div11 = n_div1+ceil((y_s_min-Hy1)/dy2);
        n_div22 = ceil((y_s_max-y_s_min)/dy2);

        c11 = m_div11-20;
        c22 = m_div11+m_div22+20;
        c33 = n_div11-20;
        c44 = n_div11+n_div22+20;
        // co-ordinate to compute force only for the rigid foil
        x_s_max=x_sf[0];
        y_s_max=y_sf[0];
        x_s_min=x_sf[0];
        y_s_min=y_sf[0];
        for (int i=0; i<Nf1; i++)
        {
		if(x_s_max<x_sf[i])
		{
		        x_s_max=x_sf[i];
		}
		if(y_s_max<y_sf[i])
		{
		        y_s_max=y_sf[i];
		}
        }
	for (int i=Nf1+N_markerf2-2; i<N_markerf; i++)
        {
		if(x_s_max<x_sf[i])
		{
		        x_s_max=x_sf[i];
		}
		if(y_s_max<y_sf[i])
		{
		        y_s_max=y_sf[i];
		}
        }
	for (int i=0; i<Nf1; i++)
        {
		if(x_s_min>x_sf[i])
		{
		        x_s_min=x_sf[i];
		}
		if(y_s_min>y_sf[i])
		{
		        y_s_min=y_sf[i];
		}
        }        
	for (int i=Nf1+N_markerf2-2; i<N_markerf; i++)
        {
		if(x_s_min>x_sf[i])
		{
		        x_s_min=x_sf[i];
		}
		if(y_s_min>y_sf[i])
		{
		        y_s_min=y_sf[i];
		}
        }

        // Boundary for flagging purpose
        m_div11 = m_div1+ceil((x_s_min-Hx1)/dx2);
        m_div22 = ceil((x_s_max-x_s_min)/dx2);
        n_div11 = n_div1+ceil((y_s_min-Hy1)/dy2);
        n_div22 = ceil((y_s_max-y_s_min)/dy2);

        c11_f = m_div11-20;
        c22_f = m_div11+m_div22-1;
        c33_f = n_div11-4;
        c44_f = n_div11+n_div22+4;
        //Generating surface points for flagging using linear interpolation
    for (int i=0; i<(Nf1*nn1)-(nn1-1); i++)
    {
		ctr1=trunc(i/nn1);
		if(i%nn1==0)
		{
			x_s[i]=x_sf[ctr1];
			y_s[i]=y_sf[ctr1];
		}
		else
		{
			dx = (x_sf[ctr1+1]-x_sf[ctr1])/(nn1);
			dy = (y_sf[ctr1+1]-y_sf[ctr1])/(nn1);
			if (atan((y_sf[ctr1+1]-y_sf[ctr1])/(x_sf[ctr1+1]-x_sf[ctr1]))==api || atan((y_sf[ctr1+1]-y_sf[ctr1])/(x_sf[ctr1+1]-x_sf[ctr1]))==-api || dx==0.0 || dx==-0.0)
			{
				x_s[i] = x_sf[ctr1];
				y_s[i] = y_sf[ctr1] + (i-ctr1*nn1)*dy;
			}
			else
			{
				x_s[i] = x_sf[ctr1] + (i-ctr1*nn1)*dx;
				y_s[i] = y_sf[ctr1] + (y_sf[ctr1+1]-y_sf[ctr1])/(x_sf[ctr1+1]-x_sf[ctr1])*(x_s[i]-x_sf[ctr1]);
			}
		}
    }
    counter_f=1;
    for (int i=(Nf1*nn1)-(nn1-1); i<((Nf2*2)*nn2+(Nf1*nn1)-(nn1-1)); i++)
    {
		ctr1=trunc((i+1-((Nf1*nn1)-(nn1-1)))/nn2);
		if((i+1-((Nf1*nn1)-(nn1-1)))%nn2==0)
		{
			x_s[i]=x_sf[(Nf1-1) + ctr1];
			y_s[i]=y_sf[(Nf1-1) + ctr1];
		}
		else
		{
			
			dx = (x_sf[(Nf1-1) + ctr1+1]-x_sf[(Nf1-1) + ctr1])/(nn2);
			dy = (y_sf[(Nf1-1) + ctr1+1]-y_sf[(Nf1-1) + ctr1])/(nn2);
			if (atan((y_sf[(Nf1-1) + ctr1+1]-y_sf[(Nf1-1) + ctr1])/(x_sf[(Nf1-1) + ctr1+1]-x_sf[(Nf1-1) + ctr1]))==api || atan((y_sf[(Nf1-1) + ctr1+1]-y_sf[(Nf1-1) + ctr1])/(x_sf[(Nf1-1) + ctr1+1]-x_sf[(Nf1-1) + ctr1]))==-api || dx==0.0 || dx==-0.0)
			{
				x_s[i] = x_sf[(Nf1-1) + ctr1];
				y_s[i] = y_sf[(Nf1-1) + ctr1] + (counter_f-ctr1*nn2)*dy;
			}
			else
			{
				x_s[i] = x_sf[(Nf1-1) + ctr1] + (counter_f-ctr1*nn2)*dx;
				y_s[i] = y_sf[(Nf1-1) + ctr1] + (y_sf[(Nf1-1) + ctr1+1]-y_sf[(Nf1-1) + ctr1])/(x_sf[(Nf1-1) + ctr1+1]-x_sf[(Nf1-1) + ctr1])*(x_s[i]-x_sf[(Nf1-1) + ctr1]);
			}
		}
		counter_f=counter_f+1;
    }
    counter_f=1;
    for (int i=((Nf2*2)*nn2+(Nf1*nn1)-(nn1-1)); i<N_marker; i++)
    {
		ctr1=trunc((i+1-((Nf2*2)*nn2+(Nf1*nn1)-(nn1-1)))/nn1);
		if((i+1-((Nf2*2)*nn2+(Nf1*nn1)-(nn1-1)))%nn1==0)
		{
			x_s[i]=x_sf[(Nf1+Nf2*2-1) + ctr1];
			y_s[i]=y_sf[(Nf1+Nf2*2-1) + ctr1];
		}
		else
		{
			if (((Nf1+Nf2*2-1) + ctr1)==N_markerf-1)
			{
				dx = (x_sf[0]-x_sf[(Nf1+Nf2*2-1) + ctr1])/(nn1);
				dy = (y_sf[0]-y_sf[(Nf1+Nf2*2-1) + ctr1])/(nn1);
			}
			else
			{
				dx = (x_sf[(Nf1+Nf2*2-1) + ctr1+1]-x_sf[(Nf1+Nf2*2-1) + ctr1])/(nn1);
				dy = (y_sf[(Nf1+Nf2*2-1) + ctr1+1]-y_sf[(Nf1+Nf2*2-1) + ctr1])/(nn1);
			}
			if (atan((y_sf[(Nf1+Nf2*2-1) + ctr1+1]-y_sf[(Nf1+Nf2*2-1) + ctr1])/(x_sf[(Nf1+Nf2*2-1) + ctr1+1]-x_sf[(Nf1+Nf2*2-1) + ctr1]))==api || atan((y_sf[(Nf1+Nf2*2-1) + ctr1+1]-y_sf[(Nf1+Nf2*2-1) + ctr1])/(x_sf[(Nf1+Nf2*2-1) + ctr1+1]-x_sf[(Nf1+Nf2*2-1) + ctr1]))==-api || dx==0.0 || dx==-0.0)
			{
				x_s[i] = x_sf[(Nf1+Nf2*2-1) + ctr1];
				y_s[i] = y_sf[(Nf1+Nf2*2-1) + ctr1] + (counter_f-ctr1*nn1)*dy;
			}
			else
			{
				x_s[i] = x_sf[(Nf1+Nf2*2-1) + ctr1] + (counter_f-ctr1*nn1)*dx;
				if (((Nf1+Nf2*2-1) + ctr1)==N_markerf-1)
				{
					y_s[i] = y_sf[(Nf1+Nf2*2-1) + ctr1] + (y_sf[0]-y_sf[(Nf1+Nf2*2-1) + ctr1])/(x_sf[0]-x_sf[(Nf1+Nf2*2-1) + ctr1])*(x_s[i]-x_sf[(Nf1+Nf2*2-1) + ctr1]);
				}
				else
				{
					y_s[i] = y_sf[(Nf1+Nf2*2-1) + ctr1] + (y_sf[(Nf1+Nf2*2-1) + ctr1+1]-y_sf[(Nf1+Nf2*2-1) + ctr1])/(x_sf[(Nf1+Nf2*2-1) + ctr1+1]-x_sf[(Nf1+Nf2*2-1) + ctr1])*(x_s[i]-x_sf[(Nf1+Nf2*2-1) + ctr1]);
				}
			}
		}
		counter_f=counter_f+1;
    }      
        for (int i=0; i<N_marker; i++)
        {
            u_c[i]=(x_s[i]-x_s_old[i])/dt; // velocity of surface point
            v_c[i]=(y_s[i]-y_s_old[i])/dt;
        }
        for (int i = 0;i<N_marker;i++)
        {
            if(i==0)
            {
                v1 = x_s[i+1] - x_s[N_marker-1];
                v2 = y_s[i+1] - y_s[N_marker-1];
            }
            else if (i==N_marker-1)
            {
                v1 = x_s[0] - x_s[i-1];
                v2 = y_s[0] - y_s[i-1];
            }
            else
            {
                v1 = x_s[i+1] - x_s[i-1];
                v2 = y_s[i+1] - y_s[i-1];
            }
            n1[i] = -v2;
            n2[i] = v1;
            n1[i] = (-v2)/sqrt((-v2)*(-v2)+v1*v1);
            n2[i] = v1/sqrt((-v2)*(-v2)+v1*v1);
            //t1[i] = v1/sqrt(v1*v1+v2*v2);
            //t2[i] = v2/sqrt(v1*v1+v2*v2);
        }
        if ( (k%writeInterval) == 0)
        {
            string fileName;
            stringstream ss;
            ss << k;
            fileName = ss.str();
            fileName = "solidBoundary"+fileName+".dat";
            ofstream file7(fileName.c_str());
            if (file7.is_open())
            {
                for(int i=0; i<N_marker; i++)
                {
                    file7 << x_s[i] << " " << y_s[i] << " " << n1[i] << " " << n2[i] << endl;
                }
                file7 << x_s[0] << " " << y_s[0] << " " << n1[0] << " " << n2[0];
                file7.close();
            }
            fileName = ss.str();
            fileName = "xf"+fileName+".dat";
            ofstream file1xf(fileName.c_str());
            if (file1xf.is_open())
            {
                for(int i=0; i<Nf; i++)
                {
                    file1xf << xf[i] << endl;
                }
                file1xf.close();
            }
            fileName = ss.str();
            fileName = "xf_old"+fileName+".dat";
            ofstream file2xf(fileName.c_str());
            if (file2xf.is_open())
            {
                for(int i=0; i<Nf; i++)
                {
                    file2xf << xf_old1[i] << endl;
                }
                file2xf.close();
            }
            fileName = ss.str();
            fileName = "yf"+fileName+".dat";
            ofstream file1yf(fileName.c_str());
            if (file1yf.is_open())
            {
                for(int i=0; i<Nf; i++)
                {
                    file1yf << yf[i] << endl;
                }
                file1yf.close();
            }
            fileName = ss.str();
            fileName = "yf_old"+fileName+".dat";
            ofstream file2yf(fileName.c_str());
            if (file2yf.is_open())
            {
                for(int i=0; i<Nf; i++)
                {
                    file2yf << yf_old1[i] << endl;
                }
                file2yf.close();
            }
            fileName = ss.str();
            fileName = "x_s"+fileName+".dat";
            ofstream filex_s(fileName.c_str());
            if (filex_s.is_open())
            {
                for(int i=0; i<N_marker; i++)
                {
                    filex_s << x_s[i] << endl;
                }
                filex_s.close();
            }
            fileName = ss.str();
            fileName = "y_s"+fileName+".dat";
            ofstream filey_s(fileName.c_str());
            if (filey_s.is_open())
            {
                for(int i=0; i<N_marker; i++)
                {
                    filey_s << y_s[i] << endl;
                }
                filey_s.close();
            }
            fileName = ss.str();
            fileName = "Zeta"+fileName+".dat";
            ofstream fileZeta(fileName.c_str());
            if (fileZeta.is_open())
            {
                for(int i=0; i<Nt; i++)
                {
                    fileZeta << Zeta[i] << endl;
                }
                fileZeta.close();
            }
        }
        if (fileb.is_open())
	{
	    fileb << fixed << setprecision(6) << time1 << " " ;
	    for(int i=0; i<Nf; i++)
	    {
	        fileb << xf[i] << " " ;
	    }
            fileb<< " " <<endl;
	    //file << xf[0][0] << " " << yf[0][0] ;
            //fileb.close();
        }
	if (filec.is_open())
	{
            filec << fixed << setprecision(6) << time1 << " " ;
            for(int i=0; i<Nf; i++)
	    {
                filec << yf[i] << " " ;
	    }
	    filec<< " " <<endl;
	    //file << xf[0][0] << " " << yf[0][0] ;
            //filec.close();
	}
        /*if (filethetaf.is_open())
	{
            for(int i=0; i<Nf; i++)
	    {
                filethetaf << thetaf[i] << " " ;
	    }
	    filethetaf<< " " <<endl;
	}
        if (filethetaf_global.is_open())
	{
            for(int i=0; i<Nf; i++)
	    {
                filethetaf_global << thetaf_global[i] << " " ;
	    }
	    filethetaf_global<< " " <<endl;
	}*/
        for (int i=0; i<N_marker; i++)
        {
            x_s_old[i]=x_s[i];
            y_s_old[i]=y_s[i];
        }
        /*c1 = m_div1 + (((xc-Hx1)-(radius_a+10*dx2))/dx2);
         c2 = m_div1+m_div2-(((Hx2-xc)-(radius_a+10*dx2))/dx2);
         c3 = n_div1 + (((yc-Hy1)-(radius_b+10*dy2))/dy2);
         c4 = n_div1+n_div2-(((Hy2-yc)-(radius_b+10*dy2))/dy2);*/
        c1 = m_div1;
        //c22_f=m_div1+ceil((Lf1-Hx1)/dx2);
        c2 = m_div1+m_div2;
        c3 = n_div1;
        c4 = n_div1+n_div2;
	//cout<< " c22 = " << c22 << " c22_f = " << c22_f <<endl;
        //#pragma omp parallel for default(shared) private(min_dis,distance,k1,temp,sig,phi) schedule(dynamic)
#pragma acc data copy(n1, n2) copy(flag_u, flag_tu, flag_v, flag_tv, flag1, flag) copy(x_s, y_s, x_u, y_u, x_v, y_v) copyin(x1, x) copyin(y1, y) copy(c1, c2, c3, c4, c11, c22, c33, c44, c11_f, c22_f, c33_f, c44_f)
{
        #pragma acc parallel loop collapse(2) private(min_dis, distance,k1,temp,sig,phi)
        for(int i=0;i<n;i++)
        {
            for (int j=0;j<m;j++)
            {
                if(j>c11 && j<c22 && i>c33 && i<c44)
                {
                    min_dis = 10.0;
                    #pragma acc loop
                    for(int k=0;k<N_marker;k++)
                    {
                        distance = sqrt((x_s[k]-x1[i][j])*(x_s[k]-x1[i][j]) + (y_s[k]-y1[i][j])*(y_s[k]-y1[i][j]));
                        if (distance < min_dis)
                        {
                            min_dis = distance;
                            k1 = k;
                        }
                    }
                    temp = (x1[i][j]-x_s[k1])*n1[k1] + (y1[i][j]-y_s[k1])*n2[k1];
                    if (temp>0)
                        sig=1;
                    else
                        sig=-1;
                    phi = sig*min_dis;
                    if (phi>0)
                        flag[i][j]=1;
                    else
                        flag[i][j]=0;
                }
                else
                    flag[i][j]=1;
            }
        }
        //#pragma omp parallel for default(shared) private(min_dis,distance,k1,temp,sig,phi) schedule(dynamic)
        #pragma acc parallel loop collapse(2) private(min_dis,distance,k1,temp,sig,phi)
        for(int i=0;i<n-1;i++)
        {
            for (int j=0;j<m-1;j++)
            {
                if(j>c11 && j<c22 && i>c33 && i<c44)
                {
                    min_dis = 10.0;
                    #pragma acc loop
                    for(int k=0;k<N_marker;k++)
                    {
                        distance = sqrt((x_s[k]-x[i][j])*(x_s[k]-x[i][j]) + (y_s[k]-y[i][j])*(y_s[k]-y[i][j]));
                        if (distance < min_dis)
                        {
                            min_dis = distance;
                            k1 = k;
                        }
                    }
                    temp = (x[i][j]-x_s[k1])*n1[k1] + (y[i][j]-y_s[k1])*n2[k1];
                    if (temp>0)
                        sig=1;
                    else
                        sig=-1;
                    phi = sig*min_dis;
                    if (phi>0)
                        flag1[i][j]=1;
                    else
                        flag1[i][j]=0;
                }
                else
                    flag1[i][j]=1;
            }
        }
        /*for(int i=0;i<n;i++)
         {
         for (int j=0;j<m;j++)
         {
         if(j>c11 && j<c22 && i>c33 && i<c44)
         {
         if (flag[i][j]==1)
         {
         if(flag[i+1][j]==1 && flag[i+1][j+1]==1 && flag[i][j+1]==1 && flag[i-1][j+1]==1 && flag[i-1][j]==1 && flag[i-1][j-1]==1 && flag[i][j-1]==1 && flag[i+1][j-1]==1)
         {
         flag_new[i][j]=1;
         }
         else
         flag_new[i][j]=3;
         }
         }
         else
         flag_new[i][j] = 1;
         }
         }*/
        //#pragma omp parallel for default(shared) private(min_dis,distance,k1,temp,sig,phi) schedule(dynamic)
        #pragma acc parallel loop collapse(2) private(min_dis,distance,k1,temp,sig,phi)
        for(int i=0;i<n_u;i++)
        {
            for (int j=0;j<m_u;j++)
            {
                if(j>c11 && j<c22 && i>c33 && i<c44)
                {
                    min_dis = 10.0;
                    #pragma acc loop
                    for(int k=0;k<N_marker;k++)
                    {
                        distance = sqrt((x_s[k]-x_u[i][j])*(x_s[k]-x_u[i][j]) + (y_s[k]-y_u[i][j])*(y_s[k]-y_u[i][j]));
                        if (distance < min_dis)
                        {
                            min_dis = distance;
                            k1 = k;
                        }
                    }
                    temp = (x_u[i][j]-x_s[k1])*n1[k1] + (y_u[i][j]-y_s[k1])*n2[k1];
                    if (temp>0)
                        sig=1;
                    else
                        sig=-1;
                    phi = sig*min_dis;
                    if (phi>0)
                        flag_tu[i][j]=1;
                    else
                        flag_tu[i][j]=0;
                }
                else
                    flag_tu[i][j]=1;
            }
        }
        #pragma acc parallel loop collapse(2)
        for(int i=0;i<n_u;i++)
        {
            for (int j=0;j<m_u;j++)
            {
                if (flag_tu[i][j]==0)
                {
                    /*if(flag_tu[i+1][j]==0 && flag_tu[i+1][j+1]==0 && flag_tu[i][j+1]==0 && flag_tu[i-1][j+1]==0 && flag_tu[i-1][j]==0 && flag_tu[i-1][j-1]==0 && flag_tu[i][j-1]==0 && flag_tu[i+1][j-1]==0)
                     {
                     flag_u[i][j]=0;
                     }
                     else*/
                    flag_u[i][j]=3;
                }
                else
                    flag_u[i][j]=1;
            }
        }
        //#pragma omp parallel for default(shared) private(min_dis,distance,k1,temp,sig,phi) schedule(dynamic)
        #pragma acc parallel loop collapse(2) private(min_dis,distance,k1,temp,sig,phi)
        for(int i=0;i<n_v;i++)
        {
            for (int j=0;j<m_v;j++)
            {
                if(j>c11 && j<c22 && i>c33 && i<c44)
                {
                    min_dis = 10.0;
                    #pragma acc loop
                    for(int k=0;k<N_marker;k++)
                    {
                        distance = sqrt((x_s[k]-x_v[i][j])*(x_s[k]-x_v[i][j]) + (y_s[k]-y_v[i][j])*(y_s[k]-y_v[i][j]));
                        if (distance < min_dis)
                        {
                            min_dis = distance;
                            k1 = k;
                        }
                    }
                    temp = (x_v[i][j]-x_s[k1])*n1[k1] + (y_v[i][j]-y_s[k1])*n2[k1];
                    if (temp>0)
                        sig=1;
                    else
                        sig=-1;
                    phi = sig*min_dis;
                    if (phi>0)
                        flag_tv[i][j]=1;
                    else
                        flag_tv[i][j]=0;
                }
                else
                    flag_tv[i][j]=1;
            }
        }
        #pragma acc parallel loop collapse(2)
        for(int i=0;i<n_v;i++)
        {
            for (int j=0;j<m_v;j++)
            {
                if (flag_tv[i][j]==0)
                {
                    /*if(flag_tv[i+1][j]==0 && flag_tv[i+1][j+1]==0 && flag_tv[i][j+1]==0 && flag_tv[i-1][j+1]==0 && flag_tv[i-1][j]==0 && flag_tv[i-1][j-1]==0 && flag_tv[i][j-1]==0 && flag_tv[i+1][j-1]==0)
                     {
                     flag_v[i][j]=0;
                     }
                     else*/
                    flag_v[i][j]=3;
                }
                else
                    flag_v[i][j]=1;
            }
        }
        #pragma acc parallel loop collapse(2)
        for(int i=0;i<n_u;i++)
        {
            for(int j=0;j<m_u;j++)
            {
                if (flag_u[i][j]==3)
                    flag_tu[i][j]=1;
                else
                    flag_tu[i][j]=0;
            }
        }
        #pragma acc parallel loop collapse(2)
        for(int i=0;i<n_v;i++)
        {
            for(int j=0;j<m_v;j++)
            {
                if (flag_v[i][j]==3)
                    flag_tv[i][j]=1;
                else
                    flag_tv[i][j]=0;
            }
        }
} // End of GPU loop
        /*if ( (k%writeInterval) == 0)
         {
         string fileName;
         stringstream ss;
         ss << k;
         fileName = ss.str();
         fileName = "grid1_"+fileName+".dat";
         ofstream file8(fileName.c_str());
         if (file8.is_open())
         {
         file8 << "ZONE T=DATA I=" << m << " " << "J=" << n << endl;
         for(int i=0;i<n;i++)
         {
         for(int j=0;j<m;j++)
         {
         file8 << x1[i][j] << " " << y1[i][j] << " " << flag[i][j] << endl;
         }
         }
         file8.close();
         }
         fileName = ss.str();
         fileName = "grid_"+fileName+".dat";
         ofstream file8a(fileName.c_str());
         if (file8a.is_open())
         {
         file8a << "ZONE T=DATA I=" << m-1 << " " << "J=" << n-1 << endl;
         for(int i=0;i<n-1;i++)
         {
         for(int j=0;j<m-1;j++)
         {
         file8a << x[i][j] << " " << y[i][j] << " " << flag1[i][j] << endl;
         }
         }
         file8a.close();
         }
         fileName = ss.str();
         fileName = "grid_u_"+fileName+".dat";
         ofstream file9(fileName.c_str());
         if (file9.is_open())
         {
         file9 << "ZONE T=DATA I=" << m_u << " " << "J=" << n_u << endl;
         for(int i=0;i<n_u;i++)
         {
         for(int j=0;j<m_u;j++)
         {
         file9 << x_u[i][j] << " " << y_u[i][j] << " " << flag_u[i][j] << endl;
         }
         }
         file9.close();
         }
         fileName = ss.str();
         fileName = "grid_v_"+fileName+".dat";
         ofstream file10(fileName.c_str());
         if (file10.is_open())
         {
         file10 << "ZONE T=DATA I=" << m_v << " " << "J=" << n_v << endl;
         for(int i=0;i<n_v;i++)
         {
         for(int j=0;j<m_v;j++)
         {
         file10 << x_v[i][j] << " " << y_v[i][j] << " " << flag_v[i][j] << endl;
         }
         }
         file10.close();
         }
         }
        c_sasv = 0.0;
        for(int i=1;i<n_u-1;i++)
        {
            c_sasv = c_sasv + u[i][m_u-1]*((y[i][m_u-1]-y[i-1][m_u-1])/sum2);
        }*/
        //#pragma omp parallel for default(shared) private(ue_old,uw_old,un_old,us_old,vn_old,vs_old,ue,uw,un,us,vn,vs) schedule(dynamic)
        for(int i=c3; i<c4; i++)
        {
            for (int j=c1; j<c2; j++)
            {
                ue_old = ((x[i][j+1]-x1[i][j+1])/(x[i][j+1]-x[i][j]))*u_old[i][j]+((x1[i][j+1]-x[i][j])/(x[i][j+1]-x[i][j]))*u_old[i][j+1];
                uw_old = ((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u_old[i][j-1]+((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u_old[i][j];
                un_old = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i][j];
                us_old = ((y[i-1][j]-y1[i-1][j])/(y1[i][j]-y1[i-1][j]))*u_old[i][j]+((y1[i][j]-y[i-1][j])/(y1[i][j]-y1[i-1][j]))*u_old[i-1][j];
                vn_old = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j+1];
                vs_old = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i-1][j+1];
                ue = ((x[i][j+1]-x1[i][j+1])/(x[i][j+1]-x[i][j]))*u[i][j]+((x1[i][j+1]-x[i][j])/(x[i][j+1]-x[i][j]))*u[i][j+1];
                uw = ((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u[i][j-1]+((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u[i][j];
                un = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u[i][j];
                us = ((y[i-1][j]-y1[i-1][j])/(y1[i][j]-y1[i-1][j]))*u[i][j]+((y1[i][j]-y[i-1][j])/(y1[i][j]-y1[i-1][j]))*u[i-1][j];
                vn = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v[i][j+1];
                vs = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v[i-1][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v[i-1][j+1];
                u_1[i][j] = u[i][j] - 1.5*((dt/(x1[i][j+1]-x1[i][j]))*(ue*ue-uw*uw)+(dt/(y[i][j]-y[i-1][j]))*(vn*un-vs*us)) + 0.5*((dt/(x1[i][j+1]-x1[i][j]))*(ue_old*ue_old-uw_old*uw_old)+(dt/(y[i][j]-y[i-1][j]))*(vn_old*un_old-vs_old*us_old)) + (dt/(Re*(x1[i][j+1]-x1[i][j])))*(((u[i][j+1]-u[i][j])/(x[i][j+1]-x[i][j]))-((u[i][j]-u[i][j-1])/(x[i][j]-x[i][j-1])))+(dt/(Re*(y[i][j]-y[i-1][j])))*(((u[i+1][j]-u[i][j])/(y1[i+1][j]-y1[i][j]))-((u[i][j]-u[i-1][j])/(y1[i][j]-y1[i-1][j]))) + (dt/(x1[i][j+1]-x1[i][j]))*(pressure[i][j]-pressure[i][j+1]);
            }
        }
        //#pragma omp parallel for default(shared) private(ue_old,uw_old,ve_old,vw_old,vn_old,vs_old,ue,uw,ve,vw,vn,vs) schedule(dynamic)
        for(int i=c3; i<c4; i++)
        {
            for (int j=c1; j<c2; j++)
            {
                ue_old = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u_old[i][j];
                uw_old = ((y[i][j-1]-y1[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i+1][j-1]+((y1[i+1][j-1]-y[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u_old[i][j-1];
                ve_old = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v_old[i][j+1];
                vw_old = ((x1[i][j]-x[i][j-1])/(x1[i][j]-x1[i][j-1]))*v_old[i][j-1]+((x[i][j-1]-x1[i][j-1])/(x1[i][j]-x1[i][j-1]))*v_old[i][j];
                vn_old = ((y1[i+1][j]-y[i][j])/(y[i+1][j]-y[i][j]))*v_old[i+1][j]+((y[i+1][j]-y1[i+1][j])/(y[i+1][j]-y[i][j]))*v_old[i][j];
                vs_old = ((y1[i][j]-y[i-1][j])/(y[i][j]-y[i-1][j]))*v_old[i][j]+((y[i][j]-y1[i][j])/(y[i][j]-y[i-1][j]))*v_old[i-1][j];
                ue = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*u[i+1][j]+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*u[i][j];
                uw = ((y[i][j-1]-y1[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u[i+1][j-1]+((y1[i+1][j-1]-y[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]))*u[i][j-1];
                ve = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*v[i][j]+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*v[i][j+1];
                vw = ((x1[i][j]-x[i][j-1])/(x1[i][j]-x1[i][j-1]))*v[i][j-1]+((x[i][j-1]-x1[i][j-1])/(x1[i][j]-x1[i][j-1]))*v[i][j];
                vn = ((y1[i+1][j]-y[i][j])/(y[i+1][j]-y[i][j]))*v[i+1][j]+((y[i+1][j]-y1[i+1][j])/(y[i+1][j]-y[i][j]))*v[i][j];
                vs = ((y1[i][j]-y[i-1][j])/(y[i][j]-y[i-1][j]))*v[i][j]+((y[i][j]-y1[i][j])/(y[i][j]-y[i-1][j]))*v[i-1][j];
                v_1[i][j] = v[i][j] - 1.5*((dt/(x[i][j]-x[i][j-1]))*(ue*ve-uw*vw)+(dt/(y1[i+1][j]-y1[i][j]))*(vn*vn-vs*vs)) + 0.5*((dt/(x[i][j]-x[i][j-1]))*(ue_old*ve_old-uw_old*vw_old)+(dt/(y1[i+1][j]-y1[i][j]))*(vn_old*vn_old-vs_old*vs_old)) + (dt/(Re*(x[i][j]-x[i][j-1])))*(((v[i][j+1]-v[i][j])/(x1[i][j+1]-x1[i][j]))-((v[i][j]-v[i][j-1])/(x1[i][j]-x1[i][j-1])))+(dt/(Re*(y1[i+1][j]-y1[i][j])))*(((v[i+1][j]-v[i][j])/(y[i+1][j]-y[i][j]))-((v[i][j]-v[i-1][j])/(y[i][j]-y[i-1][j]))) + (dt/(y1[i+1][j]-y1[i][j]))*(pressure[i][j]-pressure[i+1][j]);
            }
        }
        // interpolating the body force at the u forcing point
        for(int i=1; i<n_u-1; i++)
        {
            for (int j=1; j<m_u-1; j++)
            {
                if (flag_u[i][j]!=3)
                {
                    Bf_u[i][j] = 0.0;
                }
                else
                {
                    min_dis = 10.0;
                    for(int ctr=0;ctr<N_marker;ctr++)
                    {
                        distance = sqrt((x_s[ctr]-x_u[i][j])*(x_s[ctr]-x_u[i][j]) + (y_s[ctr]-y_u[i][j])*(y_s[ctr]-y_u[i][j]));
                        if (distance < min_dis)
                        {
                            min_dis = distance;
                            k1 = ctr;
                        }
                    }
                    if(min_dis<=0.000001)
                    {
                        if(k1>=0 && k1<(Nf1*nn1)-(nn1-1))
                        {
                                u_body = u_c_rigid ;
                        }
                        else if(k1>=((Nf2*2)*nn2+(Nf1*nn1)-(nn1-1)) && k1<N_marker)
                        {
                                u_body = u_c_rigid ;
                        }
                        else
                        {
                                u_body = u_c[k1] ;
                        }
                        //u_body = u_c[k1] ;
                        U_intp = u_body;
                        Bf_u[i][j] = (U_intp-u_1[i][j])/dt;
                    }
                    else
                    {
                        if (k1==0)
                        {
                            if((x_s[k1+1]-x_s[N_marker-1])==0)
                            {
                                x_intersect = x_s[N_marker-1];
                                y_intersect = y_u[i][j];
                            }
                            else if ((y_s[k1+1]-y_s[N_marker-1])==0)
                            {
                                x_intersect = x_u[i][j];
                                y_intersect = y_s[N_marker-1];
                            }
                            else
                            {
                                slope = (y_s[k1+1]-y_s[N_marker-1])/(x_s[k1+1]-x_s[N_marker-1]);
                                x_intersect = (y_u[i][j]-y_s[N_marker-1]+slope*x_s[N_marker-1]+(1/slope)*x_u[i][j])/(slope+(1/slope));
                                y_intersect = slope*x_intersect + y_s[N_marker-1] - slope*x_s[N_marker-1];
                            }
                        }
                        else if (k1==N_marker-1)
                        {
                            if((x_s[0]-x_s[k1-1])==0)
                            {
                                x_intersect = x_s[k1-1];
                                y_intersect = y_u[i][j];
                            }
                            else if ((y_s[0]-y_s[k1-1])==0)
                            {
                                x_intersect = x_u[i][j];
                                y_intersect = y_s[k1-1];
                            }
                            else
                            {
                                slope = (y_s[0]-y_s[k1-1])/(x_s[0]-x_s[k1-1]);
                                x_intersect = (y_u[i][j]-y_s[k1-1]+slope*x_s[k1-1]+(1/slope)*x_u[i][j])/(slope+(1/slope));
                                y_intersect = slope*x_intersect + y_s[k1-1] - slope*x_s[k1-1];
                            }
                        }
                        else
                        {
                            if((x_s[k1+1]-x_s[k1-1])==0)
                            {
                                x_intersect = x_s[k1-1];
                                y_intersect = y_u[i][j];
                            }
                            else if ((y_s[k1+1]-y_s[k1-1])==0)
                            {
                                x_intersect = x_u[i][j];
                                y_intersect = y_s[k1-1];
                            }
                            else
                            {
                                slope = (y_s[k1+1]-y_s[k1-1])/(x_s[k1+1]-x_s[k1-1]);
                                x_intersect = (y_u[i][j]-y_s[k1-1]+slope*x_s[k1-1]+(1.0/slope)*x_u[i][j])/(slope+(1.0/slope));
                                y_intersect = slope*x_intersect + y_s[k1-1] - slope*x_s[k1-1];
                            }
                        }
                        //u_body = u_c[k1] ;
                        if(k1>=0 && k1<(Nf1*nn1)-(nn1-1))
                        {
                                u_body = u_c_rigid ;
                        }
                        else if(k1>=((Nf2*2)*nn2+(Nf1*nn1)-(nn1-1))-1 && k1<N_marker)
                        {
                                u_body = u_c_rigid ;
                        }
                        else
                        {
                                theta_intp=atan((y_s[k1+1]-y_s[k1-1])/(x_s[k1+1]-x_s[k1-1])) ;
				if (theta_intp>api/2.0 || theta_intp<-api/2.0)
				{
					r_f0=0.0;
					r_f1=(y_intersect-y_s[k1-1])/sin(theta_intp) ;
					r_f2=(y_s[k1+1]-y_s[k1-1])/sin(theta_intp) ;
				}
				else
				{
					r_f0=0.0;
					r_f1=(x_intersect-x_s[k1-1])/cos(theta_intp) ;
					r_f2=(x_s[k1+1]-x_s[k1-1])/cos(theta_intp) ;
				}
				u_body=((r_f1-r_f0)/(r_f2-r_f0))*u_c[k1+1] + ((r_f2-r_f1)/(r_f2-r_f0))*u_c[k1-1] ;
                        }
                        
                                                /*if (k1==0)
						{
							theta_intp=atan((y_s[k1+1]-y_s[N_marker-1])/(x_s[k1+1]-x_s[N_marker-1])) ;
							if (theta_intp>api/2.0 || theta_intp<-api/2.0)
							{
								r_f0=0.0;
								r_f1=(y_intersect-y_s[N_marker-1])/sin(theta_intp) ;
								r_f2=(y_s[k1+1]-y_s[N_marker-1])/sin(theta_intp) ;
							}
							else
							{
								r_f0=0.0;
								r_f1=(x_intersect-x_s[N_marker-1])/cos(theta_intp) ;
								r_f2=(x_s[k1+1]-x_s[N_marker-1])/cos(theta_intp) ;
							}
							u_body=((r_f1-r_f0)/(r_f2-r_f0))*u_c[k1+1] + ((r_f2-r_f1)/(r_f2-r_f0))*u_c[N_marker-1] ;
						}
						else if (k1==N_marker-1)
						{
							theta_intp=atan((y_s[0]-y_s[k1-1])/(x_s[0]-x_s[k1-1])) ;
							if (theta_intp>api/2.0 || theta_intp<-api/2.0)
							{
								r_f0=0.0;
								r_f1=(y_intersect-y_s[k1-1])/sin(theta_intp) ;
								r_f2=(y_s[0]-y_s[k1-1])/sin(theta_intp) ;
							}
							else
							{
								r_f0=0.0;
								r_f1=(x_intersect-x_s[k1-1])/cos(theta_intp) ;
								r_f2=(x_s[0]-x_s[k1-1])/cos(theta_intp) ;
							}
							u_body=((r_f1-r_f0)/(r_f2-r_f0))*u_c[0] + ((r_f2-r_f1)/(r_f2-r_f0))*u_c[k1-1] ;
						}
						else
						{
							theta_intp=atan((y_s[k1+1]-y_s[k1-1])/(x_s[k1+1]-x_s[k1-1])) ;
							if (theta_intp>api/2.0 || theta_intp<-api/2.0)
							{
								r_f0=0.0;
								r_f1=(y_intersect-y_s[k1-1])/sin(theta_intp) ;
								r_f2=(y_s[k1+1]-y_s[k1-1])/sin(theta_intp) ;
							}
							else
							{
								r_f0=0.0;
								r_f1=(x_intersect-x_s[k1-1])/cos(theta_intp) ;
								r_f2=(x_s[k1+1]-x_s[k1-1])/cos(theta_intp) ;
							}
							u_body=((r_f1-r_f0)/(r_f2-r_f0))*u_c[k1+1] + ((r_f2-r_f1)/(r_f2-r_f0))*u_c[k1-1] ;
						}*/
                        x_mirror = x_intersect + (x_intersect-x_u[i][j]);
                        y_mirror = y_intersect + (y_intersect-y_u[i][j]);
                        n1_i = c3+floor((y_mirror-y_u[c3][c1])/dy2);
                        n1_j = c1+floor((x_mirror-x_u[c3][c1])/dx2);
                        if(flag_u[n1_i][n1_j]==1 && flag_u[n1_i][n1_j+1]==1 && flag_u[n1_i+1][n1_j+1]==1 && flag_u[n1_i+1][n1_j]==1)
                        {
                            u_mirror = (((((x_u[n1_i+1][n1_j+1]-x_mirror)*u_1[n1_i+1][n1_j]+(x_mirror-x_u[n1_i+1][n1_j])*u_1[n1_i+1][n1_j+1])/(x_u[n1_i+1][n1_j+1]-x_u[n1_i+1][n1_j]))*(y_mirror-y_u[n1_i][n1_j]))+((((x_u[n1_i][n1_j+1]-x_mirror)*u_1[n1_i][n1_j]+(x_mirror-x_u[n1_i][n1_j])*u_1[n1_i][n1_j+1])/(x_u[n1_i][n1_j+1]-x_u[n1_i][n1_j]))*(y_u[n1_i+1][n1_j]-y_mirror)))/(y_u[n1_i+1][n1_j]-y_u[n1_i][n1_j]);
                            dist_a = sqrt(pow((y_mirror-y_intersect),2)+pow((x_mirror-x_intersect),2));
                            dist_b = sqrt(pow((y_intersect-y_u[i][j]),2)+pow((x_intersect-x_u[i][j]),2));
                            U_intp = (1+(dist_b/dist_a))*u_body - (dist_b/dist_a)*u_mirror;
                            Bf_u[i][j] = (U_intp-u_1[i][j])/dt;
                        }
                        else
                        {
                            if (k1==0)
                            {
                                if((x_s[k1+1]-x_s[N_marker-1])==0)
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist;
                                        y_mirror1 = y_intersect;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist;
                                        y_mirror1 = y_intersect;
                                    }
                                }
                                else if ((y_s[k1+1]-y_s[N_marker-1])==0)
                                {
                                    if(y_mirror>y_intersect)
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect + dist;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect - dist;
                                    }
                                }
                                else
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect + dist*sin(atan((-1/slope)));
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect - dist*sin(atan((-1/slope)));
                                    }
                                }
                            }
                            else if (k1==N_marker-1)
                            {
                                if((x_s[0]-x_s[k1-1])==0)
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist;
                                        y_mirror1 = y_intersect;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist;
                                        y_mirror1 = y_intersect;
                                    }
                                }
                                else if ((y_s[0]-y_s[k1-1])==0)
                                {
                                    if(y_mirror>y_intersect)
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect + dist;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect - dist;
                                    }
                                }
                                else
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect + dist*sin(atan((-1/slope)));
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect - dist*sin(atan((-1/slope)));
                                    }
                                }
                            }
                            else
                            {
                                if((x_s[k1+1]-x_s[k1-1])==0)
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist;
                                        y_mirror1 = y_intersect;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist;
                                        y_mirror1 = y_intersect;
                                    }
                                }
                                else if ((y_s[k1+1]-y_s[k1-1])==0)
                                {
                                    if(y_mirror>y_intersect)
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect + dist;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect - dist;
                                    }
                                }
                                else
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect + dist*sin(atan((-1/slope)));
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect - dist*sin(atan((-1/slope)));
                                    }
                                }
                            }
                            n1_i = c3+floor((y_mirror1-y_u[c3][c1])/dy2);
                            n1_j = c1+floor((x_mirror1-x_u[c3][c1])/dx2);
                            if(flag_u[n1_i][n1_j]==1 && flag_u[n1_i][n1_j+1]==1 && flag_u[n1_i+1][n1_j+1]==1 && flag_u[n1_i+1][n1_j]==1)
                            {
                                u_mirror = (((((x_u[n1_i+1][n1_j+1]-x_mirror1)*u_1[n1_i+1][n1_j]+(x_mirror1-x_u[n1_i+1][n1_j])*u_1[n1_i+1][n1_j+1])/(x_u[n1_i+1][n1_j+1]-x_u[n1_i+1][n1_j]))*(y_mirror1-y_u[n1_i][n1_j]))+((((x_u[n1_i][n1_j+1]-x_mirror1)*u_1[n1_i][n1_j]+(x_mirror1-x_u[n1_i][n1_j])*u_1[n1_i][n1_j+1])/(x_u[n1_i][n1_j+1]-x_u[n1_i][n1_j]))*(y_u[n1_i+1][n1_j]-y_mirror1)))/(y_u[n1_i+1][n1_j]-y_u[n1_i][n1_j]);
                                dist_a = sqrt(pow((y_mirror1-y_intersect),2)+pow((x_mirror1-x_intersect),2));
                                dist_b = sqrt(pow((y_intersect-y_u[i][j]),2)+pow((x_intersect-x_u[i][j]),2));
                                U_intp = (1+(dist_b/dist_a))*u_body - (dist_b/dist_a)*u_mirror;
                                Bf_u[i][j] = (U_intp-u_1[i][j])/dt;
                            }
                            else
                            {
                                cout << "what will happen to me" <<endl; // --> u -->" << i << " " << j << " " << x_u[i][j] << " " << y_u[i][j] << " " << x_mirror << " " << y_mirror << " " << x_mirror1 << " " << y_mirror1 << " " << x_intersect << " " << y_intersect<< " " << x_s[k1] << " " << y_s[k1] << endl;
                                /*cout << i << " --> " << j  << " --> " << x_u[i][j] << " --> " << y_u[i][j] << " --> " << "what will happen to me" << endl;
                                 cout << x_intersect << " --> " << y_intersect  << " --> " << x_mirror << " --> " << y_mirror << " --> " << x_mirror1 << " -- >" << y_mirror1 << endl;
                                 cout << n1_i << " --> " << n1_j << endl;
                                 cout << flag_u[n1_i][n1_j] << " --> " << flag_u[n1_i][n1_j+1]  << " --> " << flag_u[n1_i+1][n1_j+1] << " --> " << flag_u[n1_i+1][n1_j] << endl;*/
                                flag_u[i][j] = 1;
                                flag_tu[i][j] = 0;
                                Bf_u[i][j] = 0.0;
                            }
                        }
                    }
                }
            }
        }

        // interpolating the body force at the v forcing point
        for(int i=1; i<n_v-1; i++)
        {
            for (int j=1; j<m_v-1; j++)
            {
                if (flag_v[i][j]!=3)
                {
                    Bf_v[i][j] = 0.0;
                }
                else
                {
                    min_dis = 10.0;
                    for(int ctr=0;ctr<N_marker;ctr++)
                    {
                        distance = sqrt((x_s[ctr]-x_v[i][j])*(x_s[ctr]-x_v[i][j]) + (y_s[ctr]-y_v[i][j])*(y_s[ctr]-y_v[i][j]));
                        if (distance < min_dis)
                        {
                            min_dis = distance;
                            k1 = ctr;
                        }
                    }

                    if(min_dis<=0.000001)
                    {
                        if(k1>=0 && k1<(Nf1*nn1)-(nn1-1))
                        {
                                v_body = v_c_rigid ;
                        }
                        else if(k1>=((Nf2*2)*nn2+(Nf1*nn1)-(nn1-1))-1 && k1<N_marker)
                        {
                                v_body = v_c_rigid ;
                        }
                        else
                        {
                              v_body=v_c[k1];
                        }
                        //v_body = v_c[k1] ;
                        V_intp = v_body;
                        Bf_v[i][j] = (V_intp-v_1[i][j])/dt;
                    }
                    else
                    {
                        if (k1==0)
                        {
                            if((x_s[k1+1]-x_s[N_marker-1])==0)
                            {
                                x_intersect = x_s[N_marker-1];
                                y_intersect = y_v[i][j];
                            }
                            else if ((y_s[k1+1]-y_s[N_marker-1])==0)
                            {
                                x_intersect = x_v[i][j];
                                y_intersect = y_s[N_marker-1];
                            }
                            else
                            {
                                slope = (y_s[k1+1]-y_s[N_marker-1])/(x_s[k1+1]-x_s[N_marker-1]);
                                x_intersect = (y_v[i][j]-y_s[N_marker-1]+slope*x_s[N_marker-1]+(1/slope)*x_v[i][j])/(slope+(1/slope));
                                y_intersect = slope*x_intersect + y_s[N_marker-1] - slope*x_s[N_marker-1];
                            }
                        }
                        else if (k1==N_marker-1)
                        {
                            if((x_s[0]-x_s[k1-1])==0)
                            {
                                x_intersect = x_s[k1-1];
                                y_intersect = y_v[i][j];
                            }
                            else if ((y_s[0]-y_s[k1-1])==0)
                            {
                                x_intersect = x_v[i][j];
                                y_intersect = y_s[k1-1];
                            }
                            else
                            {
                                slope = (y_s[0]-y_s[k1-1])/(x_s[0]-x_s[k1-1]);
                                x_intersect = (y_v[i][j]-y_s[k1-1]+slope*x_s[k1-1]+(1/slope)*x_v[i][j])/(slope+(1/slope));
                                y_intersect = slope*x_intersect + y_s[k1-1] - slope*x_s[k1-1];
                            }
                        }
                        else
                        {
                            if((x_s[k1+1]-x_s[k1-1])==0)
                            {
                                x_intersect = x_s[k1-1];
                                y_intersect = y_v[i][j];
                            }
                            else if ((y_s[k1+1]-y_s[k1-1])==0)
                            {
                                x_intersect = x_v[i][j];
                                y_intersect = y_s[k1-1];
                            }
                            else
                            {
                                slope = (y_s[k1+1]-y_s[k1-1])/(x_s[k1+1]-x_s[k1-1]);
                                x_intersect = (y_v[i][j]-y_s[k1-1]+slope*x_s[k1-1]+(1/slope)*x_v[i][j])/(slope+(1/slope));
                                y_intersect = slope*x_intersect + y_s[k1-1] - slope*x_s[k1-1];
                            }
                        }
		        //v_body = v_c[k1];
                        if(k1>=0 && k1<(Nf1*nn1)-(nn1-1))
                        {
                                v_body = v_c_rigid ;
                        }
                        else if(k1>=((Nf2*2)*nn2+(Nf1*nn1)-(nn1-1)) && k1<N_marker)
                        {
                                v_body = v_c_rigid ;
                        }
                        else
                        {
                                theta_intp=atan((y_s[k1+1]-y_s[k1-1])/(x_s[k1+1]-x_s[k1-1])) ;
				if (theta_intp>api/2.0 || theta_intp<-api/2.0)
				{
					r_f0=0.0;
					r_f1=(y_intersect-y_s[k1-1])/sin(theta_intp) ;
					r_f2=(y_s[k1+1]-y_s[k1-1])/sin(theta_intp) ;
				}
				else
				{
					r_f0=0.0;
					r_f1=(x_intersect-x_s[k1-1])/cos(theta_intp) ;
					r_f2=(x_s[k1+1]-x_s[k1-1])/cos(theta_intp) ;
				}
				v_body=((r_f1-r_f0)/(r_f2-r_f0))*v_c[k1+1] + ((r_f2-r_f1)/(r_f2-r_f0))*v_c[k1-1] ;
                        }
                                                /*if (k1==0)
						{
							theta_intp=atan((y_s[k1+1]-y_s[N_marker-1])/(x_s[k1+1]-x_s[N_marker-1])) ;
							if (theta_intp>api/2.0 || theta_intp<-pi/2.0)
							{
								r_f0=0.0;
								r_f1=(y_intersect-y_s[N_marker-1])/sin(theta_intp) ;
								r_f2=(y_s[k1+1]-y_s[N_marker-1])/sin(theta_intp) ;
							}
							else
							{
								r_f0=0.0;
								r_f1=(x_intersect-x_s[N_marker-1])/cos(theta_intp) ;
								r_f2=(x_s[k1+1]-x_s[N_marker-1])/cos(theta_intp) ;
							}
							v_body=((r_f1-r_f0)/(r_f2-r_f0))*v_c[k1+1] + ((r_f2-r_f1)/(r_f2-r_f0))*v_c[N_marker-1] ;
						}
						else if (k1==N_marker-1)
						{
							theta_intp=atan((y_s[0]-y_s[k1-1])/(x_s[0]-x_s[k1-1])) ;
							if (theta_intp>api/2.0 || theta_intp<-api/2.0)
							{
								r_f0=0.0;
								r_f1=(y_intersect-y_s[k1-1])/sin(theta_intp) ;
								r_f2=(y_s[0]-y_s[k1-1])/sin(theta_intp) ;
							}
							else
							{
								r_f0=0.0;
								r_f1=(x_intersect-x_s[k1-1])/cos(theta_intp) ;
								r_f2=(x_s[0]-x_s[k1-1])/cos(theta_intp) ;
							}
							 v_body=((r_f1-r_f0)/(r_f2-r_f0))*v_c[0] + ((r_f2-r_f1)/(r_f2-r_f0))*v_c[k1-1] ;
						}
						else
						{
							theta_intp=atan((y_s[k1+1]-y_s[k1-1])/(x_s[k1+1]-x_s[k1-1])) ;
							if (theta_intp>api/2.0 || theta_intp<-api/2.0)
							{
								r_f0=0.0;
								r_f1=(y_intersect-y_s[k1-1])/sin(theta_intp) ;
								r_f2=(y_s[k1+1]-y_s[k1-1])/sin(theta_intp) ;
							}
							else
							{
								r_f0=0.0;
								r_f1=(x_intersect-x_s[k1-1])/cos(theta_intp) ;
								r_f2=(x_s[k1+1]-x_s[k1-1])/cos(theta_intp) ;
							}
							v_body=((r_f1-r_f0)/(r_f2-r_f0))*v_c[k1+1] + ((r_f2-r_f1)/(r_f2-r_f0))*v_c[k1-1] ;
						}*/
                        x_mirror = x_intersect + (x_intersect-x_v[i][j]);
                        y_mirror = y_intersect + (y_intersect-y_v[i][j]);
                        n1_i = c3+floor((y_mirror-y_v[c3][c1])/dy2);
                        n1_j = c1+floor((x_mirror-x_v[c3][c1])/dx2);
                        if(flag_v[n1_i][n1_j]==1 && flag_v[n1_i][n1_j+1]==1 && flag_v[n1_i+1][n1_j+1]==1 && flag_v[n1_i+1][n1_j]==1)
                        {
                            v_mirror = (((((x_v[n1_i+1][n1_j+1]-x_mirror)*v_1[n1_i+1][n1_j]+(x_mirror-x_v[n1_i+1][n1_j])*v_1[n1_i+1][n1_j+1])/(x_v[n1_i+1][n1_j+1]-x_v[n1_i+1][n1_j]))*(y_mirror-y_v[n1_i][n1_j]))+((((x_v[n1_i][n1_j+1]-x_mirror)*v_1[n1_i][n1_j]+(x_mirror-x_v[n1_i][n1_j])*v_1[n1_i][n1_j+1])/(x_v[n1_i][n1_j+1]-x_v[n1_i][n1_j]))*(y_v[n1_i+1][n1_j]-y_mirror)))/(y_v[n1_i+1][n1_j]-y_v[n1_i][n1_j]);
                            dist_a = sqrt(pow((y_mirror-y_intersect),2)+pow((x_mirror-x_intersect),2));
                            dist_b = sqrt(pow((y_intersect-y_v[i][j]),2)+pow((x_intersect-x_v[i][j]),2));
                            V_intp = (1+(dist_b/dist_a))*v_body - (dist_b/dist_a)*v_mirror;
                            Bf_v[i][j] = (V_intp-v_1[i][j])/dt;
                        }
                        else
                        {
                            if (k1==0)
                            {
                                if((x_s[k1+1]-x_s[N_marker-1])==0)
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist;
                                        y_mirror1 = y_intersect;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist;
                                        y_mirror1 = y_intersect;
                                    }
                                }
                                else if ((y_s[k1+1]-y_s[N_marker-1])==0)
                                {
                                    if(y_mirror>y_intersect)
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect + dist;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect - dist;
                                    }
                                }
                                else
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect + dist*sin(atan((-1/slope)));
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect - dist*sin(atan((-1/slope)));
                                    }
                                }
                            }
                            else if (k1==N_marker-1)
                            {
                                if((x_s[0]-x_s[k1-1])==0)
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist;
                                        y_mirror1 = y_intersect;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist;
                                        y_mirror1 = y_intersect;
                                    }
                                }
                                else if ((y_s[0]-y_s[k1-1])==0)
                                {
                                    if(y_mirror>y_intersect)
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect + dist;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect - dist;
                                    }
                                }
                                else
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect + dist*sin(atan((-1/slope)));
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect - dist*sin(atan((-1/slope)));
                                    }
                                }
                            }
                            else
                            {
                                if((x_s[k1+1]-x_s[k1-1])==0)
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist;
                                        y_mirror1 = y_intersect;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist;
                                        y_mirror1 = y_intersect;
                                    }
                                }
                                else if ((y_s[k1+1]-y_s[k1-1])==0)
                                {
                                    if(y_mirror>y_intersect)
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect + dist;
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect;
                                        y_mirror1 = y_intersect - dist;
                                    }
                                }
                                else
                                {
                                    if(x_mirror>x_intersect)
                                    {
                                        x_mirror1 = x_intersect + dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect + dist*sin(atan((-1/slope)));
                                    }
                                    else
                                    {
                                        x_mirror1 = x_intersect - dist*cos(atan((-1/slope)));
                                        y_mirror1 = y_intersect - dist*sin(atan((-1/slope)));
                                    }
                                }
                            }
                            n1_i = c3+floor((y_mirror1-y_v[c3][c1])/dy2);
                            n1_j = c1+floor((x_mirror1-x_v[c3][c1])/dx2);
                            if(flag_v[n1_i][n1_j]==1 && flag_v[n1_i][n1_j+1]==1 && flag_v[n1_i+1][n1_j+1]==1 && flag_v[n1_i+1][n1_j]==1)
                            {
                                v_mirror = (((((x_v[n1_i+1][n1_j+1]-x_mirror1)*v_1[n1_i+1][n1_j]+(x_mirror1-x_v[n1_i+1][n1_j])*v_1[n1_i+1][n1_j+1])/(x_v[n1_i+1][n1_j+1]-x_v[n1_i+1][n1_j]))*(y_mirror1-y_v[n1_i][n1_j]))+((((x_v[n1_i][n1_j+1]-x_mirror1)*v_1[n1_i][n1_j]+(x_mirror1-x_v[n1_i][n1_j])*v_1[n1_i][n1_j+1])/(x_v[n1_i][n1_j+1]-x_v[n1_i][n1_j]))*(y_v[n1_i+1][n1_j]-y_mirror1)))/(y_v[n1_i+1][n1_j]-y_v[n1_i][n1_j]);
                                dist_a = sqrt(pow((y_mirror1-y_intersect),2)+pow((x_mirror1-x_intersect),2));
                                dist_b = sqrt(pow((y_intersect-y_v[i][j]),2)+pow((x_intersect-x_v[i][j]),2));
                                V_intp = (1+(dist_b/dist_a))*v_body - (dist_b/dist_a)*v_mirror;
                                Bf_v[i][j] = (V_intp-v_1[i][j])/dt;
                            }
                            else
                            {
                                cout << "what will happen to me" <<endl; // --> v -->" << i << " " << j << " " << x_v[i][j] << " " << y_v[i][j] << " " << x_mirror << " " << y_mirror << " " << x_mirror1 << " " << y_mirror1 << " " << x_intersect << " " << y_intersect<< " " << x_s[k1] << " " << y_s[k1] << endl;
                                flag_v[i][j] = 1;
                                flag_tv[i][j] = 0;
                                v[i][j] = v_1[i][j];
                                Bf_v[i][j] = 0.0;
                            }
                        }
                    }
                }
            }
        }
        // Calculation of the f_u of x-momentum equation
        //#pragma omp parallel for default(shared) private(ue_old,uw_old,un_old,us_old,vn_old,vs_old,ue,uw,un,us,vn,vs,A1,B1,C1) schedule(dynamic)
	for(int i=1; i<n_u-1; i++)
        {
            for (int j=1; j<m_u-1; j++)
            {
                ue_old = ((x[i][j+1]-x1[i][j+1])*u_old[i][j]+(x1[i][j+1]-x[i][j])*u_old[i][j+1])/(x[i][j+1]-x[i][j]);
                uw_old = ((x[i][j]-x1[i][j])*u_old[i][j-1]+(x1[i][j]-x[i][j-1])*u_old[i][j])/(x[i][j]-x[i][j-1]);
                un_old = ((y[i][j]-y1[i][j])*u_old[i+1][j]+(y1[i+1][j]-y[i][j])*u_old[i][j])/(y1[i+1][j]-y1[i][j]);
                us_old = ((y[i-1][j]-y1[i-1][j])*u_old[i][j]+(y1[i][j]-y[i-1][j])*u_old[i-1][j])/(y1[i][j]-y1[i-1][j]);
                vn_old = ((x1[i][j+1]-x[i][j])*v_old[i][j]+(x[i][j]-x1[i][j])*v_old[i][j+1])/(x1[i][j+1]-x1[i][j]);
                vs_old = ((x1[i][j+1]-x[i][j])*v_old[i-1][j]+(x[i][j]-x1[i][j])*v_old[i-1][j+1])/(x1[i][j+1]-x1[i][j]);
                ue = ((x[i][j+1]-x1[i][j+1])*u[i][j]+(x1[i][j+1]-x[i][j])*u[i][j+1])/(x[i][j+1]-x[i][j]);
                uw = ((x[i][j]-x1[i][j])*u[i][j-1]+(x1[i][j]-x[i][j-1])*u[i][j])/(x[i][j]-x[i][j-1]);
                un = ((y[i][j]-y1[i][j])*u[i+1][j]+(y1[i+1][j]-y[i][j])*u[i][j])/(y1[i+1][j]-y1[i][j]);
                us = ((y[i-1][j]-y1[i-1][j])*u[i][j]+(y1[i][j]-y[i-1][j])*u[i-1][j])/(y1[i][j]-y1[i-1][j]);
                vn = ((x1[i][j+1]-x[i][j])*v[i][j]+(x[i][j]-x1[i][j])*v[i][j+1])/(x1[i][j+1]-x1[i][j]);
                vs = ((x1[i][j+1]-x[i][j])*v[i-1][j]+(x[i][j]-x1[i][j])*v[i-1][j+1])/(x1[i][j+1]-x1[i][j]);
                A1 = 1.5*((dt/(x1[i][j+1]-x1[i][j]))*(ue*ue-uw*uw)+(dt/(y[i][j]-y[i-1][j]))*(vn*un-vs*us));
                B1 = 0.5*((dt/(x1[i][j+1]-x1[i][j]))*(ue_old*ue_old-uw_old*uw_old)+(dt/(y[i][j]-y[i-1][j]))*(vn_old*un_old-vs_old*us_old));
                C1 = (dt/(2.0*Re*(x1[i][j+1]-x1[i][j])))*(((u[i][j+1]-u[i][j])/(x[i][j+1]-x[i][j]))-((u[i][j]-u[i][j-1])/(x[i][j]-x[i][j-1])))+(dt/(2.0*Re*(y[i][j]-y[i-1][j])))*(((u[i+1][j]-u[i][j])/(y1[i+1][j]-y1[i][j]))-((u[i][j]-u[i-1][j])/(y1[i][j]-y1[i-1][j])));
                f_u[i][j] = u[i][j]-A1+B1+C1;
            }
        }
        // Calculation of the f_v term of the y-momentum equation
        //#pragma omp parallel for default(shared) private(ue_old,uw_old,ve_old,vw_old,vn_old,vs_old,ue,uw,ve,vw,vn,vs,A1,B1,C1) schedule(dynamic)
        for(int i=1; i<n_v-1; i++)
        {
            for (int j=1; j<m_v-1; j++)
            {
                ue_old = ((y[i][j]-y1[i][j])*u_old[i+1][j]+(y1[i+1][j]-y[i][j])*u_old[i][j])/(y1[i+1][j]-y1[i][j]);
                uw_old = ((y[i][j-1]-y1[i][j-1])*u_old[i+1][j-1]+(y1[i+1][j-1]-y[i][j-1])*u_old[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]);
                ve_old = ((x1[i][j+1]-x[i][j])*v_old[i][j]+(x[i][j]-x1[i][j])*v_old[i][j+1])/(x1[i][j+1]-x1[i][j]);
                vw_old = ((x1[i][j]-x[i][j-1])*v_old[i][j-1]+(x[i][j-1]-x1[i][j-1])*v_old[i][j])/(x1[i][j]-x1[i][j-1]);
                vn_old = ((y1[i+1][j]-y[i][j])*v_old[i+1][j]+(y[i+1][j]-y1[i+1][j])*v_old[i][j])/(y[i+1][j]-y[i][j]);
                vs_old = ((y1[i][j]-y[i-1][j])*v_old[i][j]+(y[i][j]-y1[i][j])*v_old[i-1][j])/(y[i][j]-y[i-1][j]);
                ue = ((y[i][j]-y1[i][j])*u[i+1][j]+(y1[i+1][j]-y[i][j])*u[i][j])/(y1[i+1][j]-y1[i][j]);
                uw = ((y[i][j-1]-y1[i][j-1])*u[i+1][j-1]+(y1[i+1][j-1]-y[i][j-1])*u[i][j-1])/(y1[i+1][j-1]-y1[i][j-1]);
                ve = ((x1[i][j+1]-x[i][j])*v[i][j]+(x[i][j]-x1[i][j])*v[i][j+1])/(x1[i][j+1]-x1[i][j]);
                vw = ((x1[i][j]-x[i][j-1])*v[i][j-1]+(x[i][j-1]-x1[i][j-1])*v[i][j])/(x1[i][j]-x1[i][j-1]);
                vn = ((y1[i+1][j]-y[i][j])*v[i+1][j]+(y[i+1][j]-y1[i+1][j])*v[i][j])/(y[i+1][j]-y[i][j]);
                vs = ((y1[i][j]-y[i-1][j])*v[i][j]+(y[i][j]-y1[i][j])*v[i-1][j])/(y[i][j]-y[i-1][j]);
                A1 = 1.5*((dt/(x[i][j]-x[i][j-1]))*(ue*ve-uw*vw)+(dt/(y1[i+1][j]-y1[i][j]))*(vn*vn-vs*vs));
                B1 = 0.5*((dt/(x[i][j]-x[i][j-1]))*(ue_old*ve_old-uw_old*vw_old)+(dt/(y1[i+1][j]-y1[i][j]))*(vn_old*vn_old-vs_old*vs_old));
                C1 = (dt/(2.0*Re*(x[i][j]-x[i][j-1])))*(((v[i][j+1]-v[i][j])/(x1[i][j+1]-x1[i][j]))-((v[i][j]-v[i][j-1])/(x1[i][j]-x1[i][j-1])))+(dt/(2.0*Re*(y1[i+1][j]-y1[i][j])))*(((v[i+1][j]-v[i][j])/(y[i+1][j]-y[i][j]))-((v[i][j]-v[i-1][j])/(y[i][j]-y[i-1][j])));
                f_v[i][j] = v[i][j]-A1+B1+C1;
            }
        }

        for(int k_f=0; k_f<Nf; k_f++)
        {
            temp_Dfc=0.0;
			minfx1=min(x_sfc[k_f],x_sfc[k_f+1]);
			minfx2=min(x_sfc[N_markerf-k_f-1],x_sfc[N_markerf-k_f-2]);
			maxfx1=max(x_sfc[k_f],x_sfc[k_f+1]);
			maxfx2=max(x_sfc[N_markerf-k_f-1],x_sfc[N_markerf-k_f-2]);
			minfy1=min(y_sfc[k_f],y_sfc[k_f+1]);
			minfy2=min(y_sfc[N_markerf-k_f-1],y_sfc[N_markerf-k_f-2]);
			maxfy1=max(y_sfc[k_f],y_sfc[k_f+1]);
			maxfy2=max(y_sfc[N_markerf-k_f-1],y_sfc[N_markerf-k_f-2]);
			cf1=m_div1+ceil((min(minfx1,minfx2)-Hx1)/dx2)-10;
			cf2=m_div1+ceil((max(maxfx1,maxfx2)-Hx1)/dx2)+10;
			cf3=n_div1+ceil((min(minfy1,minfy2)-Hy1)/dy2)-10;
			cf4=n_div1+ceil((max(maxfy1,maxfy2)-Hy1)/dy2)+10;

			if ((y_sfc[k_f]-y_sfc[N_markerf-k_f-1])>=0.0)
			{
				Afc1=(y_sfc[k_f]-y_sfc[N_markerf-k_f-1]);
				Afc2=(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2]);
				Bfc1=-(x_sfc[k_f]-x_sfc[N_markerf-k_f-1]);
				Bfc2=-(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2]);;
				Cfc1=((y_sfc[N_markerf-k_f-1]*(x_sfc[k_f]-x_sfc[N_markerf-k_f-1])) - (x_sfc[N_markerf-k_f-1]*(y_sfc[k_f]-y_sfc[N_markerf-k_f-1])));
				Cfc2=((y_sfc[N_markerf-k_f-2]*(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2])) - (x_sfc[N_markerf-k_f-2]*(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2])));
			}
			else
			{
				Afc1=-(y_sfc[k_f]-y_sfc[N_markerf-k_f-1]);
				Afc2=-(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2]);
				Bfc1=(x_sfc[k_f]-x_sfc[N_markerf-k_f-1]);
				Bfc2=(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2]);;
				Cfc1=-((y_sfc[N_markerf-k_f-1]*(x_sfc[k_f]-x_sfc[N_markerf-k_f-1])) - (x_sfc[N_markerf-k_f-1]*(y_sfc[k_f]-y_sfc[N_markerf-k_f-1])));
				Cfc2=-((y_sfc[N_markerf-k_f-2]*(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2])) - (x_sfc[N_markerf-k_f-2]*(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2])));
			}
			//cout<< cf1 << " " <<cf2 << " " << cf3 << " " << cf4 <<endl;
			for (int i=cf3;i<cf4;i++)
			{
				for (int j=cf1;j<cf2;j++)
				{
					if (flag[i][j]==0)
					{
						if (thetafc[k_f]!=-api && (xf[k_f]-xf[k_f-1])>=0.0 && (Afc1*x1[i][j]+Bfc1*y1[i][j]+Cfc1)>=0.0)
						{
							if (thetafc[k_f+1]!=-api && (xf[k_f+1]-xf[k_f])>=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)<=0.0)
							{
								uf_old = ((x1[i][j]-x[i][j-1])*u_old[i][j]+(x[i][j]-x1[i][j])*u_old[i][j-1])/(x[i][j]-x[i][j-1]);
								uf_n = ((x1[i+1][j]-x[i+1][j-1])*u_old[i+1][j]+(x[i+1][j]-x1[i+1][j])*u_old[i+1][j-1])/(x[i+1][j]-x[i+1][j-1]);
								uf_s = ((x1[i-1][j]-x[i-1][j-1])*u_old[i-1][j]+(x[i-1][j]-x1[i-1][j])*u_old[i-1][j-1])/(x[i-1][j]-x[i-1][j-1]);
								un = ((y[i][j]-y1[i][j])*uf_n+(y1[i+1][j]-y[i][j])*uf_old)/(y1[i+1][j]-y1[i][j]);
								us = ((y[i-1][j]-y1[i-1][j])*uf_old+(y1[i][j]-y[i-1][j])*uf_s)/(y1[i][j]-y1[i-1][j]);
								temp_Dfc = temp_Dfc + ((-0.5)*(((u_old[i][j]*u_old[i][j]-u_old[i][j-1]*u_old[i][j-1])/(x[i][j]-x[i][j-1]))+((un*v_old[i][j]-us*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
							}
							else if (thetafc[k_f+1]!=api && (xf[k_f+1]-xf[k_f])<=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)>=0.0)
							{
								uf_old = ((x1[i][j]-x[i][j-1])*u_old[i][j]+(x[i][j]-x1[i][j])*u_old[i][j-1])/(x[i][j]-x[i][j-1]);
								uf_n = ((x1[i+1][j]-x[i+1][j-1])*u_old[i+1][j]+(x[i+1][j]-x1[i+1][j])*u_old[i+1][j-1])/(x[i+1][j]-x[i+1][j-1]);
								uf_s = ((x1[i-1][j]-x[i-1][j-1])*u_old[i-1][j]+(x[i-1][j]-x1[i-1][j])*u_old[i-1][j-1])/(x[i-1][j]-x[i-1][j-1]);
								un = ((y[i][j]-y1[i][j])*uf_n+(y1[i+1][j]-y[i][j])*uf_old)/(y1[i+1][j]-y1[i][j]);
								us = ((y[i-1][j]-y1[i-1][j])*uf_old+(y1[i][j]-y[i-1][j])*uf_s)/(y1[i][j]-y1[i-1][j]);
								temp_Dfc = temp_Dfc + ((-0.5)*(((u_old[i][j]*u_old[i][j]-u_old[i][j-1]*u_old[i][j-1])/(x[i][j]-x[i][j-1]))+((un*v_old[i][j]-us*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
							}
							/*else
							{
								cout<<"wrong no.1d"<<endl;
							}*/
						}
						else if (thetafc[k_f]!=api && (xf[k_f]-xf[k_f-1])<=0.0 && (Afc1*x1[i][j]+Bfc1*y1[i][j]+Cfc1)<=0.0)
						{
							if (thetafc[k_f+1]!=-api && (xf[k_f+1]-xf[k_f])>=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)<=0.0)
							{
								uf_old = ((x1[i][j]-x[i][j-1])*u_old[i][j]+(x[i][j]-x1[i][j])*u_old[i][j-1])/(x[i][j]-x[i][j-1]);
								uf_n = ((x1[i+1][j]-x[i+1][j-1])*u_old[i+1][j]+(x[i+1][j]-x1[i+1][j])*u_old[i+1][j-1])/(x[i+1][j]-x[i+1][j-1]);
								uf_s = ((x1[i-1][j]-x[i-1][j-1])*u_old[i-1][j]+(x[i-1][j]-x1[i-1][j])*u_old[i-1][j-1])/(x[i-1][j]-x[i-1][j-1]);
								un = ((y[i][j]-y1[i][j])*uf_n+(y1[i+1][j]-y[i][j])*uf_old)/(y1[i+1][j]-y1[i][j]);
								us = ((y[i-1][j]-y1[i-1][j])*uf_old+(y1[i][j]-y[i-1][j])*uf_s)/(y1[i][j]-y1[i-1][j]);
								temp_Dfc = temp_Dfc + ((-0.5)*(((u_old[i][j]*u_old[i][j]-u_old[i][j-1]*u_old[i][j-1])/(x[i][j]-x[i][j-1]))+((un*v_old[i][j]-us*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
							}
							else if (thetafc[k_f+1]!=api && (xf[k_f+1]-xf[k_f])<=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)>=0.0)
							{
								uf_old = ((x1[i][j]-x[i][j-1])*u_old[i][j]+(x[i][j]-x1[i][j])*u_old[i][j-1])/(x[i][j]-x[i][j-1]);
								uf_n = ((x1[i+1][j]-x[i+1][j-1])*u_old[i+1][j]+(x[i+1][j]-x1[i+1][j])*u_old[i+1][j-1])/(x[i+1][j]-x[i+1][j-1]);
								uf_s = ((x1[i-1][j]-x[i-1][j-1])*u_old[i-1][j]+(x[i-1][j]-x1[i-1][j])*u_old[i-1][j-1])/(x[i-1][j]-x[i-1][j-1]);
								un = ((y[i][j]-y1[i][j])*uf_n+(y1[i+1][j]-y[i][j])*uf_old)/(y1[i+1][j]-y1[i][j]);
								us = ((y[i-1][j]-y1[i-1][j])*uf_old+(y1[i][j]-y[i-1][j])*uf_s)/(y1[i][j]-y1[i-1][j]);
								temp_Dfc = temp_Dfc + ((-0.5)*(((u_old[i][j]*u_old[i][j]-u_old[i][j-1]*u_old[i][j-1])/(x[i][j]-x[i][j-1]))+((un*v_old[i][j]-us*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
							}
							/*else
							{
								cout<<"wrong no.2d"<<endl;
							}*/
						}
					}
				}
			}
            temp_Df[k_f]=temp_Dfc;
        }
        temp_D=0.0;
        //#pragma omp parallel for default(shared) private(uf_old,uf_n,uf_s,un,us) reduction(+:temp_D) schedule(dynamic)
        for(int i=c33_f;i<c44_f;i++)
        {
            for (int j=c11_f;j<c22_f;j++)
            {
				if(flag[i][j]==0)
				{
					uf_old = ((x1[i][j]-x[i][j-1])*u_old[i][j]+(x[i][j]-x1[i][j])*u_old[i][j-1])/(x[i][j]-x[i][j-1]);
					uf_n = ((x1[i+1][j]-x[i+1][j-1])*u_old[i+1][j]+(x[i+1][j]-x1[i+1][j])*u_old[i+1][j-1])/(x[i+1][j]-x[i+1][j-1]);
					uf_s = ((x1[i-1][j]-x[i-1][j-1])*u_old[i-1][j]+(x[i-1][j]-x1[i-1][j])*u_old[i-1][j-1])/(x[i-1][j]-x[i-1][j-1]);
					un = ((y[i][j]-y1[i][j])*uf_n+(y1[i+1][j]-y[i][j])*uf_old)/(y1[i+1][j]-y1[i][j]);
					us = ((y[i-1][j]-y1[i-1][j])*uf_old+(y1[i][j]-y[i-1][j])*uf_s)/(y1[i][j]-y1[i-1][j]);
					temp_D += ((-0.5)*(((u_old[i][j]*u_old[i][j]-u_old[i][j-1]*u_old[i][j-1])/(x[i][j]-x[i][j-1]))+((un*v_old[i][j]-us*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
				}
            }
        }

        for(int k_f=0; k_f<Nf; k_f++)
        {
            temp_Lfc=0.0;
			minfx1=min(x_sfc[k_f],x_sfc[k_f+1]);
			minfx2=min(x_sfc[N_markerf-k_f-1],x_sfc[N_markerf-k_f-2]);
			maxfx1=max(x_sfc[k_f],x_sfc[k_f+1]);
			maxfx2=max(x_sfc[N_markerf-k_f-1],x_sfc[N_markerf-k_f-2]);
			minfy1=min(y_sfc[k_f],y_sfc[k_f+1]);
			minfy2=min(y_sfc[N_markerf-k_f-1],y_sfc[N_markerf-k_f-2]);
			maxfy1=max(y_sfc[k_f],y_sfc[k_f+1]);
			maxfy2=max(y_sfc[N_markerf-k_f-1],y_sfc[N_markerf-k_f-2]);
			cf1=m_div1+ceil((min(minfx1,minfx2)-Hx1)/dx2)-10;
			cf2=m_div1+ceil((max(maxfx1,maxfx2)-Hx1)/dx2)+10;
			cf3=n_div1+ceil((min(minfy1,minfy2)-Hy1)/dy2)-10;
			cf4=n_div1+ceil((max(maxfy1,maxfy2)-Hy1)/dy2)+10;

			if ((y_sfc[k_f]-y_sfc[N_markerf-k_f-1])>=0.0)
				{
					Afc1=(y_sfc[k_f]-y_sfc[N_markerf-k_f-1]);
					Afc2=(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2]);
					Bfc1=-(x_sfc[k_f]-x_sfc[N_markerf-k_f-1]);
					Bfc2=-(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2]);;
					Cfc1=((y_sfc[N_markerf-k_f-1]*(x_sfc[k_f]-x_sfc[N_markerf-k_f-1])) - (x_sfc[N_markerf-k_f-1]*(y_sfc[k_f]-y_sfc[N_markerf-k_f-1])));
					Cfc2=((y_sfc[N_markerf-k_f-2]*(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2])) - (x_sfc[N_markerf-k_f-2]*(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2])));
				}
				else
				{
					Afc1=-(y_sfc[k_f]-y_sfc[N_markerf-k_f-1]);
					Afc2=-(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2]);
					Bfc1=(x_sfc[k_f]-x_sfc[N_markerf-k_f-1]);
					Bfc2=(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2]);;
					Cfc1=-((y_sfc[N_markerf-k_f-1]*(x_sfc[k_f]-x_sfc[N_markerf-k_f-1])) - (x_sfc[N_markerf-k_f-1]*(y_sfc[k_f]-y_sfc[N_markerf-k_f-1])));
					Cfc2=-((y_sfc[N_markerf-k_f-2]*(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2])) - (x_sfc[N_markerf-k_f-2]*(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2])));
				}
				//cout<< cf1 << " " <<cf2 << " " << cf3 << " " << cf4 <<endl;
			for (int i=cf3;i<cf4;i++)
			{
				for (int j=cf1;j<cf2;j++)
				{
					if (flag[i][j]==0)
					{
						if (thetafc[k_f]!=-api && (xf[k_f]-xf[k_f-1])>=0.0 && (Afc1*x1[i][j]+Bfc1*y1[i][j]+Cfc1)>=0.0)
						{
							if (thetafc[k_f+1]!=-api && (xf[k_f+1]-xf[k_f])>=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)<=0.0)
							{
								vf_old = ((y1[i][j]-y[i-1][j])*v_old[i][j]+(y[i][j]-y1[i][j])*v_old[i-1][j])/(y[i][j]-y[i-1][j]);
								vf_e = ((y1[i][j+1]-y[i-1][j+1])*v_old[i][j+1]+(y[i][j+1]-y1[i][j+1])*v_old[i-1][j+1])/(y[i][j+1]-y[i-1][j+1]);
								vf_w = ((y1[i][j-1]-y[i-1][j-1])*v_old[i][j-1]+(y[i][j-1]-y1[i][j-1])*v_old[i-1][j-1])/(y[i][j-1]-y[i-1][j-1]);
								ve = ((x1[i][j+1]-x[i][j])*vf_old+(x[i][j]-x1[i][j])*vf_e)/(x1[i][j+1]-x1[i][j]);
								vw = ((x1[i][j]-x[i][j-1])*vf_w+(x[i][j-1]-x1[i][j-1])*vf_old)/(x1[i][j]-x1[i][j-1]);
								temp_Lfc = temp_Lfc + ((-0.5)*(((u_old[i][j]*ve-u_old[i][j-1]*vw)/(x[i][j]-x[i][j-1]))+((v_old[i][j]*v_old[i][j]-v_old[i-1][j]*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
							}
							else if (thetafc[k_f+1]!=api && (xf[k_f+1]-xf[k_f])<=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)>=0.0)
							{
								vf_old = ((y1[i][j]-y[i-1][j])*v_old[i][j]+(y[i][j]-y1[i][j])*v_old[i-1][j])/(y[i][j]-y[i-1][j]);
								vf_e = ((y1[i][j+1]-y[i-1][j+1])*v_old[i][j+1]+(y[i][j+1]-y1[i][j+1])*v_old[i-1][j+1])/(y[i][j+1]-y[i-1][j+1]);
								vf_w = ((y1[i][j-1]-y[i-1][j-1])*v_old[i][j-1]+(y[i][j-1]-y1[i][j-1])*v_old[i-1][j-1])/(y[i][j-1]-y[i-1][j-1]);
								ve = ((x1[i][j+1]-x[i][j])*vf_old+(x[i][j]-x1[i][j])*vf_e)/(x1[i][j+1]-x1[i][j]);
								vw = ((x1[i][j]-x[i][j-1])*vf_w+(x[i][j-1]-x1[i][j-1])*vf_old)/(x1[i][j]-x1[i][j-1]);
								temp_Lfc = temp_Lfc + ((-0.5)*(((u_old[i][j]*ve-u_old[i][j-1]*vw)/(x[i][j]-x[i][j-1]))+((v_old[i][j]*v_old[i][j]-v_old[i-1][j]*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
							}
							/*else
							{
								cout<<"wrong no.1L"<<endl;
							}*/
						}
						else if (thetafc[k_f]!=api && (xf[k_f]-xf[k_f-1])<=0.0 && (Afc1*x1[i][j]+Bfc1*y1[i][j]+Cfc1)<=0.0)
						{
							if (thetafc[k_f+1]!=-api && (xf[k_f+1]-xf[k_f])>=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)<=0.0)
							{
								vf_old = ((y1[i][j]-y[i-1][j])*v_old[i][j]+(y[i][j]-y1[i][j])*v_old[i-1][j])/(y[i][j]-y[i-1][j]);
								vf_e = ((y1[i][j+1]-y[i-1][j+1])*v_old[i][j+1]+(y[i][j+1]-y1[i][j+1])*v_old[i-1][j+1])/(y[i][j+1]-y[i-1][j+1]);
								vf_w = ((y1[i][j-1]-y[i-1][j-1])*v_old[i][j-1]+(y[i][j-1]-y1[i][j-1])*v_old[i-1][j-1])/(y[i][j-1]-y[i-1][j-1]);
								ve = ((x1[i][j+1]-x[i][j])*vf_old+(x[i][j]-x1[i][j])*vf_e)/(x1[i][j+1]-x1[i][j]);
								vw = ((x1[i][j]-x[i][j-1])*vf_w+(x[i][j-1]-x1[i][j-1])*vf_old)/(x1[i][j]-x1[i][j-1]);
								temp_Lfc = temp_Lfc + ((-0.5)*(((u_old[i][j]*ve-u_old[i][j-1]*vw)/(x[i][j]-x[i][j-1]))+((v_old[i][j]*v_old[i][j]-v_old[i-1][j]*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
							}
							else if (thetafc[k_f+1]!=api && (xf[k_f+1]-xf[k_f])<=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)>=0.0)
							{
								vf_old = ((y1[i][j]-y[i-1][j])*v_old[i][j]+(y[i][j]-y1[i][j])*v_old[i-1][j])/(y[i][j]-y[i-1][j]);
								vf_e = ((y1[i][j+1]-y[i-1][j+1])*v_old[i][j+1]+(y[i][j+1]-y1[i][j+1])*v_old[i-1][j+1])/(y[i][j+1]-y[i-1][j+1]);
								vf_w = ((y1[i][j-1]-y[i-1][j-1])*v_old[i][j-1]+(y[i][j-1]-y1[i][j-1])*v_old[i-1][j-1])/(y[i][j-1]-y[i-1][j-1]);
								ve = ((x1[i][j+1]-x[i][j])*vf_old+(x[i][j]-x1[i][j])*vf_e)/(x1[i][j+1]-x1[i][j]);
								vw = ((x1[i][j]-x[i][j-1])*vf_w+(x[i][j-1]-x1[i][j-1])*vf_old)/(x1[i][j]-x1[i][j-1]);
								temp_Lfc = temp_Lfc + ((-0.5)*(((u_old[i][j]*ve-u_old[i][j-1]*vw)/(x[i][j]-x[i][j-1]))+((v_old[i][j]*v_old[i][j]-v_old[i-1][j]*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
							}
							/*else
							{
								cout<<"wrong no.2L"<<endl;
							}*/
						}
					}
				}
			}
            temp_Lf[k_f]=temp_Lfc;
        }
        temp_L=0.0;
        //#pragma omp parallel for default(shared) private(vf_old,vf_e,vf_w,ve,vw) reduction(+:temp_L) schedule(dynamic)
        for(int i=c33_f;i<c44_f;i++)
        {
            for (int j=c11_f;j<c22_f;j++)
            {
				if(flag[i][j]==0)
				{
					vf_old = ((y1[i][j]-y[i-1][j])*v_old[i][j]+(y[i][j]-y1[i][j])*v_old[i-1][j])/(y[i][j]-y[i-1][j]);
					vf_e = ((y1[i][j+1]-y[i-1][j+1])*v_old[i][j+1]+(y[i][j+1]-y1[i][j+1])*v_old[i-1][j+1])/(y[i][j+1]-y[i-1][j+1]);
					vf_w = ((y1[i][j-1]-y[i-1][j-1])*v_old[i][j-1]+(y[i][j-1]-y1[i][j-1])*v_old[i-1][j-1])/(y[i][j-1]-y[i-1][j-1]);
					ve = ((x1[i][j+1]-x[i][j])*vf_old+(x[i][j]-x1[i][j])*vf_e)/(x1[i][j+1]-x1[i][j]);
					vw = ((x1[i][j]-x[i][j-1])*vf_w+(x[i][j-1]-x1[i][j-1])*vf_old)/(x1[i][j]-x1[i][j-1]);
					temp_L += ((-0.5)*(((u_old[i][j]*ve-u_old[i][j-1]*vw)/(x[i][j]-x[i][j-1]))+((v_old[i][j]*v_old[i][j]-v_old[i-1][j]*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
				}
            }
        }
        //storing the values of u and v
        for (int i=0; i<n_u; i++)
        {
            for(int j=0; j<m_u; j++)
            {
                u_old[i][j]= u[i][j];
            }
        }
        for (int i=0; i<n_v; i++)
        {
            for(int j=0; j<m_v; j++)
            {
                v_old[i][j]= v[i][j];
            }
        }
        /*
         ofstream filew6("a_u.txt");
         if (filew6.is_open())
         {
         for(int i=0; i<n_u; i++)
         {
         for(int j=0; j<m_u; j++)
         {
         filew6 << a_u[i][j] << " ";
         }
         filew6 << endl;
         }
         filew6.close();
         }
         ofstream filew7("b_u.txt");
         if (filew7.is_open())
         {
         for(int i=0; i<n_u; i++)
         {
         for(int j=0; j<m_u; j++)
         {
         filew7 << b_u[i][j] << " ";
         }
         filew7 << endl;
         }
         filew7.close();
         }
         ofstream filew8("c_u.txt");
         if (filew8.is_open())
         {
         for(int i=0; i<n_u; i++)
         {
         for(int j=0; j<m_u; j++)
         {
         filew8 << c_u[i][j] << " ";
         }
         filew8 << endl;
         }
         filew8.close();
         }
         ofstream filew9("d_u.txt");
         if (filew9.is_open())
         {
         for(int i=0; i<n_u; i++)
         {
         for(int j=0; j<m_u; j++)
         {
         filew9 << d_u[i][j] << " ";
         }
         filew9 << endl;
         }
         filew9.close();
         }
         ofstream filew10("e_u.txt");
         if (filew10.is_open())
         {
         for(int i=0; i<n_u; i++)
         {
         for(int j=0; j<m_u; j++)
         {
         filew10 << e_u[i][j] << " ";
         }
         filew10 << endl;
         }
         filew10.close();
         }
         ofstream filew11("f_u.txt");
         if (filew11.is_open())
         {
         for(int i=0; i<n_u; i++)
         {
         for(int j=0; j<m_u; j++)
         {
         filew11 << f_u[i][j] << " ";
         }
         filew11 << endl;
         }
         filew11.close();
         }
         ofstream filew17("f_v.txt");
         if (filew17.is_open())
         {
         for(int i=0; i<n_v; i++)
         {
         for(int j=0; j<m_v; j++)
         {
         filew17 << f_v[i][j] << " ";
         }
         filew17 << endl;
         }
         filew17.close();
         }
         */
#pragma acc data copy(u_red, u_black, u_red_old, u_black_old, uA_old,v_red, v_black, v_red_old, v_black_old, vA_old) copy(u,v, pressure) copyin(a_u, b_u, c_u, d_u, e_u, Bf_u, f_u,x1, x, x_s) copyin(a_v, b_v, c_v, d_v, e_v, Bf_v, f_v,y1, y, y_s)  copyin(flag_u, flag_tu, flag_v, flag_tv) copyin(u_c,v_c,u_c_rigid,v_c_rigid) copy(/*u_red, u_black, u_red_old, u_black_old, uA_old,v_red, v_black, v_red_old, v_black_old, vA_old,*/ f_p, p_c, p_c_red, p_c_red_old, p_c_black, p_c_black_old, pA_old) copyin(a_p, b_p, c_p, d_p, e_p) copy(c1, c2, c3, c4, c11, c22, c33, c44, c11_f, c22_f, c33_f, c44_f)
{
        sum = 0.0;
        #pragma acc parallel loop reduction(+:sum)
        for(int i=1;i<n_u-1;i++)
        {
            sum = sum + (u[i][0]*(y[i][m_u-1]-y[i-1][m_u-1]));
        }
        // Solution of the x and y momentum equation
        for (int p=0; p<GSite; p++)
        {
            #pragma acc parallel loop collapse(2)
            for(int i=0;i<n_u;i++)
            {
                for(int j=0;j<m_u;j++)
                {
                    uA_old[i][j]=u[i][j];
                    u_red_old[i][j]=u_red[i][j];
                    u_black_old[i][j]=u_black[i][j];
                }
            }
            //#pragma omp parallel for default(shared) schedule(dynamic)
            #pragma acc parallel loop
            for(int i=0; i<n_u; i++)
            {
                int j1;
                if(i%2==0)
                     j1=1;
                else
                     j1=2;
                #pragma acc loop
                for (int j=j1; j<m_u; j+=2)
                {
                    if(i==0)
                    {
                        u_red[i][j] = u_black[i+1][j];
                    }
                    else if(i==n_u-1)
                    {
                        u_red[i][j] = u_black[i-1][j];
                    }
                    else
                    {
                        if(j==m_u-1)
                        /*u[i][j] = (u_old[i][j]+((c_sasv*dt)/(x[0][j]-x[0][j-1]))*u[i][j-1])/(1+((c_sasv*dt)/(x[0][j]-x[0][j-1])));*/
                            u_red[i][j]=u_black[i][j-1];
                        else
                            u_red[i][j] = alpha_SOR*((b_u[i][j]*u_black[i+1][j] + c_u[i][j]*u_black[i-1][j] + d_u[i][j]*u_black[i][j+1] + e_u[i][j]*u_black[i][j-1] + f_u[i][j] + (dt/(x1[i][j+1]-x1[i][j]))*(pressure[i][j]-pressure[i][j+1]) + dt*Bf_u[i][j])/a_u[i][j])+(1-alpha_SOR)*u_red_old[i][j];
                    }
                }
            }
            //#pragma omp parallel for default(shared) schedule(dynamic)
            #pragma acc parallel loop
            for(int i=0; i<n_u; i++)
            {
                int j1;
                if(i%2==0)
                     j1=2;
                else
                     j1=1;
                #pragma acc loop
                for (int j=j1; j<m_u; j+=2)
                {
                    if(i==0)
                    {
                        u_black[i][j] = u_red[i+1][j];
                    }
                    else if(i==n_u-1)
                    {
                        u_black[i][j] = u_red[i-1][j];
                    }
                    else
                    {
                        if(j==m_u-1)
                        /*u[i][j] = (u_old[i][j]+((c_sasv*dt)/(x[0][j]-x[0][j-1]))*u[i][j-1])/(1+((c_sasv*dt)/(x[0][j]-x[0][j-1])));*/
                            u_black[i][j]=u_red[i][j-1];
                        else
                            u_black[i][j] = alpha_SOR*((b_u[i][j]*u_red[i+1][j] + c_u[i][j]*u_red[i-1][j] + d_u[i][j]*u_red[i][j+1] + e_u[i][j]*u_red[i][j-1] + f_u[i][j] + (dt/(x1[i][j+1]-x1[i][j]))*(pressure[i][j]-pressure[i][j+1]) + dt*Bf_u[i][j])/a_u[i][j])+(1-alpha_SOR)*u_black_old[i][j];
                    }
                }
            }

            #pragma acc parallel loop
            for(int i=0;i<n_u;i++)
            {
                int j1;
                if(i%2==0)
                     j1=1;
                else
                     j1=2;
                #pragma acc loop
                for(int j=j1;j<m_u;j+=2)
                {
                    u[i][j]=u_red[i][j];
                }
                if(i%2==0)
                     j1=2;
                else
                     j1=1;
                #pragma acc loop
                for(int j=j1;j<m_u;j+=2)
                {
                    u[i][j]=u_black[i][j];
                }
            }
            for(int ctr=1;ctr<=40;ctr++)
            {
                sum1 = 0.0;
                #pragma acc parallel loop reduction(+:sum1)
                for(int i=1;i<n_u-1;i++)
                {
                    sum1 = sum1 + (u[i][m_u-1]*(y[i][m_u-1]-y[i-1][m_u-1]));
                }
                #pragma acc parallel loop
                for(int i=1; i<n_u-1; i++)
                {
                    u[i][m_u-1] = u[i][m_u-1]+((sum-sum1)/sum2);
                }
            }

            big=0;
            #pragma acc parallel loop reduction(max:big)
            for(int i=0;i<n_u;i++)
            {
                for(int j=0;j<m_u;j++)
                {
                    error=abs(u[i][j]-uA_old[i][j]);
                    if (error>big)
                        big = error;
                }
            }
            max_error1[p]=big;

            #pragma acc parallel loop collapse(2)
            for(int i=0;i<n_v;i++)
            {
                for(int j=0;j<m_v;j++)
                {
                    vA_old[i][j]=v[i][j];
                    v_red_old[i][j]=v_red[i][j];
                    v_black_old[i][j]=v_black[i][j];
                }
            }
            //#pragma omp parallel for default(shared) schedule(dynamic)
            #pragma acc parallel loop
            for(int i=1; i<n_v-1; i++)
            {
                int j1;
                if(i%2==0)
                    j1=1;
                else
                    j1=2;
                #pragma acc loop
                for (int j=j1; j<=m_v-1; j+=2)
                {
                    if(j==m_v-1)
                    {
                        /* v[i][j] = (v_old[i][j]+((c_sasv*dt)/(x[0][j]-x[0][j-1]))*v[i][j-1])/(1+((c_sasv*dt)/(x[0][j]-x[0][j-1])));*/
                        v_red[i][j]=v_black[i][j-1];
                    }
                    else
                    {
                        v_red[i][j] = alpha_SOR*((b_v[i][j]*v_black[i+1][j] + c_v[i][j]*v_black[i-1][j] + d_v[i][j]*v_black[i][j+1] + e_v[i][j]*v_black[i][j-1] + f_v[i][j] + (dt/(y1[i+1][j]-y1[i][j]))*(pressure[i][j]-pressure[i+1][j]) + dt*Bf_v[i][j])/a_v[i][j])+(1-alpha_SOR)*v_red_old[i][j];
                    }
                }
            }
            //#pragma omp parallel for default(shared) schedule(dynamic)
            #pragma acc parallel loop
            for(int i=1; i<n_v-1; i++)
            {
                int j1;
                if(i%2==0)
                    j1=2;
                else
                    j1=1;
                #pragma acc loop
                for (int j=j1; j<=m_v-1; j+=2)
                {
                    if(j==m_v-1)
                    {
                        /* v[i][j] = (v_old[i][j]+((c_sasv*dt)/(x[0][j]-x[0][j-1]))*v[i][j-1])/(1+((c_sasv*dt)/(x[0][j]-x[0][j-1])));*/
                        v_black[i][j]=v_red[i][j-1];
                    }
                    else
                    {
                        v_black[i][j] = alpha_SOR*((b_v[i][j]*v_red[i+1][j] + c_v[i][j]*v_red[i-1][j] + d_v[i][j]*v_red[i][j+1] + e_v[i][j]*v_red[i][j-1] + f_v[i][j] + (dt/(y1[i+1][j]-y1[i][j]))*(pressure[i][j]-pressure[i+1][j]) + dt*Bf_v[i][j])/a_v[i][j])+(1-alpha_SOR)*v_black_old[i][j];
                    }
                }
            }
            #pragma acc parallel loop
            for(int i=0;i<n_v;i++)
            {
                int j1;
                if(i%2==0)
                    j1=1;
                else
                    j1=2;
                #pragma acc loop
                for(int j=j1;j<m_v;j+=2)
                {
                    v[i][j]=v_red[i][j];
                }
                if(i%2==0)
                    j1=2;
                else
                    j1=1;
                #pragma acc loop
                for(int j=j1;j<m_v;j+=2)
                {
                    v[i][j]=v_black[i][j];
                }
            }

            big=0;
            #pragma acc parallel loop reduction(max:big)
            for(int i=0;i<n_v;i++)
            {
                for(int j=0;j<m_v;j++)
                {
                    error=abs(v[i][j]-vA_old[i][j]);
                    if (error>big)
                        big = error;
                }
            }
            max_error2[p]=big;

            if (max_error1[p] <= tol && max_error2[p] <= tol)
            {
                break;
            }
        }
        for(int ctr=1;ctr<=40;ctr++)
        {
            sum1 = 0.0;
            #pragma acc parallel loop reduction(+:sum1)
            for(int i=1;i<n_u-1;i++)
            {
                sum1 = sum1 + (u[i][m_u-1]*(y[i][m_u-1]-y[i-1][m_u-1]));
            }
            #pragma acc parallel loop
            for(int i=1; i<n_u-1; i++)
            {
                u[i][m_u-1] = u[i][m_u-1]+((sum-sum1)/sum2);
            }
        }
        /*sum1 = 0.0;
         for(int i=1;i<n_u-1;i++)
         {
         sum1 = sum1 + (u[i][m_u-1]*(y[i][m_u-1]-y[i-1][m_u-1]));
         }
         cout << sum << "  " << (sum-sum1) << "  " << sum1 << endl << endl;*/
        //Calculation of f_p term of the pressure correction equation
//        #pragma omp parallel for default(shared) private(min_dis,k1,u_body,v_body,distance,x_intersect,y_intersect,slope) schedule(dynamic)
        #pragma acc parallel loop collapse(2) private(k1, min_dis, slope)
        for(int i=1; i<n-1; i++)
        {
            for (int j=1; j<m-1; j++)
            {
                if(flag_u[i][j]==1 && flag_u[i][j-1]==1 && flag_v[i][j]==1 && flag_v[i-1][j]==1)
                    f_p[i][j] = -(1.0/dt)*(((u[i][j]-u[i][j-1])/(x[i][j]-x[i][j-1]))+((v[i][j]-v[i-1][j])/(y[i][j]-y[i-1][j])));
                else if (flag_u[i][j]==0 && flag_u[i][j-1]==0 && flag_v[i][j]==0 && flag_v[i-1][j]==0)
                    f_p[i][j] = 0.0;
                else
                {
                    min_dis = 10.0;
                    //#pragma acc loop
                    for(int ctr=0;ctr<N_marker;ctr++)
                    {
                        distance = sqrt((x_s[ctr]-x1[i][j])*(x_s[ctr]-x1[i][j]) + (y_s[ctr]-y1[i][j])*(y_s[ctr]-y1[i][j]));
                        if (distance < min_dis)
                        {
                            min_dis = distance;
                            k1 = ctr;
                        }
                    }
                    if (k1==0)
                    {
                        if((x_s[k1+1]-x_s[N_marker-1])==0)
                        {
                            x_intersect = x_s[N_marker-1];
                            y_intersect = y1[i][j];
                        }
                        else if ((y_s[k1+1]-y_s[N_marker-1])==0)
                        {
                            x_intersect = x1[i][j];
                            y_intersect = y_s[N_marker-1];
                        }
                        else
                        {
                            slope = (y_s[k1+1]-y_s[N_marker-1])/(x_s[k1+1]-x_s[N_marker-1]);
                            x_intersect = (y1[i][j]-y_s[N_marker-1]+slope*x_s[N_marker-1]+(1/slope)*x1[i][j])/(slope+(1/slope));
                            y_intersect = slope*x_intersect + y_s[N_marker-1] - slope*x_s[N_marker-1];
                        }
                    }
                    else if (k1==N_marker-1)
                    {
                        if((x_s[0]-x_s[k1-1])==0)
                        {
                            x_intersect = x_s[k1-1];
                            y_intersect = y1[i][j];
                        }
                        else if ((y_s[0]-y_s[k1-1])==0)
                        {
                            x_intersect = x1[i][j];
                            y_intersect = y_s[k1-1];
                        }
                        else
                        {
                            slope = (y_s[0]-y_s[k1-1])/(x_s[0]-x_s[k1-1]);
                            x_intersect = (y1[i][j]-y_s[k1-1]+slope*x_s[k1-1]+(1/slope)*x1[i][j])/(slope+(1/slope));
                            y_intersect = slope*x_intersect + y_s[k1-1] - slope*x_s[k1-1];
                        }
                    }
                    else
                    {
                        if((x_s[k1+1]-x_s[k1-1])==0)
                        {
                            x_intersect = x_s[k1-1];
                            y_intersect = y1[i][j];
                        }
                        else if ((y_s[k1+1]-y_s[k1-1])==0)
                        {
                            x_intersect = x1[i][j];
                            y_intersect = y_s[k1-1];
                        }
                        else
                        {
                            slope = (y_s[k1+1]-y_s[k1-1])/(x_s[k1+1]-x_s[k1-1]);
                            x_intersect = (y1[i][j]-y_s[k1-1]+slope*x_s[k1-1]+(1/slope)*x1[i][j])/(slope+(1/slope));
                            y_intersect = slope*x_intersect + y_s[k1-1] - slope*x_s[k1-1];
                        }
                    }
                    //u_body = u_c[k1] ;
                    //v_body = v_c[k1] ;
                        if(k1>=0 && k1<(Nf1*nn1)-(nn1-1))
                        {
                                u_body = u_c_rigid ;
                                v_body = v_c_rigid ;
                        }
                        else if(k1>=((Nf2*2)*nn2+(Nf1*nn1)-(nn1-1))-1 && k1<N_marker)
                        {
                                u_body = u_c_rigid ;
                                v_body = v_c_rigid ;
                        }
                        else
                        {
                                theta_intp=atan((y_s[k1+1]-y_s[k1-1])/(x_s[k1+1]-x_s[k1-1])) ;
				if (theta_intp>api/2.0 || theta_intp<-api/2.0)
				{
					r_f0=0.0;
					r_f1=(y_intersect-y_s[k1-1])/sin(theta_intp) ;
					r_f2=(y_s[k1+1]-y_s[k1-1])/sin(theta_intp) ;
				}
				else
				{
					r_f0=0.0;
					r_f1=(x_intersect-x_s[k1-1])/cos(theta_intp) ;
					r_f2=(x_s[k1+1]-x_s[k1-1])/cos(theta_intp) ;
				}
				u_body=((r_f1-r_f0)/(r_f2-r_f0))*u_c[k1+1] + ((r_f2-r_f1)/(r_f2-r_f0))*u_c[k1-1] ;
				v_body=((r_f1-r_f0)/(r_f2-r_f0))*v_c[k1+1] + ((r_f2-r_f1)/(r_f2-r_f0))*v_c[k1-1] ;
                        }
                    
                   	                /*if (k1==0)
					{
						theta_intp=atan((y_s[k1+1]-y_s[N_marker-1])/(x_s[k1+1]-x_s[N_marker-1])) ;
						if (theta_intp>api/2.0 || theta_intp<-api/2.0)
						{
							r_f0=0.0;
							r_f1=(y_intersect-y_s[N_marker-1])/sin(theta_intp) ;
							r_f2=(y_s[k1+1]-y_s[N_marker-1])/sin(theta_intp) ;
						}
						else
						{
							r_f0=0.0;
							r_f1=(x_intersect-x_s[N_marker-1])/cos(theta_intp) ;
							r_f2=(x_s[k1+1]-x_s[N_marker-1])/cos(theta_intp) ;
						}
						u_body=((r_f1-r_f0)/(r_f2-r_f0))*u_c[k1+1] + ((r_f2-r_f1)/(r_f2-r_f0))*u_c[N_marker-1] ;
						v_body=((r_f1-r_f0)/(r_f2-r_f0))*v_c[k1+1] + ((r_f2-r_f1)/(r_f2-r_f0))*v_c[N_marker-1] ;
					}
					else if (k1==N_marker-1)
					{
						theta_intp=atan((y_s[0]-y_s[k1-1])/(x_s[0]-x_s[k1-1])) ;
						if (theta_intp>api/2.0 || theta_intp<-api/2.0)
						{
							r_f0=0.0;
							r_f1=(y_intersect-y_s[k1-1])/sin(theta_intp) ;
							r_f2=(y_s[0]-y_s[k1-1])/sin(theta_intp) ;
						}
						else
						{
							r_f0=0.0;
							r_f1=(x_intersect-x_s[k1-1])/cos(theta_intp) ;
							r_f2=(x_s[0]-x_s[k1-1])/cos(theta_intp) ;
						}
						u_body=((r_f1-r_f0)/(r_f2-r_f0))*u_c[0] + ((r_f2-r_f1)/(r_f2-r_f0))*u_c[k1-1] ;
						v_body=((r_f1-r_f0)/(r_f2-r_f0))*v_c[0] + ((r_f2-r_f1)/(r_f2-r_f0))*v_c[k1-1] ;
					}
					else
					{
						theta_intp=atan((y_s[k1+1]-y_s[k1-1])/(x_s[k1+1]-x_s[k1-1])) ;
						if (theta_intp>api/2.0 || theta_intp<-api/2.0)
						{
							r_f0=0.0;
							r_f1=(y_intersect-y_s[k1-1])/sin(theta_intp) ;
							r_f2=(y_s[k1+1]-y_s[k1-1])/sin(theta_intp) ;
						}
						else
						{
							r_f0=0.0;
							r_f1=(x_intersect-x_s[k1-1])/cos(theta_intp) ;
							r_f2=(x_s[k1+1]-x_s[k1-1])/cos(theta_intp) ;
						}
						u_body=((r_f1-r_f0)/(r_f2-r_f0))*u_c[k1+1] + ((r_f2-r_f1)/(r_f2-r_f0))*u_c[k1-1] ;
						v_body=((r_f1-r_f0)/(r_f2-r_f0))*v_c[k1+1] + ((r_f2-r_f1)/(r_f2-r_f0))*v_c[k1-1] ;
					}*/
                    f_p[i][j] = -(1.0/dt)*((((u[i][j]-u[i][j-1])/(x[i][j]-x[i][j-1]))+((v[i][j]-v[i-1][j])/(y[i][j]-y[i-1][j]))) - (((flag_tu[i][j]*(u[i][j]-u_body))/(x[i][j]-x[i][j-1]))+(-(flag_tu[i][j-1]*(u[i][j-1]-u_body))/(x[i][j]-x[i][j-1]))+((flag_tv[i][j]*(v[i][j]-v_body))/(y[i][j]-y[i-1][j]))+(-(flag_tv[i-1][j]*(v[i-1][j]-v_body))/(y[i][j]-y[i-1][j]))));
                }
            }
        }
        /*   ofstream filew23("f_p.txt");
         if (filew23.is_open())
         {
         for(int i=0; i<n; i++)
         {
         for(int j=0; j<m; j++)
         {
         filew23 << f_p[i][j]*dt << " ";
         }
         filew23 << endl;
         }
         filew23.close();
         }
         ofstream filew1("Cavity_u.txt");
         if (filew1.is_open())
         {
         for(int i=0; i<n_u; i++)
         {
         for(int j=0; j<m_u; j++)
         {
         filew1 << u[i][j] << " ";
         }
         filew1 << endl;
         }
         filew1.close();
         }
         ofstream filew2("Cavity_v.txt");
         if (filew2.is_open())
         {
         for(int i=0; i<n_v; i++)
         {
         for(int j=0; j<m_v; j++)
         {
         filew2 << v[i][j] << " ";
         }
         filew2 << endl;
         }
         filew2.close();
         }*/
        /*big = 0;
        #pragma acc parallel loop reduction(max:big)
        for (int i=1;i<n-1;i++)
        {
            for (int j=1;j<m-1;j++)
            {
                error=abs(f_p[i][j]);
                if (error>big)
                    big = error;
            }
        }
        residual_flux=big*dt;*/

        for (int p=0; p<GSite2; p++)
        {
            #pragma acc parallel loop collapse(2)
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<m;j++)
                {
                    pA_old[i][j]=p_c[i][j];
                    p_c_red_old[i][j]=p_c_red[i][j];
                    p_c_black_old[i][j]=p_c_black[i][j];
                }
            }
            //#pragma omp parallel for default(shared) schedule(dynamic)
            #pragma acc parallel loop
            for(int i=0; i<n; i++)
            {
                int j1;
                if(i%2==0)
                    j1=0;
                else
                    j1=1;
                #pragma acc loop
                for (int j=j1; j<m; j+=2)
                {
                    if (i==0)//bottom row
                    {
                        p_c_red[i][j] = p_c_black[i+1][j];
                    }
                    else if (i>0 && i<n-1 && j==0)//left column
                    {
                        p_c_red[i][j] = p_c_black[i][j+1];
                    }
                    else if (i>0 && i<n-1 && j==m-1)//right column
                    {
                        //p_c[i][j] = 0.0f;
                        p_c_red[i][j] = p_c_black[i][j-1];
                    }
                    else if (i==n-1)//top row
                    {
                        p_c_red[i][j] = p_c_black[i-1][j];
                    }
                    else
                    {
                        p_c_red[i][j] = alpha_SOR*((b_p[i][j]*p_c_black[i+1][j] + c_p[i][j]*p_c_black[i-1][j] + d_p[i][j]*p_c_black[i][j+1] + e_p[i][j]*p_c_black[i][j-1] + f_p[i][j])/a_p[i][j]) + (1-alpha_SOR)*p_c_red_old[i][j];
                    }
                }
            }
            //#pragma omp parallel for default(shared) schedule(dynamic)
            #pragma acc parallel loop
            for(int i=0; i<n; i++)
            {
                int j1;
                if(i%2==0)
                    j1=1;
                else
                    j1=0;
                #pragma acc loop
                for (int j=j1; j<m; j+=2)
                {
                    if (i==0)//bottom row
                    {
                        p_c_black[i][j] = p_c_red[i+1][j];
                    }
                    else if (i>0 && i<n-1 && j==0)//left column
                    {
                        p_c_black[i][j] = p_c_red[i][j+1];
                    }
                    else if (i>0 && i<n-1 && j==m-1)//right column
                    {
                        //p_c[i][j] = 0.0f;
                        p_c_black[i][j] = p_c_red[i][j-1];
                    }
                    else if (i==n-1)//top row
                    {
                        p_c_black[i][j] = p_c_red[i-1][j];
                    }
                    else
                    {
                        p_c_black[i][j] = alpha_SOR*((b_p[i][j]*p_c_red[i+1][j] + c_p[i][j]*p_c_red[i-1][j] + d_p[i][j]*p_c_red[i][j+1] + e_p[i][j]*p_c_red[i][j-1] + f_p[i][j])/a_p[i][j]) + (1-alpha_SOR)*p_c_black_old[i][j];
                    }
                }
            }
            #pragma acc parallel loop
            for(int i=0;i<n;i++)
            {
                int j1;
                if(i%2==0)
                    j1=0;
                else
                    j1=1;
                #pragma acc loop
                for (int j=j1; j<m; j+=2)
                {
                    p_c[i][j]=p_c_red[i][j];
                }
                if(i%2==0)
                    j1=1;
                else
                    j1=0;
                #pragma acc loop
                for (int j=j1; j<m; j+=2)
                {
                    p_c[i][j]=p_c_black[i][j];
                }
            }

            big=0;
            #pragma acc parallel loop collapse(2) reduction(max:big)
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<m;j++)
                {
                    error=abs(p_c[i][j]-pA_old[i][j]);
                    if (error>big)
                        big = error;
                }
            }
            max_error3[p]=big;

            if (max_error3[p] <= tol)
            {
                break;
            }
        }
        //pressure correction
        #pragma acc parallel loop collapse(2)
        for(int i=1; i<n-1; i++)
        {
            for (int j=1; j<m-1; j++)
            {
                pressure[i][j] = pressure[i][j]+ 0.0001*(p_c[i][j] - (dt/(2*Re))*(((((p_c[i][j+1]-p_c[i][j])/(x1[i][j+1]-x1[i][j]))-((p_c[i][j]-p_c[i][j-1])/(x1[i][j]-x1[i][j-1])))/(x[i][j]-x[i][j-1]))+((((p_c[i+1][j]-p_c[i][j])/(y1[i+1][j]-y1[i][j]))-((p_c[i][j]-p_c[i-1][j])/(y1[i][j]-y1[i-1][j])))/(y[i][j]-y[i-1][j]))));
            }
        }
        //u velocity correction
        #pragma acc parallel loop collapse(2)
        for(int i=1; i<n_u-1; i++)
        {
            for (int j=1; j<m_u-1; j++)
            {
                u[i][j] = u[i][j]-dt*((p_c[i][j+1]-p_c[i][j])/(x1[i][j+1]-x1[i][j]));
            }
        }
        //v velocity correction
        #pragma acc parallel loop collapse(2)
        for(int i=1; i<n_v-1; i++)
        {
            for (int j=1; j<m_v-1; j++)
            {
                v[i][j] = v[i][j]-dt*((p_c[i+1][j]-p_c[i][j])/(y1[i+1][j]-y1[i][j]));
            }
        }


        }  // End of GPU
        for (int p=0; p<4; p++)
        {
            for(int i=0; i<n_u; i++)
            {
                for (int j=1; j<m_u; j++)
                {
                    if(i==0)
                    {
                        u[i][j] = u[i+1][j];
                    }
                    else if(i==n_u-1)
                    {
                        u[i][j] = u[i-1][j];
                    }
                    else if(j==m_u-1)
                    {
                        /*u[i][j] = (u_old[i][j]+((c_sasv*dt)/(x[0][j]-x[0][j-1]))*u[i][j-1])/(1+((c_sasv*dt)/(x[0][j]-x[0][j-1])));*/
                        u[i][j]=u[i][j-1];
                    }
                }
            }
            for(int ctr=1;ctr<=30;ctr++)
            {
                sum1 = 0.0;
                for(int i=1;i<n_u-1;i++)
                {
                    sum1 = sum1 + (u[i][m_u-1]*(y[i][m_u-1]-y[i-1][m_u-1]));
                }
                for(int i=1; i<n_u-1; i++)
                {
                    u[i][m_u-1] = u[i][m_u-1]+((sum-sum1)/sum2);
                }
            }
            for(int i=1; i<n_v-1; i++)
            {
                for (int j=1; j<=m_v-1; j++)
                {
                    if(j==m_v-1)
                    {
                        /*v[i][j] = (v_old[i][j]+((c_sasv*dt)/(x[0][j]-x[0][j-1]))*v[i][j-1])/(1+((c_sasv*dt)/(x[0][j]-x[0][j-1])));*/
                       v[i][j]=v[i][j-1];
                    }
                }
            }
        }
        /*sum=0;
         for(int i=1;i<n_u-1;i++)
         {
         sum = sum + u[i][m_u-1]*(y[i][0]-y[i-1][0]);
         }
         sum1=0;
         for(int i=1;i<n_u-1;i++)
         {
         sum1 = sum1 + u[i][0]*(y[i][0]-y[i-1][0]);
         }
         cout << c_sasv << "  " << sum << "  " << sum1 << endl;*/

		for(int k_f=0; k_f<Nf; k_f++)
        {
            tempfc=0.0;
			minfx1=min(x_sfc[k_f],x_sfc[k_f+1]);
			minfx2=min(x_sfc[N_markerf-k_f-1],x_sfc[N_markerf-k_f-2]);
			maxfx1=max(x_sfc[k_f],x_sfc[k_f+1]);
			maxfx2=max(x_sfc[N_markerf-k_f-1],x_sfc[N_markerf-k_f-2]);
			minfy1=min(y_sfc[k_f],y_sfc[k_f+1]);
			minfy2=min(y_sfc[N_markerf-k_f-1],y_sfc[N_markerf-k_f-2]);
			maxfy1=max(y_sfc[k_f],y_sfc[k_f+1]);
			maxfy2=max(y_sfc[N_markerf-k_f-1],y_sfc[N_markerf-k_f-2]);
            cf1=m_div1+ceil((min(minfx1,minfx2)-Hx1)/dx2)-10;
			cf2=m_div1+ceil((max(maxfx1,maxfx2)-Hx1)/dx2)+10;
			cf3=n_div1+ceil((min(minfy1,minfy2)-Hy1)/dy2)-10;
			cf4=n_div1+ceil((max(maxfy1,maxfy2)-Hy1)/dy2)+10;

			if ((y_sfc[k_f]-y_sfc[N_markerf-k_f-1])>=0.0)
            {
                Afc1=(y_sfc[k_f]-y_sfc[N_markerf-k_f-1]);
                Afc2=(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2]);
                Bfc1=-(x_sfc[k_f]-x_sfc[N_markerf-k_f-1]);
                Bfc2=-(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2]);;
                Cfc1=((y_sfc[N_markerf-k_f-1]*(x_sfc[k_f]-x_sfc[N_markerf-k_f-1])) - (x_sfc[N_markerf-k_f-1]*(y_sfc[k_f]-y_sfc[N_markerf-k_f-1])));
                Cfc2=((y_sfc[N_markerf-k_f-2]*(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2])) - (x_sfc[N_markerf-k_f-2]*(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2])));
            }
            else
            {
                Afc1=-(y_sfc[k_f]-y_sfc[N_markerf-k_f-1]);
                Afc2=-(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2]);
                Bfc1=(x_sfc[k_f]-x_sfc[N_markerf-k_f-1]);
                Bfc2=(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2]);;
                Cfc1=-((y_sfc[N_markerf-k_f-1]*(x_sfc[k_f]-x_sfc[N_markerf-k_f-1])) - (x_sfc[N_markerf-k_f-1]*(y_sfc[k_f]-y_sfc[N_markerf-k_f-1])));
                Cfc2=-((y_sfc[N_markerf-k_f-2]*(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2])) - (x_sfc[N_markerf-k_f-2]*(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2])));
            }
			//cout<< cf1 << " " <<cf2 << " " << cf3 << " " << cf4 <<endl;
			for (int i=cf3;i<cf4;i++)
			{
				for (int j=cf1;j<cf2;j++)
				{
					if (flag[i][j]==0)
					{
						if (thetafc[k_f]!=-api && (xf[k_f]-xf[k_f-1])>=0.0 && (Afc1*x1[i][j]+Bfc1*y1[i][j]+Cfc1)>=0.0)
						{
							if (thetafc[k_f+1]!=-api && (xf[k_f+1]-xf[k_f])>=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)<=0.0)
							{
								uf = ((x1[i][j]-x[i][j-1])*u[i][j]+(x[i][j]-x1[i][j])*u[i][j-1])/(x[i][j]-x[i][j-1]);
								uf_old = ((x1[i][j]-x[i][j-1])*u_old[i][j]+(x[i][j]-x1[i][j])*u_old[i][j-1])/(x[i][j]-x[i][j-1]);
								uf_n = ((x1[i+1][j]-x[i+1][j-1])*u_old[i+1][j]+(x[i+1][j]-x1[i+1][j])*u_old[i+1][j-1])/(x[i+1][j]-x[i+1][j-1]);
								uf_s = ((x1[i-1][j]-x[i-1][j-1])*u_old[i-1][j]+(x[i-1][j]-x1[i-1][j])*u_old[i-1][j-1])/(x[i-1][j]-x[i-1][j-1]);
								un = ((y[i][j]-y1[i][j])*uf_n+(y1[i+1][j]-y[i][j])*uf_old)/(y1[i+1][j]-y1[i][j]);
								us = ((y[i-1][j]-y1[i-1][j])*uf_old+(y1[i][j]-y[i-1][j])*uf_s)/(y1[i][j]-y1[i-1][j]);
								tempfc = tempfc+ (((uf-uf_old)/dt)+1.5*(((u_old[i][j]*u_old[i][j]-u_old[i][j-1]*u_old[i][j-1])/(x[i][j]-x[i][j-1]))+((un*v_old[i][j]-us*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
							}
							else if (thetafc[k_f+1]!=api && (xf[k_f+1]-xf[k_f])<=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)>=0.0)
							{
								uf = ((x1[i][j]-x[i][j-1])*u[i][j]+(x[i][j]-x1[i][j])*u[i][j-1])/(x[i][j]-x[i][j-1]);
								uf_old = ((x1[i][j]-x[i][j-1])*u_old[i][j]+(x[i][j]-x1[i][j])*u_old[i][j-1])/(x[i][j]-x[i][j-1]);
								uf_n = ((x1[i+1][j]-x[i+1][j-1])*u_old[i+1][j]+(x[i+1][j]-x1[i+1][j])*u_old[i+1][j-1])/(x[i+1][j]-x[i+1][j-1]);
								uf_s = ((x1[i-1][j]-x[i-1][j-1])*u_old[i-1][j]+(x[i-1][j]-x1[i-1][j])*u_old[i-1][j-1])/(x[i-1][j]-x[i-1][j-1]);
								un = ((y[i][j]-y1[i][j])*uf_n+(y1[i+1][j]-y[i][j])*uf_old)/(y1[i+1][j]-y1[i][j]);
								us = ((y[i-1][j]-y1[i-1][j])*uf_old+(y1[i][j]-y[i-1][j])*uf_s)/(y1[i][j]-y1[i-1][j]);
								tempfc = tempfc+ (((uf-uf_old)/dt)+1.5*(((u_old[i][j]*u_old[i][j]-u_old[i][j-1]*u_old[i][j-1])/(x[i][j]-x[i][j-1]))+((un*v_old[i][j]-us*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
							}
							/*else
							{
								cout<<"wrong no.1Df"<<endl;
							}*/
						}
						else if (thetafc[k_f]!=api && (xf[k_f]-xf[k_f-1])<=0.0 && (Afc1*x1[i][j]+Bfc1*y1[i][j]+Cfc1)<=0.0)
						{
							if (thetafc[k_f+1]!=-api && (xf[k_f+1]-xf[k_f])>=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)<=0.0)
							{
								uf = ((x1[i][j]-x[i][j-1])*u[i][j]+(x[i][j]-x1[i][j])*u[i][j-1])/(x[i][j]-x[i][j-1]);
								uf_old = ((x1[i][j]-x[i][j-1])*u_old[i][j]+(x[i][j]-x1[i][j])*u_old[i][j-1])/(x[i][j]-x[i][j-1]);
								uf_n = ((x1[i+1][j]-x[i+1][j-1])*u_old[i+1][j]+(x[i+1][j]-x1[i+1][j])*u_old[i+1][j-1])/(x[i+1][j]-x[i+1][j-1]);
								uf_s = ((x1[i-1][j]-x[i-1][j-1])*u_old[i-1][j]+(x[i-1][j]-x1[i-1][j])*u_old[i-1][j-1])/(x[i-1][j]-x[i-1][j-1]);
								un = ((y[i][j]-y1[i][j])*uf_n+(y1[i+1][j]-y[i][j])*uf_old)/(y1[i+1][j]-y1[i][j]);
								us = ((y[i-1][j]-y1[i-1][j])*uf_old+(y1[i][j]-y[i-1][j])*uf_s)/(y1[i][j]-y1[i-1][j]);
								tempfc = tempfc+ (((uf-uf_old)/dt)+1.5*(((u_old[i][j]*u_old[i][j]-u_old[i][j-1]*u_old[i][j-1])/(x[i][j]-x[i][j-1]))+((un*v_old[i][j]-us*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
							}
							else if (thetafc[k_f+1]!=api && (xf[k_f+1]-xf[k_f])<=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)>=0.0)
							{
								uf = ((x1[i][j]-x[i][j-1])*u[i][j]+(x[i][j]-x1[i][j])*u[i][j-1])/(x[i][j]-x[i][j-1]);
								uf_old = ((x1[i][j]-x[i][j-1])*u_old[i][j]+(x[i][j]-x1[i][j])*u_old[i][j-1])/(x[i][j]-x[i][j-1]);
								uf_n = ((x1[i+1][j]-x[i+1][j-1])*u_old[i+1][j]+(x[i+1][j]-x1[i+1][j])*u_old[i+1][j-1])/(x[i+1][j]-x[i+1][j-1]);
								uf_s = ((x1[i-1][j]-x[i-1][j-1])*u_old[i-1][j]+(x[i-1][j]-x1[i-1][j])*u_old[i-1][j-1])/(x[i-1][j]-x[i-1][j-1]);
								un = ((y[i][j]-y1[i][j])*uf_n+(y1[i+1][j]-y[i][j])*uf_old)/(y1[i+1][j]-y1[i][j]);
								us = ((y[i-1][j]-y1[i-1][j])*uf_old+(y1[i][j]-y[i-1][j])*uf_s)/(y1[i][j]-y1[i-1][j]);
								tempfc = tempfc+ (((uf-uf_old)/dt)+1.5*(((u_old[i][j]*u_old[i][j]-u_old[i][j-1]*u_old[i][j-1])/(x[i][j]-x[i][j-1]))+((un*v_old[i][j]-us*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
							}
							/*else
							{
								cout<<"wrong no.2Df"<<endl;
							}*/
						}
					}
				}
			}
            tempf[k_f]=tempfc;
        }
        temp=0.0;
        //#pragma omp parallel for default(shared) private(uf,uf_old,uf_n,uf_s,un,us) reduction(+:temp) schedule(dynamic)
        for(int i=c33_f;i<c44_f;i++)
        {
            for (int j=c11_f;j<c22_f;j++)
            {
				if(flag[i][j]==0)
				{
					uf = ((x1[i][j]-x[i][j-1])*u[i][j]+(x[i][j]-x1[i][j])*u[i][j-1])/(x[i][j]-x[i][j-1]);
					uf_old = ((x1[i][j]-x[i][j-1])*u_old[i][j]+(x[i][j]-x1[i][j])*u_old[i][j-1])/(x[i][j]-x[i][j-1]);
					uf_n = ((x1[i+1][j]-x[i+1][j-1])*u_old[i+1][j]+(x[i+1][j]-x1[i+1][j])*u_old[i+1][j-1])/(x[i+1][j]-x[i+1][j-1]);
					uf_s = ((x1[i-1][j]-x[i-1][j-1])*u_old[i-1][j]+(x[i-1][j]-x1[i-1][j])*u_old[i-1][j-1])/(x[i-1][j]-x[i-1][j-1]);
					un = ((y[i][j]-y1[i][j])*uf_n+(y1[i+1][j]-y[i][j])*uf_old)/(y1[i+1][j]-y1[i][j]);
					us = ((y[i-1][j]-y1[i-1][j])*uf_old+(y1[i][j]-y[i-1][j])*uf_s)/(y1[i][j]-y1[i-1][j]);
					temp += (((uf-uf_old)/dt)+1.5*(((u_old[i][j]*u_old[i][j]-u_old[i][j-1]*u_old[i][j-1])/(x[i][j]-x[i][j-1]))+((un*v_old[i][j]-us*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
				}
            }
        }
		for(int k_f=0; k_f<Nf; k_f++)
        {
            temp1fc=0.0;
			minfx1=min(x_sfc[k_f],x_sfc[k_f+1]);
			minfx2=min(x_sfc[N_markerf-k_f-1],x_sfc[N_markerf-k_f-2]);
			maxfx1=max(x_sfc[k_f],x_sfc[k_f+1]);
			maxfx2=max(x_sfc[N_markerf-k_f-1],x_sfc[N_markerf-k_f-2]);
			minfy1=min(y_sfc[k_f],y_sfc[k_f+1]);
			minfy2=min(y_sfc[N_markerf-k_f-1],y_sfc[N_markerf-k_f-2]);
			maxfy1=max(y_sfc[k_f],y_sfc[k_f+1]);
			maxfy2=max(y_sfc[N_markerf-k_f-1],y_sfc[N_markerf-k_f-2]);
            cf1=m_div1+ceil((min(minfx1,minfx2)-Hx1)/dx2)-10;
			cf2=m_div1+ceil((max(maxfx1,maxfx2)-Hx1)/dx2)+10;
			cf3=n_div1+ceil((min(minfy1,minfy2)-Hy1)/dy2)-10;
			cf4=n_div1+ceil((max(maxfy1,maxfy2)-Hy1)/dy2)+10;

			if ((y_sfc[k_f]-y_sfc[N_markerf-k_f-1])>=0.0)
            {
                Afc1=(y_sfc[k_f]-y_sfc[N_markerf-k_f-1]);
                Afc2=(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2]);
                Bfc1=-(x_sfc[k_f]-x_sfc[N_markerf-k_f-1]);
                Bfc2=-(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2]);;
                Cfc1=((y_sfc[N_markerf-k_f-1]*(x_sfc[k_f]-x_sfc[N_markerf-k_f-1])) - (x_sfc[N_markerf-k_f-1]*(y_sfc[k_f]-y_sfc[N_markerf-k_f-1])));
                Cfc2=((y_sfc[N_markerf-k_f-2]*(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2])) - (x_sfc[N_markerf-k_f-2]*(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2])));
            }
            else
            {
                Afc1=-(y_sfc[k_f]-y_sfc[N_markerf-k_f-1]);
                Afc2=-(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2]);
                Bfc1=(x_sfc[k_f]-x_sfc[N_markerf-k_f-1]);
                Bfc2=(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2]);;
                Cfc1=-((y_sfc[N_markerf-k_f-1]*(x_sfc[k_f]-x_sfc[N_markerf-k_f-1])) - (x_sfc[N_markerf-k_f-1]*(y_sfc[k_f]-y_sfc[N_markerf-k_f-1])));
                Cfc2=-((y_sfc[N_markerf-k_f-2]*(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2])) - (x_sfc[N_markerf-k_f-2]*(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2])));
            }
			//cout<< cf1 << " " <<cf2 << " " << cf3 << " " << cf4 <<endl;
			for (int i=cf3;i<cf4;i++)
			{
				for (int j=cf1;j<cf2;j++)
				{
					if (flag_u[i][j] == 3)
					{
						if (thetafc[k_f]!=-api && (xf[k_f]-xf[k_f-1])>=0.0 && (Afc1*x1[i][j]+Bfc1*y1[i][j]+Cfc1)>=0.0)
						{
							if (thetafc[k_f+1]!=-api && (xf[k_f+1]-xf[k_f])>=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)<=0.0)
							{
								temp1fc = temp1fc + Bf_u[i][j]*(x1[i][j+1]-x1[i][j])*(y[i][j]-y[i-1][j]);
							}
							else if (thetafc[k_f+1]!=api && (xf[k_f+1]-xf[k_f])<=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)>=0.0)
							{
								temp1fc = temp1fc + Bf_u[i][j]*(x1[i][j+1]-x1[i][j])*(y[i][j]-y[i-1][j]);
							}
							/*else
							{
								cout<<"wrong no.1f"<<endl;
							}*/
						}
						else if (thetafc[k_f]!=api && (xf[k_f]-xf[k_f-1])<=0.0 && (Afc1*x1[i][j]+Bfc1*y1[i][j]+Cfc1)<=0.0)
						{
							if (thetafc[k_f+1]!=-api && (xf[k_f+1]-xf[k_f])>=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)<=0.0)
							{
								temp1fc = temp1fc + Bf_u[i][j]*(x1[i][j+1]-x1[i][j])*(y[i][j]-y[i-1][j]);
							}
							else if (thetafc[k_f+1]!=api && (xf[k_f+1]-xf[k_f])<=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)>=0.0)
							{
								temp1fc = temp1fc + Bf_u[i][j]*(x1[i][j+1]-x1[i][j])*(y[i][j]-y[i-1][j]);
							}
							/*else
							{
								cout<<"wrong no.2f"<<endl;
							}*/
						}
					}
				}
			}
            temp1f[k_f]=temp1fc;
        }
        temp1 = 0.0;
        for(int i=c33_f;i<c44_f;i++)
        {
            for(int j=c11_f;j<c22_f;j++)
            {
                if (flag_u[i][j] == 3)
				{
                    temp1 = temp1 + Bf_u[i][j]*(x1[i][j+1]-x1[i][j])*(y[i][j]-y[i-1][j]);
				}
            }
        }
        for(int k_f=0; k_f<Nf; k_f++)
        {
            drag3f[k_f] = 2.0*(-temp1f[k_f]+tempf[k_f]+temp_Df[k_f]);
        }
        drag3 = 2.0*(-temp1+temp+temp_D); //gives drag coeff Cd = D/(0.5*rho*u^2*(d*1))
		for(int k_f=0; k_f<Nf; k_f++)
        {
            tempfc=0.0;
			minfx1=min(x_sfc[k_f],x_sfc[k_f+1]);
			minfx2=min(x_sfc[N_markerf-k_f-1],x_sfc[N_markerf-k_f-2]);
			maxfx1=max(x_sfc[k_f],x_sfc[k_f+1]);
			maxfx2=max(x_sfc[N_markerf-k_f-1],x_sfc[N_markerf-k_f-2]);
			minfy1=min(y_sfc[k_f],y_sfc[k_f+1]);
			minfy2=min(y_sfc[N_markerf-k_f-1],y_sfc[N_markerf-k_f-2]);
			maxfy1=max(y_sfc[k_f],y_sfc[k_f+1]);
			maxfy2=max(y_sfc[N_markerf-k_f-1],y_sfc[N_markerf-k_f-2]);
            cf1=m_div1+ceil((min(minfx1,minfx2)-Hx1)/dx2)-10;
			cf2=m_div1+ceil((max(maxfx1,maxfx2)-Hx1)/dx2)+10;
			cf3=n_div1+ceil((min(minfy1,minfy2)-Hy1)/dy2)-10;
			cf4=n_div1+ceil((max(maxfy1,maxfy2)-Hy1)/dy2)+10;

			if ((y_sfc[k_f]-y_sfc[N_markerf-k_f-1])>=0.0)
            {
                Afc1=(y_sfc[k_f]-y_sfc[N_markerf-k_f-1]);
                Afc2=(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2]);
                Bfc1=-(x_sfc[k_f]-x_sfc[N_markerf-k_f-1]);
                Bfc2=-(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2]);;
                Cfc1=((y_sfc[N_markerf-k_f-1]*(x_sfc[k_f]-x_sfc[N_markerf-k_f-1])) - (x_sfc[N_markerf-k_f-1]*(y_sfc[k_f]-y_sfc[N_markerf-k_f-1])));
                Cfc2=((y_sfc[N_markerf-k_f-2]*(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2])) - (x_sfc[N_markerf-k_f-2]*(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2])));
            }
            else
            {
                Afc1=-(y_sfc[k_f]-y_sfc[N_markerf-k_f-1]);
                Afc2=-(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2]);
                Bfc1=(x_sfc[k_f]-x_sfc[N_markerf-k_f-1]);
                Bfc2=(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2]);;
                Cfc1=-((y_sfc[N_markerf-k_f-1]*(x_sfc[k_f]-x_sfc[N_markerf-k_f-1])) - (x_sfc[N_markerf-k_f-1]*(y_sfc[k_f]-y_sfc[N_markerf-k_f-1])));
                Cfc2=-((y_sfc[N_markerf-k_f-2]*(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2])) - (x_sfc[N_markerf-k_f-2]*(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2])));
            }
			//cout<< cf1 << " " <<cf2 << " " << cf3 << " " << cf4 <<endl;
			for (int i=cf3;i<cf4;i++)
			{
				for (int j=cf1;j<cf2;j++)
				{
					if (flag[i][j]==0)
					{
						if (thetafc[k_f]!=-api && (xf[k_f]-xf[k_f-1])>=0.0 && (Afc1*x1[i][j]+Bfc1*y1[i][j]+Cfc1)>=0.0)
						{
							if (thetafc[k_f+1]!=-api && (xf[k_f+1]-xf[k_f])>=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)<=0.0)
							{
								vf = ((y1[i][j]-y[i-1][j])*v[i][j]+(y[i][j]-y1[i][j])*v[i-1][j])/(y[i][j]-y[i-1][j]);
								vf_old = ((y1[i][j]-y[i-1][j])*v_old[i][j]+(y[i][j]-y1[i][j])*v_old[i-1][j])/(y[i][j]-y[i-1][j]);
								vf_e = ((y1[i][j+1]-y[i-1][j+1])*v_old[i][j+1]+(y[i][j+1]-y1[i][j+1])*v_old[i-1][j+1])/(y[i][j+1]-y[i-1][j+1]);
								vf_w = ((y1[i][j-1]-y[i-1][j-1])*v_old[i][j-1]+(y[i][j-1]-y1[i][j-1])*v_old[i-1][j-1])/(y[i][j-1]-y[i-1][j-1]);
								ve = ((x1[i][j+1]-x[i][j])*vf_old+(x[i][j]-x1[i][j])*vf_e)/(x1[i][j+1]-x1[i][j]);
								vw = ((x1[i][j]-x[i][j-1])*vf_w+(x[i][j-1]-x1[i][j-1])*vf_old)/(x1[i][j]-x1[i][j-1]);
								tempfc = tempfc + (((vf-vf_old)/dt)+1.5*(((u_old[i][j]*ve-u_old[i][j-1]*vw)/(x[i][j]-x[i][j-1]))+((v_old[i][j]*v_old[i][j]-v_old[i-1][j]*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
							}
							else if (thetafc[k_f+1]!=api && (xf[k_f+1]-xf[k_f])<=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)>=0.0)
							{
								vf = ((y1[i][j]-y[i-1][j])*v[i][j]+(y[i][j]-y1[i][j])*v[i-1][j])/(y[i][j]-y[i-1][j]);
								vf_old = ((y1[i][j]-y[i-1][j])*v_old[i][j]+(y[i][j]-y1[i][j])*v_old[i-1][j])/(y[i][j]-y[i-1][j]);
								vf_e = ((y1[i][j+1]-y[i-1][j+1])*v_old[i][j+1]+(y[i][j+1]-y1[i][j+1])*v_old[i-1][j+1])/(y[i][j+1]-y[i-1][j+1]);
								vf_w = ((y1[i][j-1]-y[i-1][j-1])*v_old[i][j-1]+(y[i][j-1]-y1[i][j-1])*v_old[i-1][j-1])/(y[i][j-1]-y[i-1][j-1]);
								ve = ((x1[i][j+1]-x[i][j])*vf_old+(x[i][j]-x1[i][j])*vf_e)/(x1[i][j+1]-x1[i][j]);
								vw = ((x1[i][j]-x[i][j-1])*vf_w+(x[i][j-1]-x1[i][j-1])*vf_old)/(x1[i][j]-x1[i][j-1]);
								tempfc = tempfc + (((vf-vf_old)/dt)+1.5*(((u_old[i][j]*ve-u_old[i][j-1]*vw)/(x[i][j]-x[i][j-1]))+((v_old[i][j]*v_old[i][j]-v_old[i-1][j]*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
									}
							/*else
							{
								cout<<"wrong no.1LF"<<endl;
							}*/
						}
						else if (thetafc[k_f]!=api && (xf[k_f]-xf[k_f-1])<=0.0 && (Afc1*x1[i][j]+Bfc1*y1[i][j]+Cfc1)<=0.0)
						{
							if (thetafc[k_f+1]!=-api && (xf[k_f+1]-xf[k_f])>=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)<=0.0)
							{
								vf = ((y1[i][j]-y[i-1][j])*v[i][j]+(y[i][j]-y1[i][j])*v[i-1][j])/(y[i][j]-y[i-1][j]);
								vf_old = ((y1[i][j]-y[i-1][j])*v_old[i][j]+(y[i][j]-y1[i][j])*v_old[i-1][j])/(y[i][j]-y[i-1][j]);
								vf_e = ((y1[i][j+1]-y[i-1][j+1])*v_old[i][j+1]+(y[i][j+1]-y1[i][j+1])*v_old[i-1][j+1])/(y[i][j+1]-y[i-1][j+1]);
								vf_w = ((y1[i][j-1]-y[i-1][j-1])*v_old[i][j-1]+(y[i][j-1]-y1[i][j-1])*v_old[i-1][j-1])/(y[i][j-1]-y[i-1][j-1]);
								ve = ((x1[i][j+1]-x[i][j])*vf_old+(x[i][j]-x1[i][j])*vf_e)/(x1[i][j+1]-x1[i][j]);
								vw = ((x1[i][j]-x[i][j-1])*vf_w+(x[i][j-1]-x1[i][j-1])*vf_old)/(x1[i][j]-x1[i][j-1]);
								tempfc = tempfc + (((vf-vf_old)/dt)+1.5*(((u_old[i][j]*ve-u_old[i][j-1]*vw)/(x[i][j]-x[i][j-1]))+((v_old[i][j]*v_old[i][j]-v_old[i-1][j]*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
							}
							else if (thetafc[k_f+1]!=api && (xf[k_f+1]-xf[k_f])<=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)>=0.0)
							{
								vf = ((y1[i][j]-y[i-1][j])*v[i][j]+(y[i][j]-y1[i][j])*v[i-1][j])/(y[i][j]-y[i-1][j]);
								vf_old = ((y1[i][j]-y[i-1][j])*v_old[i][j]+(y[i][j]-y1[i][j])*v_old[i-1][j])/(y[i][j]-y[i-1][j]);
								vf_e = ((y1[i][j+1]-y[i-1][j+1])*v_old[i][j+1]+(y[i][j+1]-y1[i][j+1])*v_old[i-1][j+1])/(y[i][j+1]-y[i-1][j+1]);
								vf_w = ((y1[i][j-1]-y[i-1][j-1])*v_old[i][j-1]+(y[i][j-1]-y1[i][j-1])*v_old[i-1][j-1])/(y[i][j-1]-y[i-1][j-1]);
								ve = ((x1[i][j+1]-x[i][j])*vf_old+(x[i][j]-x1[i][j])*vf_e)/(x1[i][j+1]-x1[i][j]);
								vw = ((x1[i][j]-x[i][j-1])*vf_w+(x[i][j-1]-x1[i][j-1])*vf_old)/(x1[i][j]-x1[i][j-1]);
								tempfc = tempfc + (((vf-vf_old)/dt)+1.5*(((u_old[i][j]*ve-u_old[i][j-1]*vw)/(x[i][j]-x[i][j-1]))+((v_old[i][j]*v_old[i][j]-v_old[i-1][j]*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
							}
							/*else
							{
								cout<<"wrong no.2LF"<<endl;
							}*/
						}
					}
				}
			}
            tempf[k_f]=tempfc;
        }
        temp=0.0;
        //#pragma omp parallel for default(shared) private(vf,vf_old,vf_e,vf_w,ve,vw) reduction(+:temp) schedule(dynamic)
        for(int i=c33_f;i<c44_f;i++)
        {
            for (int j=c11_f;j<c22_f;j++)
            {
				if(flag[i][j]==0)
				{
					vf = ((y1[i][j]-y[i-1][j])*v[i][j]+(y[i][j]-y1[i][j])*v[i-1][j])/(y[i][j]-y[i-1][j]);
					vf_old = ((y1[i][j]-y[i-1][j])*v_old[i][j]+(y[i][j]-y1[i][j])*v_old[i-1][j])/(y[i][j]-y[i-1][j]);
					vf_e = ((y1[i][j+1]-y[i-1][j+1])*v_old[i][j+1]+(y[i][j+1]-y1[i][j+1])*v_old[i-1][j+1])/(y[i][j+1]-y[i-1][j+1]);
					vf_w = ((y1[i][j-1]-y[i-1][j-1])*v_old[i][j-1]+(y[i][j-1]-y1[i][j-1])*v_old[i-1][j-1])/(y[i][j-1]-y[i-1][j-1]);
					ve = ((x1[i][j+1]-x[i][j])*vf_old+(x[i][j]-x1[i][j])*vf_e)/(x1[i][j+1]-x1[i][j]);
					vw = ((x1[i][j]-x[i][j-1])*vf_w+(x[i][j-1]-x1[i][j-1])*vf_old)/(x1[i][j]-x1[i][j-1]);
					temp += (((vf-vf_old)/dt)+1.5*(((u_old[i][j]*ve-u_old[i][j-1]*vw)/(x[i][j]-x[i][j-1]))+((v_old[i][j]*v_old[i][j]-v_old[i-1][j]*v_old[i-1][j])/(y[i][j]-y[i-1][j]))))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
				}
            }
        }
		for(int k_f=0; k_f<Nf; k_f++)
        {
            temp1fc=0.0;
			minfx1=min(x_sfc[k_f],x_sfc[k_f+1]);
			minfx2=min(x_sfc[N_markerf-k_f-1],x_sfc[N_markerf-k_f-2]);
			maxfx1=max(x_sfc[k_f],x_sfc[k_f+1]);
			maxfx2=max(x_sfc[N_markerf-k_f-1],x_sfc[N_markerf-k_f-2]);
			minfy1=min(y_sfc[k_f],y_sfc[k_f+1]);
			minfy2=min(y_sfc[N_markerf-k_f-1],y_sfc[N_markerf-k_f-2]);
			maxfy1=max(y_sfc[k_f],y_sfc[k_f+1]);
			maxfy2=max(y_sfc[N_markerf-k_f-1],y_sfc[N_markerf-k_f-2]);
            cf1=m_div1+ceil((min(minfx1,minfx2)-Hx1)/dx2)-10;
			cf2=m_div1+ceil((max(maxfx1,maxfx2)-Hx1)/dx2)+10;
			cf3=n_div1+ceil((min(minfy1,minfy2)-Hy1)/dy2)-10;
			cf4=n_div1+ceil((max(maxfy1,maxfy2)-Hy1)/dy2)+10;

			if ((y_sfc[k_f]-y_sfc[N_markerf-k_f-1])>=0.0)
            {
                Afc1=(y_sfc[k_f]-y_sfc[N_markerf-k_f-1]);
                Afc2=(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2]);
                Bfc1=-(x_sfc[k_f]-x_sfc[N_markerf-k_f-1]);
                Bfc2=-(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2]);;
                Cfc1=((y_sfc[N_markerf-k_f-1]*(x_sfc[k_f]-x_sfc[N_markerf-k_f-1])) - (x_sfc[N_markerf-k_f-1]*(y_sfc[k_f]-y_sfc[N_markerf-k_f-1])));
                Cfc2=((y_sfc[N_markerf-k_f-2]*(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2])) - (x_sfc[N_markerf-k_f-2]*(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2])));
            }
            else
            {
                Afc1=-(y_sfc[k_f]-y_sfc[N_markerf-k_f-1]);
                Afc2=-(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2]);
                Bfc1=(x_sfc[k_f]-x_sfc[N_markerf-k_f-1]);
                Bfc2=(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2]);;
                Cfc1=-((y_sfc[N_markerf-k_f-1]*(x_sfc[k_f]-x_sfc[N_markerf-k_f-1])) - (x_sfc[N_markerf-k_f-1]*(y_sfc[k_f]-y_sfc[N_markerf-k_f-1])));
                Cfc2=-((y_sfc[N_markerf-k_f-2]*(x_sfc[k_f+1]-x_sfc[N_markerf-k_f-2])) - (x_sfc[N_markerf-k_f-2]*(y_sfc[k_f+1]-y_sfc[N_markerf-k_f-2])));
            }
			//cout<< cf1 << " " <<cf2 << " " << cf3 << " " << cf4 <<endl;
			for (int i=cf3;i<cf4;i++)
			{
				for (int j=cf1;j<cf2;j++)
				{
					if (flag_v[i][j] == 3)
					{
						if (thetafc[k_f]!=-api && (xf[k_f]-xf[k_f-1])>=0.0 && (Afc1*x1[i][j]+Bfc1*y1[i][j]+Cfc1)>=0.0)
						{
							if (thetafc[k_f+1]!=-api && (xf[k_f+1]-xf[k_f])>=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)<=0.0)
							{
								temp1fc = temp1fc + Bf_v[i][j]*(x[i][j]-x[i][j-1])*(y1[i+1][j]-y1[i][j]);
							}
							else if (thetafc[k_f+1]!=api && (xf[k_f+1]-xf[k_f])<=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)>=0.0)
							{
								temp1fc = temp1fc + Bf_v[i][j]*(x[i][j]-x[i][j-1])*(y1[i+1][j]-y1[i][j]);
							}
							/*else
							{
								cout<<"wrong no.1Lf"<<endl;
							}*/
						}
						else if (thetafc[k_f]!=api && (xf[k_f]-xf[k_f-1])<=0.0 && (Afc1*x1[i][j]+Bfc1*y1[i][j]+Cfc1)<=0.0)
						{
							if (thetafc[k_f+1]!=-api && (xf[k_f+1]-xf[k_f])>=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)<=0.0)
							{
								temp1fc = temp1fc + Bf_v[i][j]*(x[i][j]-x[i][j-1])*(y1[i+1][j]-y1[i][j]);
							}
							else if (thetafc[k_f+1]!=api && (xf[k_f+1]-xf[k_f])<=0.0 && (Afc2*x1[i][j]+Bfc2*y1[i][j]+Cfc2)>=0.0)
							{
								temp1fc = temp1fc + Bf_v[i][j]*(x[i][j]-x[i][j-1])*(y1[i+1][j]-y1[i][j]);
							}
							/*else
							{
								cout<<"wrong no.2Lf"<<endl;
							}*/
						}
					}
				}
			}
            temp1f[k_f]=temp1fc;
        }
        temp1 = 0.0;
        for(int i=c33_f;i<c44_f;i++)
        {
            for(int j=c11_f;j<c22_f;j++)
            {
                if (flag_v[i][j] == 3)
				{
                    temp1 = temp1 + Bf_v[i][j]*(x[i][j]-x[i][j-1])*(y1[i+1][j]-y1[i][j]);
				}
            }
        }
        for(int k_f=0; k_f<Nf; k_f++)
        {
            lift3f[k_f] = 2.0*(-temp1f[k_f]+tempf[k_f]+temp_Lf[k_f]);
        }
        lift3 = 2.0*(-temp1+temp+temp_L); //gives lift coeff Cl = L/(0.5*rho*u^2*(d*1))

        /*        temp=0.0;
         for(int i=0;i<n;i++)
         {
         for (int j=0;j<m;j++)
         {
         if(j>c1 && j<c2 && i>c3 && i<c4)
         {
         if(flag[i][j]==0)
         {
         uf = ((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u[i][j]+((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u[i][j-1];
         uf_n = ((x1[i+1][j]-x[i+1][j-1])/(x[i+1][j]-x[i+1][j-1]))*u[i+1][j]+((x[i+1][j]-x1[i+1][j])/(x[i+1][j]-x[i+1][j-1]))*u[i+1][j-1];
         uf_s = ((x1[i-1][j]-x[i-1][j-1])/(x[i-1][j]-x[i-1][j-1]))*u[i-1][j]+((x[i-1][j]-x1[i-1][j])/(x[i-1][j]-x[i-1][j-1]))*u[i-1][j-1];
         uf_old = ((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u_old[i][j]+((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u_old[i][j-1];
         un = ((y[i][j]-y1[i][j])/(y1[i+1][j]-y1[i][j]))*uf_n+((y1[i+1][j]-y[i][j])/(y1[i+1][j]-y1[i][j]))*uf;
         us = ((y[i-1][j]-y1[i-1][j])/(y1[i][j]-y1[i-1][j]))*uf+((y1[i][j]-y[i-1][j])/(y1[i][j]-y1[i-1][j]))*uf_s;
         temp = temp + (((uf-uf_old)/dt)+((u[i][j]*u[i][j]-u[i][j-1]*u[i][j-1])/(x[i][j]-x[i][j-1]))+((un*v[i][j]-us*v[i-1][j])/(y[i][j]-y[i-1][j])))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
         }
         }
         }
         }
         temp1 = 0.0;
         for(int i=0;i<n_u;i++)
         {
         for(int j=0;j<m_u;j++)
         {
         if (flag_u[i][j] == 3)
         temp1 = temp1 + Bf_u[i][j]*(x1[i][j+1]-x1[i][j])*(y[i][j]-y[i-1][j]);
         }
         }
         drag1 = 2*(-temp1+temp); //gives drag coeff Cd = D/(0.5*rho*u^2*(d*1))
         temp=0.0;
         for(int i=0;i<n;i++)
         {
         for (int j=0;j<m;j++)
         {
         if(j>c1 && j<c2 && i>c3 && i<c4)
         {
         if(flag[i][j]==0)
         {
         vf = ((y1[i][j]-y[i-1][j])/(y[i][j]-y[i-1][j]))*v[i][j]+((y[i][j]-y1[i][j])/(y[i][j]-y[i-1][j]))*v[i-1][j];
         vf_e = ((y1[i][j+1]-y[i-1][j+1])/(y[i][j+1]-y[i-1][j+1]))*v[i][j+1]+((y[i][j+1]-y1[i][j+1])/(y[i][j+1]-y[i-1][j+1]))*v[i-1][j+1];
         vf_w = ((y1[i][j-1]-y[i-1][j-1])/(y[i][j-1]-y[i-1][j-1]))*v[i][j-1]+((y[i][j-1]-y1[i][j-1])/(y[i][j-1]-y[i-1][j-1]))*v[i-1][j-1];
         vf_old = ((y1[i][j]-y[i-1][j])/(y[i][j]-y[i-1][j]))*v_old[i][j]+((y[i][j]-y1[i][j])/(y[i][j]-y[i-1][j]))*v_old[i-1][j];
         ve = ((x1[i][j+1]-x[i][j])/(x1[i][j+1]-x1[i][j]))*vf+((x[i][j]-x1[i][j])/(x1[i][j+1]-x1[i][j]))*vf_e;
         vw = ((x1[i][j]-x[i][j-1])/(x1[i][j]-x1[i][j-1]))*vf_w+((x[i][j-1]-x1[i][j-1])/(x1[i][j]-x1[i][j-1]))*vf;
         temp = temp + (((vf-vf_old)/dt)+((u[i][j]*ve-u[i][j-1]*vw)/(x[i][j]-x[i][j-1]))+((v[i][j]*v[i][j]-v[i-1][j]*v[i-1][j])/(y[i][j]-y[i-1][j])))*((x[i][j]-x[i][j-1])*(y[i][j]-y[i-1][j]));
         }
         }
         }
         }
         temp1 = 0.0;
         for(int i=0;i<n_v;i++)
         {
         for(int j=0;j<m_v;j++)
         {
         if (flag_v[i][j] == 3)
         temp1 = temp1 + Bf_v[i][j]*(x[i][j]-x[i][j-1])*(y1[i+1][j]-y1[i][j]);
         }
         }
         lift1 = 2*(-temp1+temp); //gives lift coeff Cl = L/(0.5*rho*u^2*(d*1))

         if (file5a.is_open())
         {
         file5a << time1 << " " << drag1 << " " << lift1 << endl;
         //file5a.close();
         }*/
        if (file5c.is_open())
        {
            file5c << fixed << setprecision(6) << time1 << " " << drag3 << " " << lift3 << endl;
            //file5c.close();
        }
        /*if (file6.is_open())
        {
            file6 << residual_flux << endl;
            //file6.close();
        }*/
        if (file6a.is_open())
        {
            file6a << time1 << " " << ((x1[probe_i][probe_j]-x[probe_i][probe_j-1])/(x[probe_i][probe_j]-x[probe_i][probe_j-1]))*u[probe_i][probe_j]+((x[probe_i][probe_j]-x1[probe_i][probe_j])/(x[probe_i][probe_j]-x[probe_i][probe_j-1]))*u[probe_i][probe_j-1] << " " << ((y1[probe_i][probe_j]-y[probe_i-1][probe_j])/(y[probe_i][probe_j]-y[probe_i-1][probe_j]))*v[probe_i][probe_j]+((y[probe_i][probe_j]-y1[probe_i][probe_j])/(y[probe_i][probe_j]-y[probe_i-1][probe_j]))*v[probe_i-1][probe_j] << endl;
            //file6a.close();
        }
        cout << k << endl;
        if ( (k%writeInterval) == 0)
        {
            /*ofstream filew1("Cavity_u.txt");
             if (filew1.is_open())
             {
             for(int i=0; i<n_u; i++)
             {
             for(int j=0; j<m_u; j++)
             {
             filew1 << u[i][j] << " ";
             }
             filew1 << endl;
             }
             filew1.close();
             }

             ofstream filew2("Cavity_v.txt");
             if (filew2.is_open())
             {
             for(int i=0; i<n_v; i++)
             {
             for(int j=0; j<m_v; j++)
             {
             filew2 << v[i][j] << " ";
             }
             filew2 << endl;
             }
             filew2.close();
             }
             ofstream filew3("pressure.txt");
             if (filew3.is_open())
             {
             for(int i=0; i<n; i++)
             {
             for(int j=0; j<m; j++)
             {
             filew3 << pressure[i][j] << " ";
             }
             filew3 << endl;
             }
             filew3.close();
             }
             */
             string fileName;
            stringstream ss;
            ss << k;
            fileName = ss.str();
            fileName = "Cavity_u"+fileName+".txt";
            ofstream filew1(fileName.c_str());
            if (filew1.is_open())
            {
                for(int i=0; i<n_u; i++)
                {
                    for(int j=0; j<m_u; j++)
                    {
                        filew1 << u[i][j] << endl;
                    }
                }
                filew1.close();
            }
            fileName = ss.str();
            fileName = "Cavity_v"+fileName+".txt";
            ofstream filew2(fileName.c_str());
            if (filew2.is_open())
            {
                for(int i=0; i<n_v; i++)
                {
                    for(int j=0; j<m_v; j++)
                    {
                        filew2 << v[i][j] << endl;
                    }
                }
                filew2.close();
            }
            fileName = ss.str();
            fileName = "Cavity_u_old"+fileName+".txt";
            ofstream filew1a(fileName.c_str());
            if (filew1a.is_open())
            {
                for(int i=0; i<n_u; i++)
                {
                    for(int j=0; j<m_u; j++)
                    {
                        filew1a << u_old[i][j] << endl;
                    }
                }
                filew1.close();
            }
            fileName = ss.str();
            fileName = "Cavity_v_old"+fileName+".txt";
            ofstream filew2a(fileName.c_str());
            if (filew2a.is_open())
            {
                for(int i=0; i<n_v; i++)
                {
                    for(int j=0; j<m_v; j++)
                    {
                        filew2a << v_old[i][j] << endl;
                    }
                }
                filew2a.close();
            }
            fileName = ss.str();
            fileName = "pressure"+fileName+".txt";
            ofstream filew3(fileName.c_str());
            if (filew3.is_open())
            {
                for(int i=0; i<n; i++)
                {
                    for(int j=0; j<m; j++)
                    {
                        filew3 << pressure[i][j] << endl;
                    }
                }
                filew3.close();
            }
            fileName = ss.str();
            fileName = "p_c"+fileName+".txt";
            ofstream filew3a(fileName.c_str());
            if (filew3a.is_open())
            {
                for(int i=0; i<n; i++)
                {
                    for(int j=0; j<m; j++)
                    {
                        filew3a << p_c[i][j] << endl;
                    }
                }
                filew3a.close();
            }
            fileName = ss.str();
            fileName = "Cavity_u_red"+fileName+".txt";
            ofstream filew1f(fileName.c_str());
            if (filew1f.is_open())
            {
                for(int i=0; i<n_u; i++)
                {
                    for(int j=0; j<m_u; j++)
                    {
                        filew1f << u_red[i][j] << endl;
                    }
                }
                filew1f.close();
            }
            fileName = ss.str();
            fileName = "Cavity_v_red"+fileName+".txt";
            ofstream filew2f(fileName.c_str());
            if (filew2f.is_open())
            {
                for(int i=0; i<n_v; i++)
                {
                    for(int j=0; j<m_v; j++)
                    {
                        filew2f << v_red[i][j] << endl;
                    }
                }
                filew2f.close();
            }
            fileName = ss.str();
            fileName = "Cavity_u_black"+fileName+".txt";
            ofstream filew1fa(fileName.c_str());
            if (filew1fa.is_open())
            {
                for(int i=0; i<n_u; i++)
                {
                    for(int j=0; j<m_u; j++)
                    {
                        filew1fa << u_black[i][j] << endl;
                    }
                }
                filew1fa.close();
            }
            fileName = ss.str();
            fileName = "Cavity_v_black"+fileName+".txt";
            ofstream filew2fa(fileName.c_str());
            if (filew2fa.is_open())
            {
                for(int i=0; i<n_v; i++)
                {
                    for(int j=0; j<m_v; j++)
                    {
                        filew2fa << v_black[i][j] << endl;
                    }
                }
                filew2fa.close();
            }
            fileName = ss.str();
            fileName = "p_c_red"+fileName+".txt";
            ofstream filew3f(fileName.c_str());
            if (filew3f.is_open())
            {
                for(int i=0; i<n; i++)
                {
                    for(int j=0; j<m; j++)
                    {
                        filew3f << p_c_red[i][j] << endl;
                    }
                }
                filew3f.close();
            }
            fileName = ss.str();
            fileName = "p_c_black"+fileName+".txt";
            ofstream filew3fa(fileName.c_str());
            if (filew3fa.is_open())
            {
                for(int i=0; i<n; i++)
                {
                    for(int j=0; j<m; j++)
                    {
                        filew3fa << p_c_black[i][j] << endl;
                    }
                }
                filew3fa.close();
            }
            fileName = ss.str();
            fileName = "lift3f"+fileName+".txt";
            ofstream filew3lf(fileName.c_str());
            if (filew3lf.is_open())
            {
                for(int i=0; i<Nf; i++)
                {
                     filew3lf << lift3f[i] << endl;
                }
                filew3lf.close();
            }
            fileName = ss.str();
            fileName = "drag3f"+fileName+".txt";
            ofstream filew3df(fileName.c_str());
            if (filew3df.is_open())
            {
                for(int i=0; i<Nf; i++)
                {
                     filew3df << drag3f[i] << endl;
                }
                filew3df.close();
            }
            for (int i=0;i<n;i++)
            {
                for (int j=0;j<m;j++)
                {
                    if (i==0 && j==0)//bottom row
                    {
                        uf1[i][j] = u[i][j];
                        vf1[i][j] = v[i][j];
                        pressure[i][j] = 2*pressure[i+1][j]-pressure[i+2][j];
                    }
                    else if (i==0 && j>0 && j <m-1)
                    {
                        uf1[i][j] = ((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u[i][j]+((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u[i][j-1];
                        vf1[i][j] = v[i][j];
                        pressure[i][j] = 2*pressure[i+1][j]-pressure[i+2][j];
                    }
                    else if (i==0 && j==m-1)
                    {
                        uf1[i][j] = u[i][j-1];
                        vf1[i][j] = v[i][j];
                        pressure[i][j] = 2*pressure[i+1][j]-pressure[i+2][j];
                    }
                    else if (i>0 && i<n-1 && j==0)//left column
                    {
                        uf1[i][j] = u[i][j];
                        vf1[i][j] = ((y1[i][0]-y[i-1][0])/(y[i][0]-y[i-1][0]))*v[i][j]+((y[i][0]-y1[i][0])/(y[i][0]-y[i-1][0]))*v[i-1][j];
                        pressure[i][j] = 2*pressure[i][j+1]-pressure[i][j+2];
                    }
                    else if (i>0 && i<n-1 && j==m-1)//right column
                    {
                        uf1[i][j] = u[i][j-1];
                        vf1[i][j] = ((y1[i][0]-y[i-1][0])/(y[i][0]-y[i-1][0]))*v[i][j]+((y[i][0]-y1[i][0])/(y[i][0]-y[i-1][0]))*v[i-1][j];
                        pressure[i][j] = 2*pressure[i][j-1]-pressure[i][j-2];
                    }
                    else if (i==n-1 && j==0)
                    {
                        uf1[i][j] = u[i][j];
                        vf1[i][j] = v[i-1][j];
                        pressure[i][j] = 2*pressure[i-1][j]-pressure[i-2][j];
                    }
                    else if (i==n-1 && j>0 && j <m-1)//top row
                    {
                        uf1[i][j] = ((x1[0][j]-x[0][j-1])/(x[0][j]-x[0][j-1]))*u[i][j]+((x[0][j]-x1[0][j])/(x[0][j]-x[0][j-1]))*u[i][j-1];
                        vf1[i][j] = v[i-1][j];
                        pressure[i][j] = 2*pressure[i-1][j]-pressure[i-2][j];
                    }
                    else if (i==n-1 && j==m-1)
                    {
                        uf1[i][j] = u[i][j-1];
                        vf1[i][j] = v[i-1][j];
                        pressure[i][j] = 2*pressure[i-1][j]-pressure[i-2][j];
                    }
                    else
                    {
                        uf1[i][j] = ((x1[i][j]-x[i][j-1])/(x[i][j]-x[i][j-1]))*u[i][j]+((x[i][j]-x1[i][j])/(x[i][j]-x[i][j-1]))*u[i][j-1];
                        vf1[i][j] = ((y1[i][j]-y[i-1][j])/(y[i][j]-y[i-1][j]))*v[i][j]+((y[i][j]-y1[i][j])/(y[i][j]-y[i-1][j]))*v[i-1][j];
                        pressure[i][j] = pressure[i][j];
                    }
                }
            }
            fileName = ss.str();
            fileName = "AllinOne"+fileName+".dat";
            ofstream file7(fileName.c_str());
            if (file7.is_open())
            {
                file7 << "ZONE T=DATA I=" << m << " J=" << n << endl;
                for(int i=0; i<n; i++)
                {
                    for(int j=0; j<m; j++)
                    {
                        file7 << x1[i][j] << " " << y1[i][j] << " " << uf1[i][j] << " " << vf1[i][j] << " " << pressure[i][j] << endl;
                    }
                }
                file7.close();
            }
            /*ofstream file8("error_u.txt");
             if (file8.is_open())
             {
             for(int i=0; i<GSite; i++)
             {
             file8 << max_error1[i] << endl;
             }
             file8.close();
             }
             ofstream file9("error_v.txt");
             if (file9.is_open())
             {
             for(int i=0; i<GSite; i++)
             {
             file9 << max_error2[i] << endl;
             }
             file9.close();
             }
             ofstream file10("error_p_c.txt");
             if (file10.is_open())
             {
             for(int i=0; i<GSite2; i++)
             {
             file10 << max_error3[i] << endl;
             }
             file10.close();
             }*/
        }
    }
    /*ofstream file8("error_u.txt");
     if (file8.is_open())
     {
     for(int i=0; i<GSite; i++)
     {
     file8 << max_error1[i] << endl;
     }
     file8.close();
     }
     ofstream file9("error_v.txt");
     if (file9.is_open())
     {
     for(int i=0; i<GSite; i++)
     {
     file9 << max_error2[i] << endl;
     }
     file9.close();
     }
     ofstream file10("error_p_c.txt");
     if (file10.is_open())
     {
     for(int i=0; i<GSite2; i++)
     {
     file10 << max_error3[i] << endl;
     }
     file10.close();
     }*/
    /*for (int ctr=0; ctr<N_marker; ctr++)
     {
     min_dis=10;
     for (int i=1;i<n-1;i++)
     {
     for (int j=1;j<m-1;j++)
     {
     if(j>c1 && j<c2 && i>c3 && i<c4)
     {
     if(flag[i][j]==1)
     {
     distance = sqrt(pow((x_s[ctr]-x1[i][j]),2)+pow((y_s[ctr]-y1[i][j]),2));
     if (distance<min_dis)
     {
     min_dis = distance;
     i_min = i;
     j_min = j;
     }
     }
     }
     }
     }
     surface_pr[ctr] = pressure[i_min][j_min];
     }
     ofstream file11("surface_pressure.txt");
     if (file11.is_open())
     {
     for(int i=0; i<N_marker; i++)
     {
     file11 << 2*surface_pr[i] << endl;
     }
     file11 << 2*surface_pr[0] << endl;
     file11.close();
     }*/
    cout << "m = " << m << "   n = " << n << endl;
    cout << "m_u = " << m_u << "   n_u = " << n_u << "   m_v = " << m_v << "   n_v = " << n_v << endl;
    cout << "dx1 = " << dx1 << "   dx2 = " << dx2 << endl;
    cout << "dy1 = " << dy1 << "   dy2 = " << dy2 << endl;
    //cout << "iterationN = " << iterationN << endl;
    return 0;
}
