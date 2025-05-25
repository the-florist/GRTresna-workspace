/* GRTresna
 * Copyright 2022 The GRTL collaboration.
 */

 #ifndef RANDOMSCALARFIELD_HPP_
 #define RANDOMSCALARFIELD_HPP_

 #include "MatterParams.hpp"
 #include "Grids.hpp"
 #include "fftw3.h"
 #include <random>


class RandomScalarField
{
    public:
        RandomScalarField(MatterParams::params_t a_matter_params, const Real a_V0, const RealVect &a_domainLength, const RealVect &a_dx) 
        : m_matter_params(a_matter_params), V0(a_V0), m_domainLength(a_domainLength), m_dx(a_dx)
        {
            // Error trap (all Ns must be equal for FFTW to work)
            for(int l=0; l<2; l++)
            { 
                if(m_domainLength[l] != m_domainLength[l+1] || m_dx[l] != m_dx[l+1])
                {
                    MayDay::Error("Grid resolution (N) must be cubic when RandomScalarField initial data is used");
                }
            }

            H0 = std::sqrt(8. * M_PI * m_matter_params.G_Newton * (0.5 * std::pow(m_matter_params.pi_0, 2.) + V0) / 3.);
        }

        void set_random_scalar_field(FArrayBox &a_multigrid_vars_box, Box &a_ghosted_box, const Real meff)
        {   
            // Set domain parameters
            Real L = m_domainLength[0];
            int N = L/m_dx[0];

            // Set the random generator
	        default_random_engine generator(142564253);
	        uniform_real_distribution<double> distribution(0, 1);

            // Allocate memory to store random draws
            // (which makes the random draws MPI safe)
            double** random_draws;
            random_draws = (double**) malloc(sizeof(double*) * 4);
            for(int l=0; l<4; l++)
            {
                random_draws[l] = (double*) malloc(sizeof(double) * N * N * N);
            }

            // Create random draws (independent of the MPI box)
            for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<=N/2; k++)
            {
                int offset = k + N*(j + N*i);
                for(int ofs=0; ofs<offset; ofs++) { distribution(generator); }
                for(int l=0; l<4; l++) { random_draws[l][offset] = distribution(generator); }
            }

            // Allocate memory to store the field and field derivative values
            fftw_complex** phi_k;
            double** phi_x;
            phi_k = (fftw_complex**) fftw_malloc(sizeof(fftw_complex*) * 2);
            phi_x = (double**) malloc(sizeof(double*) * 2);
            fftw_plan phi_plan[2];

            for(int l=0; l<2; l++)
            {
                phi_k[l] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N * (N/2+1));
                phi_x[l] = (double*) malloc(sizeof(double) * N * N * N);
                phi_plan[l] = fftw_plan_dft_c2r_3d(N, N, N, phi_k[l], phi_x[l], FFTW_ESTIMATE);
            }

            for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<=N/2; k++)
            {
                int offset = k + (N/2+1)*(j + N*i);

                // Put the 0 mode in the right spot in memory
                int I = (i < N/2) ? i : N-i;
                int J = (j < N/2) ? j : N-j;

                // Find |k|
                double kmag = (double)(pow((pow((double)I, 2.0) + pow((double)J, 2.0) + pow((double)k, 2.0))*4.*M_PI*M_PI/L/L, 0.5));

                for(int s=0; s<2; s++) 
                {
                    // Set the modulus with a Rayleigh draw
                    phi_k[0][offset][s] = std::sqrt(-2. * log(random_draws[s][offset]) 
                                            / (2. * std::sqrt(std::pow(kmag, 2.) + std::pow(meff, 2.))));
                    phi_k[1][offset][s] = std::sqrt(-2. * log(random_draws[s][offset]) 
                                            *  std::sqrt(std::pow(kmag, 2.) + std::pow(meff, 2.))/2.);

                    // Set the "physical" phase
                    double kpr = kmag/H0;
                    std::vector<double> ms_args{atan2((cos(kpr) + kpr*sin(kpr)), (kpr*cos(kpr) - sin(kpr))), -atan2(cos(kpr), sin(kpr))};   

                    for (int l=0; l<2; l++)
                    {
                        if(s==0) { phi_k[l][offset][s] *= cos(ms_args[l]); }
                        else if(s==1) { phi_k[l][offset][s] *= sin(ms_args[l]); }
                    }

                    // Set the stochastic phase
                    double rand_arg = 2. * M_PI * random_draws[2+s][offset];
                    for (int l=0; l<2; l++)
                    {
                        if(s==0) 
                        { 
                            double temp = phi_k[l][offset][s] * cos(rand_arg) - phi_k[l][offset][s-1] * sin(rand_arg);
                            phi_k[l][offset][s] = temp;
                        }
                        else if(s==1) 
                        { 
                            double temp = phi_k[l][offset][s-1] * sin(rand_arg) + phi_k[l][offset][s] * cos(rand_arg);
                            phi_k[l][offset][s] = temp;
                        }
                    }
                }

                // Create random draw (Bunch-Davies approximation, no MS phase)
                /*for(int s=0; s<2; s++) 
                {
                    phi_k[0][offset][s] = std::sqrt(-2. * log(random_draws[s][offset]) 
                                            / (2. * std::sqrt(std::pow(kmag, 2.) + std::pow(m_matter_params.scalar_mass, 2.))));
                    phi_k[1][offset][s] = std::sqrt(-2. * log(random_draws[s][offset]) 
                                            *  std::sqrt(std::pow(kmag, 2.) + std::pow(m_matter_params.scalar_mass, 2.))/2.);
                }
                phi_k[0][offset][0] *= cos(2. * M_PI * random_draws[2][offset]);
                phi_k[0][offset][1] *= sin(2. * M_PI * random_draws[3][offset]);

                //FIXME I think the velocity needs a different phase but idk what this is yet.
                phi_k[1][offset][0] *= cos(2. * M_PI * random_draws[2][offset]);
                phi_k[1][offset][1] *= sin(2. * M_PI * random_draws[3][offset]);*/
            }
            
            free(random_draws);

            // Apply Hermitian symmetry conditions
            for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<=N/2; k++)
            {
                int offset = k + (N/2+1)*(j + N*i);

                // Put the 0 mode in the right spot in memory
                int I = (i < N/2) ? i : N-i;
                int J = (j < N/2) ? j : N-j;

                // Nyquist points
                if ((i == 0 || i == N/2) && (j == 0 || j == N/2) && (k == 0 || k == N/2))
                {
                    for(int s=0; s<2; s++) { phi_k[s][k + (N/2+1)*(j + N*i)][1] = 0.; }
                }
                if (k==0 || k==N/2) 
                {
                    if((i>N/2 && j==N/2) || (i==0 && j>N/2) || (i>N/2 && j==0) || (i==N/2 && j>N/2)) // Nyquist lines
                    {
                        for(int s=0; s<2; s++) 
                        {
                            phi_k[s][k + (N/2+1)*(j + N*i)][0] = phi_k[s][k + (N/2+1)*(J + N*I)][0];
                            phi_k[s][k + (N/2+1)*(j + N*i)][1] = -phi_k[s][k + (N/2+1)*(J + N*I)][1];
                        }
                    }
                    else if(j > N/2) // Nyquist planes
                    {
                        for(int s=0; s<2; s++) 
                        {
                            phi_k[s][k + (N/2+1)*(j + N*i)][0] = phi_k[s][k + (N/2+1)*(J + N*std::abs(N-i))][0];
                            phi_k[s][k + (N/2+1)*(j + N*i)][1] = -phi_k[s][k + (N/2+1)*(J + N*std::abs(N-i))][1];
                        }
                    }
                }
            }

            for(int s=0; s<2; s++) 
            { 
                fftw_execute(phi_plan[s]); 
                fftw_destroy_plan(phi_plan[s]); 
            }

            BoxIterator bit(a_ghosted_box);
            for (bit.begin(); bit.ok(); ++bit)
            {
                // work out location on the grid
                IntVect iv = bit();
		
                int i = iv[0];
                int j = iv[1];
                int k = iv[2];

                if(i < 0) { i += N; }
                else if(i >= N) { i -= N; }
                if(j < 0) { j += N; }
                else if(j >= N) { j -= N; }
                if(k < 0) { k += N; }
                else if(k >= N) { k -= N; }	
                
		        int offset = k + N*(j + N*i);
                a_multigrid_vars_box(iv, c_phi_0) = phi_x[0][offset];
                a_multigrid_vars_box(iv, c_Pi_0) = phi_x[1][offset];
            }

            free(phi_x);
        }

    protected:
        MatterParams::params_t m_matter_params;
        const Real V0;
        const RealVect m_domainLength;
        const RealVect m_dx;
        Real H0;

};

#endif /* RANDOMSCALARFIELD_HPP_ */
