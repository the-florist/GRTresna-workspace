/* GRTresna
 * Copyright 2022 The GRTL collaboration.
 */

 #ifndef RANDOMSCALARFIELD_HPP_
 #define RANDOMSCALARFIELD_HPP_

 #include "fftw3.h"
 #include <random>


class RandomScalarField
{
    public:
        RandomScalarField(double mass) : m_mass(mass)
        {}

        void set_random_scalar_field(FArrayBox &a_multigrid_vars_box, Box &a_ghosted_box, 
            const RealVect &a_dx, const std::array<double, SpaceDim> center)
        {
            // Set the random generator
            //std::mt19937 generator;
            //std::uniform_real_distribution<double> distribution(0.0, 1.0);
	    default_random_engine generator(142564253);
	    uniform_real_distribution<double> distribution(0, 1);

            int N = 64;
            double L = 3000.;

	    pout() << "Starting random draws.\n";

            double** random_draws;
            random_draws = (double**) malloc(sizeof(double*) * 4);
            for(int l=0; l<4; l++)
            {
                random_draws[l] = (double*) malloc(sizeof(double) * N * N * N);
            }

            //FArrayBox random_numbers(a_ghosted_box, 2);
            for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<=N/2; k++)
            {
                int offset = k + N*(j + N*i);
                for(int ofs=0; ofs<offset; ofs++) { distribution(generator); }
                for(int l=0; l<4; l++) { random_draws[l][offset] = distribution(generator); }
            }

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

	    std::cout << "Random draws complete and all memory allocated.\n";
	    std::cout << "Starting main grid loop.\n";

            for(int i=0; i<N; i++) for(int j=0; j<N; j++) for(int k=0; k<=N/2; k++)
            {
                // Put the 0 mode in the right spot in memory
                int offset = k + (N/2+1)*(j + N*i);

                int I = (i < N/2) ? i : N-i;
                int J = (j < N/2) ? j : N-j;

                // Find |k|
                double kmag = (double)(pow((pow((double)I, 2.0) + pow((double)J, 2.0) + pow((double)k, 2.0))*4.*M_PI*M_PI/L/L, 0.5));
                for(int s=0; s<2; s++) 
                {
                    phi_k[0][offset][s] = std::sqrt(-2. * log(random_draws[s][offset])/(2. * std::sqrt(std::pow(kmag, 2.) + std::pow(m_mass, 2.))));
                    phi_k[1][offset][s] = std::sqrt(-2. * log(random_draws[s][offset]) *  std::sqrt(std::pow(kmag, 2.) + std::pow(m_mass, 2.))/2.);
                }
                phi_k[0][offset][0] *= cos(2. * M_PI * random_draws[2][offset]);
                phi_k[0][offset][1] *= sin(2. * M_PI * random_draws[3][offset]);

                if ((i == 0 || i == N/2) && (j == 0 || j == N/2) && (k == 0 || k == N/2))
                {
                    for(int s=0; s<2; s++) { phi_k[s][k + (N/2+1)*(j + N*i)][1] = 0.; }
                }
                if (k==0 || k==N/2) 
                {
                    if((i>N/2 && j==N/2) || (i==0 && j>N/2) || (i>N/2 && j==0) || (i==N/2 && j>N/2)) // Special lines
                    {
                        for(int s=0; s<2; s++) 
                        {
                            phi_k[s][k + (N/2+1)*(j + N*i)][0] = phi_k[s][k + (N/2+1)*(J + N*I)][0];
                            phi_k[s][k + (N/2+1)*(j + N*i)][1] = -phi_k[s][k + (N/2+1)*(J + N*I)][1];
                        }
                    }
                    else if(j > N/2) // Special plane bulk
                    {
                        for(int s=0; s<2; s++) 
                        {
                            phi_k[s][k + (N/2+1)*(j + N*i)][0] = phi_k[s][k + (N/2+1)*(J + N*std::abs(N-i))][0];
                            phi_k[s][k + (N/2+1)*(j + N*i)][1] = -phi_k[s][k + (N/2+1)*(J + N*std::abs(N-i))][1];
                        }
                    }
                }
            }

	    std::cout << "Main grid loop complete, starting fourier transform.\n";

            for(int s=0; s<2; s++) { fftw_execute(phi_plan[s]); }

            //MayDay::Error("RandomScalarField::set_random_scalar_field: not implemented yet");

	    std::cout << "FT complete, assigning grid values.\n";

            BoxIterator bit(a_ghosted_box);
            for (bit.begin(); bit.ok(); ++bit)
            {
                // work out location on the grid
                IntVect iv = bit();
		
		int i = iv[0];
		int j = iv[1];
		int k = iv[2];;
		if(i < 0) { i += N; }
		else if(i >= N) { i -= N; }
		if(j < 0) { j += N; }
		else if(j >= N) { j -= N; }
		if(k < 0) { k += N; }
		else if(k >= N) { k -= N; }	
                
		int offset = k + N*(j + N*i);
		//std::cout << iv << "\n";
                a_multigrid_vars_box(iv, c_phi_0) = phi_x[0][offset];
                a_multigrid_vars_box(iv, c_Pi_0) = phi_x[1][offset];
            }
        }

    protected:
        double m_mass;

};

#endif /* RANDOMSCALARFIELD_HPP_ */
