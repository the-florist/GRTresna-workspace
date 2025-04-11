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

        void set_random_scalar_field(FArrayBox &a_multigrid_vars_box, Box &a_ghosted_box)
        {
            // Set the random generator
            std::mt19937 generator;
            std::uniform_real_distribution<double> distribution(0.0, 1.0);

            FArrayBox random_numbers(a_ghosted_box, 2);
            for (BoxIterator bit(a_ghosted_box); bit.ok(); ++bit)
            {
                IntVect iv = bit();
                int offset = a_multigrid_vars_box.minIndex() * 2;

                /*std::cout << iv << ": ";
                std::cout << a_multigrid_vars_box.min() << ", ";
                std::cout << a_multigrid_vars_box.minIndex() << "\n";*/
            }

            MayDay::Error("RandomScalarField::set_random_scalar_field: not implemented yet");

            BoxIterator bit(a_ghosted_box);
            for (bit.begin(); bit.ok(); ++bit)
            {

                // work out location on the grid
                IntVect iv = bit();
                a_multigrid_vars_box(iv, c_phi_0) = 0.;//random_fields(iv, 0);
                a_multigrid_vars_box(iv, c_Pi_0) = 0.;//random_fields(iv, 1);
            }
        }

};