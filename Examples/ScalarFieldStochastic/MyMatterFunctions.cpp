/* GRTresna
 * Copyright 2024 The GRTL collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "ScalarField.hpp"

Real ScalarField::my_potential_function(const Real &phi_here) const
{
    Real V;

    if(phi_here < 0) { V = m_matter_params.V0; }
    else if (phi_here <= m_matter_params.phi_c) 
    { 
        V = m_matter_params.V0 * (1. - std::pow(phi_here/m_matter_params.mu, 4.)); 
    }
    else 
    { 
        V = std::pow(m_matter_params.scalar_mass * (phi_here - m_matter_params.phi_s), 2.); 
    }

    return V;
    //return 0.5 * pow(m_matter_params.scalar_mass * phi_here, 2.0);
}

Real ScalarField::my_potential_deriv1(const Real &phi_here) const
{
    Real dV;

    if(phi_here < 0) { dV = 0.; }
    else if (phi_here <= m_matter_params.phi_c) 
    { 
        dV = -4. * m_matter_params.V0 * std::pow(phi_here/m_matter_params.mu, 3.) / m_matter_params.mu; 
    }
    else 
    { 
        dV = 2. * phi_here * std::pow(m_matter_params.scalar_mass, 2.) 
            - 2. * m_matter_params.phi_s * std::pow(m_matter_params.scalar_mass, 2.); 
    }

    return dV;
}

Real ScalarField::my_potential_deriv2(const Real &phi_here) const
{
    Real ddV;

    if(phi_here < 0) { ddV = 0.; }
    else if (phi_here <= m_matter_params.phi_c) 
    { 
        ddV = -12. * m_matter_params.V0 * std::pow(phi_here/m_matter_params.mu, 2.) / std::pow(m_matter_params.mu, 2.); 
    }
    else 
    { 
        ddV = 2. * std::pow(m_matter_params.scalar_mass, 2.); 
    }

    return ddV;
}

Real ScalarField::my_phi_function(const RealVect &loc) const
{
    Real rr = sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
    Real L = domainLength[0];
    Real dphi_value = m_matter_params.dphi / 3. *
                      (sin(2 * M_PI * loc[0] / L) + sin(2 * M_PI * loc[1] / L) +
                       sin(2 * M_PI * loc[2] / L));
    return m_matter_params.phi_0 + dphi_value;
}

Real ScalarField::my_Pi_function(const RealVect &loc) const
{
    Real rr = sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
    Real L = domainLength[0];
    Real dpi_value = m_matter_params.dpi / 3. *
                     (sin(2 * M_PI * loc[0] / L) + sin(2 * M_PI * loc[1] / L) +
                      sin(2 * M_PI * loc[2] / L));
    return m_matter_params.pi_0 + dpi_value;
}
