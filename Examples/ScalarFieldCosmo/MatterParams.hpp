/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef MATTERPARAMS_HPP_
#define MATTERPARAMS_HPP_

#include "GRParmParse.hpp"
#include "REAL.H"

namespace MatterParams
{

struct params_t
{
    Real phi_0;
    Real dphi;
    Real pi_0;
    Real dpi;
    Real scalar_mass;
    Real G_Newton;

    Real phi_s;
    Real phi_c;
    Real mu;
    Real V0;

    int use_random_field = 0;
    Real domain_length;
};

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    pp.get("phi_0", matter_params.phi_0);
    pp.load("dphi", matter_params.dphi);
    pp.load("pi_0", matter_params.pi_0);
    pp.load("dpi", matter_params.dpi);
    pp.get("scalar_mass", matter_params.scalar_mass);
    pp.get("G_Newton", matter_params.G_Newton);

    pp.load("phi_s", matter_params.phi_s);
    pp.load("phi_c", matter_params.phi_c);
    pp.load("mu", matter_params.mu);
    pp.load("V0", matter_params.V0);

    pp.get("use_random_field", matter_params.use_random_field);
    pp.get("L", matter_params.domain_length);
}

}; // namespace MatterParams

#endif