/*
 * Author: Nathaniel Rupprecht
 * Start Data: May 19, 2017
 *
 */

#ifndef __PARSING_TOKENS_HPP__
#define __PARSING_TOKENS_HPP__

#include <string>
using std::string;

const char Comment_Tok   = '/';
const string Variable_Tok= "let";
const string Bool_True   = "true";
const string Bool_False  = "false";

const string WrapX_Tok   = "wrapX";
const string WrapY_Tok   = "wrapY";
const string Bounds_Tok  = "bounds";
const string Gravity_Tok = "gravity";

const string Number_Tok  = "number";
const string Phi_Tok     = "phi";

const string Position_Tok = "pos";

const string Sigma_Tok   = "sigma";
const string Dispersion_Tok = "dispersion";

const string Normal_Velocity_Tok = "normal_velocity";
const string Normal_KE_Tok = "normal_ke";
const string Velocity_Tok = "velocity";

const string Omega_Tok = "omega";

const string Repulsion_Tok = "repulsion";
const string Dissipation_Tok = "dissipation";
const string Coeff_Tok = "coeff";
const string Interaction_Tok = "interaction";

const string Region_Tok  = "region";
const string Wall_Tok = "wall";
const string Particle_Tok = "particle";
const string End_Tok     = "end";

const string Dt_Tok    = "Dt";
const string MinDt_Tok = "MinDt";
const string MaxDt_Tok = "MaxDt";
const string AdjustDt_Tok = "AdjustDt";
const string AdjustDelay_Tok = "AdjustDelay";

#endif // __PARSING_TOKENS_HPP__
