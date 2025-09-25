// SHG.tpp - Template implementations for Spherical Harmonic Gravity
// Part of the Spherical-Harmonic-Gravity (SHG) library

#include <iostream>
#include <cmath>

namespace SHG {

// Template for gravitational acceleration
// T should be a type that can be initialized with {aI, aJ, aK}
template <typename T>
T g(double r, double phi, double lambda, int l_max, int m_max, const std::vector<std::vector<double>> &C, const std::vector<std::vector<double>> &S, double a, double GM)
{
    using namespace SHG;
    std::vector<std::vector<double>> P = Plm_bar(l_max, m_max, phi);
    std::vector<double> mTan = recursive_tangent(m_max, phi);
    auto [sinL, cosL] = recursive_sine_cosine(m_max, lambda);
    double dudr = 1.0;
    double dudphi = 0.0;
    double dudlambda = 0.0;
    double r_ratio_l = a/r; // (a/r)^l will be computed iteratively
    for (int l = 2; l <= l_max; ++l)
    {
        r_ratio_l = r_ratio_l * (a / r); // Update (a/r)^l
        for (int m = 0; m <= l; ++m)
        {
            double SF;
            if (l==m){
                SF = 0.0;
            }
            else if  (m==0){
                SF = sqrt((l+1)*l/2);
            }
            else {
                SF = sqrt((l+m+1)*(l-m));
            }
            double C_lm = C[l][m];
            double S_lm = S[l][m];
            double P_lm = P[l][m];
            double P_lmp1 = P[l][m+1];
            double cos_m_lambda = cosL[m];
            double sin_m_lambda = sinL[m];
            double mtan_phi = mTan[m];
            double C_plus_S = C_lm * cos_m_lambda + S_lm * sin_m_lambda;
            double S_minus_C = S_lm * cos_m_lambda - C_lm * sin_m_lambda;
            dudr += r_ratio_l * (l + 1) * P_lm * C_plus_S;
            dudphi += r_ratio_l * (P_lmp1*SF - mtan_phi * P_lm) * C_plus_S;
            dudlambda += r_ratio_l * m * P_lm * S_minus_C;
            

            if (std::isnan(dudr) || std::isnan(dudphi) || std::isnan(dudlambda))
            {
                std::cerr << "NaN encountered in gravitational acceleration computation at l=" << l << ", m=" << m << std::endl;
            }
        }
    }



    dudr *= -GM / (r * r);
    dudphi *= GM / r;
    dudlambda *= GM / r;

    double ri = r * cos(phi) * cos(lambda);
    double rj = r * cos(phi) * sin(lambda);
    double rk = r * sin(phi);
    double r2 = r * r;
    double ri_rj2 = ri * ri + rj * rj;
    double sqrt_ri_rj2 = sqrt(ri_rj2);
    
    // Handle coordinate singularity at poles (phi = ±90°)
    // When cos(phi) ≈ 0, we need to handle the coordinate transformation carefully
    const double pole_tolerance = 1e-10;
    double cos_phi = cos(phi);
    double aI, aJ, aK;
    
    if (abs(cos_phi) < pole_tolerance) {
        // At the poles: use L'Hôpital's rule and limiting behavior
        // The horizontal accelerations at poles are determined by the gradient in the tangent plane
        // For spherical harmonics, this involves derivatives of Legendre polynomials
        
        // At poles, ri = rj ≈ 0, so horizontal accelerations come from dudlambda/cos(phi)
        // Using L'Hôpital's rule: lim(phi→±π/2) dudlambda/cos(phi) = d(dudlambda)/dphi at pole
        
        // The proper limit gives us the horizontal components based on the azimuthal derivative
        // For numerical stability, we compute this using the derivative relationship
        double sign_phi = (phi > 0) ? 1.0 : -1.0;  // North or South pole
        
        // At poles, the azimuthal acceleration components are:
        aI = -sign_phi * dudlambda / r * sin(lambda);
        aJ = sign_phi * dudlambda / r * cos(lambda);
        aK = dudr / r * rk + sign_phi * dudphi / r;
    } else {
        // Normal case: away from poles
        aI = (1 / r * dudr - rk / (r2 * sqrt_ri_rj2) * dudphi) * ri - (1 / ri_rj2 * dudlambda) * rj;
        aJ = (1 / r * dudr - rk / (r2 * sqrt_ri_rj2) * dudphi) * rj + (1 / ri_rj2 * dudlambda) * ri;
        aK = (1 / r * dudr * rk + sqrt_ri_rj2 / r2 * dudphi);
    }
    T result = {aI, aJ, aK};
    return result;
}

template <typename T>
T gravitational_acceleration_from_cartesian(double x, double y, double z, int l_max, int m_max, const std::vector<std::vector<double>>& C, const std::vector<std::vector<double>>& S, double a, double GM) {
    auto [r, phi, lambda] = cartesian_to_geocentric(x, y, z);
    return SHG::g<T>(r, phi, lambda, l_max, m_max, C, S, a, GM);
}

} // namespace SHG
