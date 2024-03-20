#include <cmath>
#include "vapor/vapor.hpp"




/**
 * A relativistic envelope based on the BNS merger scenario of Andrei Beloborodov
 */
struct envelope_t
{
    HD double shell_time_m(double m) const
    {
        return m / mdot_wind();
    }
    HD double shell_time_mprime(double m) const
    {
        return 1.0 / mdot_wind();
    }
    HD double shell_gamma_beta_m(double m) const
    {
        return u_wind + pow(m / m1, -psi);
    }
    HD double shell_gamma_beta_mprime(double m) const
    {
        return -pow(m / m1, -psi) * psi / m;
    }
    HD double shell_speed_m(double m) const
    {
        auto u = shell_gamma_beta_m(m);
        return u / sqrt(1 + u * u);
    }
    HD double shell_speed_mprime(double m) const
    {
        auto u = shell_gamma_beta_m(m);
        auto du_dm = shell_gamma_beta_mprime(m);
        auto dv_du = pow(1 + u * u, -1.5);
        return dv_du * du_dm;
    }
    HD double shell_radius_mt(double m, double t) const
    {
        auto v = shell_speed_m(m);
        auto t0 = shell_time_m(m);
        return v * (t - t0);
    }
    HD double shell_density_mt(double m, double t) const
    {
        auto t0 = shell_time_m(m);
        auto t0_prime = shell_time_mprime(m);
        auto u = shell_gamma_beta_m(m);
        auto u_prime = shell_gamma_beta_mprime(m);
        auto mdot_inverse = t0_prime - u_prime / u * (t - t0) / (1 + u * u);
        auto r = shell_radius_mt(m, t);
        return 1.0 / (4 * M_PI * r * r * u * mdot_inverse);
    }
    HD double shell_mass_rt(double r, double t) const
    {
        auto f = [this, r, t] (double m)
        {
            return r - shell_radius_mt(m, t);
        };
        auto g = [this, t] (double m)
        {
            auto v = shell_speed_m(m);
            auto t0 = shell_time_m(m);
            auto dv = shell_speed_mprime(m);
            auto dt0 = shell_time_mprime(m);
            return -(dv * (t - t0) - v * dt0);
        };
        auto m = 1e-12;
        auto n = 0;
        while (true)
        {
            auto fm = f(m);
            auto gm = g(m);

            m -= fm / gm;
            if (fabs(fm) < 1e-10) {
                return m;
            }
            if (n > 100) {
                printf("[warning envelope_t::shell_mass_rt] could not find mass coordinate\n");
                return m;
            }
            n += 1;
        }
    }
    HD double mdot_wind() const
    {
        return m_cloud / t_delay;
    }
    double t_delay = 1.0;
    double m1 = 1.0;
    double psi = 0.25;
    double m_cloud = 1e5;
    double u_wind = 0.1;
};
