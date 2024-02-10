/**
================================================================================
Copyright 2024, Jonathan Zrake

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
================================================================================
*/

#include <cmath>
#include <functional>
#include "vapor/vapor.hpp"




/**
 * 
 */
using namespace vapor;
using cons_t = vec_t<double, 3>;
using prim_t = vec_t<double, 3>;
using cons_array_t = memory_backed_array_t<1, cons_t, ref_counted_ptr_t>;
using prim_array_t = memory_backed_array_t<1, prim_t, ref_counted_ptr_t>;
using Product = memory_backed_array_t<1, double, ref_counted_ptr_t>;

#define gamma (4.0 / 3.0)
#define min2(a, b) ((a) < (b) ? (a) : (b))
#define max2(a, b) ((a) > (b) ? (a) : (b))
#define min3(a, b, c) min2(a, min2(b, c))
#define max3(a, b, c) max2(a, max2(b, c))
#define sign(x) copysign(1.0, x)
#define minabs(a, b, c) min3(fabs(a), fabs(b), fabs(c))

HD static inline double plm_minmod(
    double yl,
    double yc,
    double yr,
    double plm_theta)
{
    double a = (yc - yl) * plm_theta;
    double b = (yr - yl) * 0.5;
    double c = (yr - yc) * plm_theta;
    return 0.25 * fabs(sign(a) + sign(b)) * (sign(a) + sign(c)) * minabs(a, b, c);
}




/**
 * 
 */
HD static auto gamma_beta_squared(prim_t p) -> double
{
    return p[1] * p[1];
}

HD static auto momentum_squared(cons_t u) -> double
{
    return u[1] * u[1];
}

HD static auto lorentz_factor(prim_t p) -> double
{
    return sqrt(1.0 + gamma_beta_squared(p));
}

HD static auto beta_component(prim_t p, int axis) -> double
{
    return p[axis + 1] / lorentz_factor(p);
}

HD static auto enthalpy_density(prim_t p) -> double
{
    auto rho = p[0];
    auto pre = p[2];
    return rho + pre * (1.0 + 1.0 / (gamma - 1.0));
}

HD static auto prim_to_cons(prim_t p) -> cons_t
{
    auto rho = p[0];
    auto pre = p[2];
    auto w = lorentz_factor(p);
    auto h = enthalpy_density(p) / rho;
    auto m = rho * w;
    auto u = cons_t{};
    u[0] = m;
    u[1] = m * (h * p[1]);
    u[2] = m * (h * w - 1.0) - pre;
    return u;
}

HD static auto cons_to_prim(cons_t cons, double p=0.0) -> optional_t<prim_t>
{
    auto newton_iter_max = 50;
    auto error_tolerance = 1e-12 * (cons[0] + cons[2]);
    auto gm = gamma;
    auto m = cons[0];
    auto tau = cons[2];
    auto ss = momentum_squared(cons);
    auto n = 0;
    auto w0 = 0.0;

    while (true)
    {
        auto et = tau + p + m;
        auto b2 = min2(ss / et / et, 1.0 - 1e-10);
        auto w2 = 1.0 / (1.0 - b2);
        auto w = sqrt(w2);
        auto d = m / w;
        auto de = (tau + m * (1.0 - w) + p * (1.0 - w2)) / w2;
        auto dh = d + de + p;
        auto a2 = gm * p / dh;
        auto g = b2 * a2 - 1.0;
        auto f = de * (gm - 1.0) - p;

        if (n == newton_iter_max)
        {
            // WARNING: none values are silently unwrapped in the code below,
            // because the cache_unwrap facility in vapor is currently broken.
            return none<prim_t>();
        }
        if (fabs(f) < error_tolerance)
        {
            w0 = w;
            break;
        }        
        p -= f / g;
        n += 1;
    }
    return some(prim_t{m / w0, w0 * cons[1] / (tau + m + p), p});
}

HD static auto prim_and_cons_to_flux(prim_t p, cons_t u, int axis) -> cons_t
{
    auto pre = p[2];
    auto vn = beta_component(p, axis);
    auto f = cons_t{};
    f[0] = vn * u[0];
    f[1] = vn * u[1] + pre * (axis == 0);
    f[2] = vn * u[2] + pre * vn;
    return f;
}

HD static auto sound_speed_squared(prim_t p) -> double
{
    auto pre = p[2];
    auto rho_h = enthalpy_density(p);
    return gamma * pre / rho_h;
}

HD static auto outer_wavespeeds(prim_t p, int axis) -> dvec_t<2>
{
    auto a2 = sound_speed_squared(p);
    auto uu = gamma_beta_squared(p);
    auto vn = beta_component(p, axis);
    auto v2 = vn * vn;
    // These are the Marti & Muller versions of the eigenvalue formula:
    // 
    // auto vv = uu / (1.0 + uu);
    // auto k0 = sqrt(a2 * (1.0 - vv) * (1.0 - vv * a2 - v2 * (1.0 - a2)));
    // auto xx = vec(
    //     (vn * (1.0 - a2) - k0) / (1.0 - vv * a2),
    //     (vn * (1.0 - a2) + k0) / (1.0 - vv * a2)
    // );
    // These are the Mignone & Bodo (2005) versions of the eigenvalue formula:
    auto g2 = 1.0 + uu;
    auto s2 = a2 / g2 / (1.0 - a2);
    auto k0 = sqrt(s2 * (1.0 - v2 + s2));
    return vec(vn - k0, vn + k0) / (1.0 + s2);
}

HD static auto riemann_hlle(prim_t pl, prim_t pr, cons_t ul, cons_t ur) -> cons_t
{
    auto fl = prim_and_cons_to_flux(pl, ul, 0);
    auto fr = prim_and_cons_to_flux(pr, ur, 0);
    // These wave speed estimates are technically valid in 3d:
    // 
    // auto al = outer_wavespeeds(pl, 0);
    // auto ar = outer_wavespeeds(pr, 0);
    // auto am = min3(0.0, al[0], ar[0]);
    // auto ap = max3(0.0, al[1], ar[1]);
    // 
    // These wave speed estimates should be equivalent to those above, for 1d:
    auto cl = sqrt(sound_speed_squared(pl));
    auto cr = sqrt(sound_speed_squared(pr));
    auto vl = beta_component(pl, 0);
    auto vr = beta_component(pr, 0);
    auto alm = (vl - cl) / (1.0 - vl * cl);
    auto alp = (vl + cl) / (1.0 + vl * cl);
    auto arm = (vr - cr) / (1.0 - vr * cr);
    auto arp = (vr + cr) / (1.0 + vr * cr);
    auto am = min3(0.0, alm, arm);
    auto ap = max3(0.0, alp, arp);
    (void)outer_wavespeeds; // silence unused function warning
    // These wave speed estimates are from Schneider et al. (1993):
    //
    // auto vb = 0.5 * (vl + vr);
    // auto cs = 0.5 * (cl + cr);
    // auto am = min2(0.0, (vb - cs) / (1.0 - vb * cs));
    // auto ap = max2(0.0, (vb + cs) / (1.0 + vb * cs));
    return (fl * ap - fr * am - (ul - ur) * ap * am) / (ap - am);
};

HD static auto spherical_geometry_source_terms(prim_t p, double r0, double r1)
{
    // Eqn A8 in Zhang & MacFadyen (2006), integrated over the spherical shell
    // between r0 and r1, and specializing to radial velocity only.
    auto pg = p[2];
    auto dr2 = pow(r1, 2) - pow(r0, 2);
    auto srdot = pg * dr2;
    return cons_t{0.0, srdot, 0.0};
}




template<class X, class A, class V, class S>
struct grid_geometry_t
{
    array_t<1, X> face_position;
    array_t<1, A> face_area;
    array_t<1, V> cell_volume;
    S geometric_source_terms;
};

template<class X, class A, class V, class S>
auto grid_geometry(
    array_t<1, X> face_position,
    array_t<1, A> face_area,
    array_t<1, V> cell_volume,
    S geometric_source_terms)
{
    return grid_geometry_t<X, A, V, S>{
        face_position,
        face_area,
        cell_volume,
        geometric_source_terms,
    };
}




enum class Setup
{
    uniform,
    sod,
    mm96p1,
    wind,
    bmk,
    thermal_bomb,
    fast_shell,
};
static Setup setup_from_string(const std::string& name)
{
    if (name == "uniform") return Setup::uniform;
    if (name == "sod") return Setup::sod;
    if (name == "mm96p1") return Setup::mm96p1;
    if (name == "wind") return Setup::wind;
    if (name == "bmk") return Setup::bmk;
    if (name == "thermal_bomb") return Setup::thermal_bomb;
    if (name == "fast_shell") return Setup::fast_shell;
    throw std::runtime_error("unknown setup " + name);
}

enum class CoordinateSystem
{
    planar,
    spherical,
};
static CoordinateSystem coords_from_string(const std::string& name)
{
    if (name == "planar") return CoordinateSystem::planar;
    if (name == "spherical") return CoordinateSystem::spherical;
    throw std::runtime_error("unknown coords " + name);
}




/**
 * 
 */
struct Config
{
    int fold = 50;
    int rk = 2;
    double theta = 1.5;
    double tfinal = 0.0;
    double cfl = 0.4;
    double cpi = 0.0;
    double spi = 0.0;
    double tsi = 0.0;
    double dx = 1e-2;
    double bmk_gamma_shock = 5.0;
    vec_t<double, 2> domain = {{1.0, 10.0}}; // x0, x1
    vec_t<char, 2> bc = {{'f', 'f'}};
    std::vector<uint> sp = {0, 1, 2, 3};
    std::vector<uint> ts;
    std::string outdir = ".";
    std::string method = "plm";
    std::string setup = "uniform";
    std::string coords = "spherical"; // or planar
};
VISITABLE_STRUCT(Config,
    fold,
    rk,
    theta,
    tfinal,
    cfl,
    cpi,
    spi,
    tsi,
    dx,
    bmk_gamma_shock,
    domain,
    bc,
    sp,
    ts,
    outdir,
    method,
    setup,
    coords
);




/**
 * 
 */
struct State
{
    double time;
    double iter;
    cons_array_t cons;
};
VISITABLE_STRUCT(State, time, iter, cons);




static State average(const State& a, const State& b, double x)
{
    return x == 1.0 ? b : State{
        (a.time * (1.0 - x) + b.time * x),
        (a.iter * (1.0 - x) + b.iter * x),
        (a.cons * (1.0 - x) + b.cons * x).cache()
    };
}




static void update_prim(const State& state, prim_array_t& p)
{
    auto u = state.cons;

    if (p.space() != u.space())
    {
        p = zeros<prim_t>(u.space()).cache();
    }
    p = range(u.shape()[0]).map([p, u] HD (int i)
    {
        auto ui = u[i];
        auto pi = p[i];
        return cons_to_prim(ui, pi[2]).get();
    }).cache();
}




template<class F>
auto set_bc(const array_t<1, F> &u, const Config& config, int ng)
{
    auto bcl = config.bc[0];
    auto bcr = config.bc[1];
    auto il = index_space(ivec(0), uvec(ng));
    auto ir = index_space(ivec(u.size() - ng), uvec(ng));
    auto ul = cons_t{};
    auto ur = cons_t{};

    if (bcl == 'f') {
        ul = u[ng];
    } else {
        ul = u[ng - 1];
    }
    if (bcr == 'f') {
        ur = u[u.size() - ng];
    } else {
        ur = u[u.size() - ng - 1];
    }
    if (bcl == 'r') {
        ul[1] *= -1.0;
    }
    if (bcr == 'r') {
        ur[1] *= -1.0;
    }
    return u.at(il).set(ul).at(ir).set(ur);
}




template<class G>
static State next_pcm(const State& state, const G& geom, const Config& config, prim_array_t& p, double dt, int prim_dirty)
{
    auto u = state.cons;
    auto rf = geom.face_position;
    auto da = geom.face_area;
    auto dv = geom.cell_volume;
    auto st = geom.geometric_source_terms;
    auto ic = range(dv.space());
    auto iv = range(rf.space());
    auto interior_cells = ic.space().contract(2);
    auto interior_faces = iv.space().contract(2);

    if (prim_dirty) {
        update_prim(state, p);
    }

    auto fhat = iv[interior_faces].map([p, u] HD (int i)
    {
        auto ul = u[i - 1];
        auto ur = u[i];
        auto pl = p[i - 1];
        auto pr = p[i];
        return riemann_hlle(pl, pr, ul, ur);
    }).cache();

    auto du = ic[interior_cells].map([rf, da, dv, st, fhat, p] HD (int i)
    {
        auto rm = rf[i + 0];
        auto rp = rf[i + 1];
        auto am = da[i + 0];
        auto ap = da[i + 1];
        auto fm = fhat[i + 0];
        auto fp = fhat[i + 1];
        auto udot = st(p[i], rm, rp);
        return (fm * am - fp * ap + udot) / dv[i];
    }) * dt;

    return State{
        state.time + dt,
        state.iter + 1.0,
        set_bc(u.at(interior_cells) + du, config, 1).cache(),
    };
    return state;
}




template<class G>
static State next_plm(const State& state, const G& geom, const Config& config, prim_array_t& p, double dt, int prim_dirty)
{
    auto u = state.cons;
    auto rf = geom.face_position;
    auto da = geom.face_area;
    auto dv = geom.cell_volume;
    auto st = geom.geometric_source_terms;
    auto ic = range(dv.space());
    auto iv = range(rf.space());
    auto gradient_cells = ic.space().contract(1);
    auto interior_cells = ic.space().contract(2);
    auto interior_faces = iv.space().contract(2);
    auto plm_theta = config.theta;

    if (prim_dirty) {
        update_prim(state, p);
    }

    auto grad = ic[gradient_cells].map([p, plm_theta] HD (int i)
    {
        auto pl = p[i - 1];
        auto pc = p[i + 0];
        auto pr = p[i + 1];
        auto gc = prim_t{};

        for (int n = 0; n < 3; ++n)
        {
            gc[n] = plm_minmod(pl[n], pc[n], pr[n], plm_theta);
        }
        return gc;
    }).cache();

    auto fhat = iv[interior_faces].map([p, grad] HD (int i)
    {
        auto pl = p[i - 1] + grad[i - 1] * 0.5;
        auto pr = p[i + 0] - grad[i + 0] * 0.5;
        auto ul = prim_to_cons(pl);
        auto ur = prim_to_cons(pr);
        return riemann_hlle(pl, pr, ul, ur);
    }).cache();

    auto du = ic[interior_cells].map([rf, da, dv, st, fhat, p] HD (int i)
    {
        auto rm = rf[i + 0];
        auto rp = rf[i + 1];
        auto am = da[i + 0];
        auto ap = da[i + 1];
        auto fm = fhat[i + 0];
        auto fp = fhat[i + 1];
        auto udot = st(p[i], rm, rp);
        return (fm * am - fp * ap + udot) / dv[i];
    }) * dt;

    return State{
        state.time + dt,
        state.iter + 1.0,
        set_bc(u.at(interior_cells) + du, config, 2).cache(),
    };
}




template<class Geometry>
static void update_state(State& state, const Geometry& geom, const Config& config)
{
    static prim_array_t p;
    auto next = std::function<State(State&, const Geometry&, const Config&, prim_array_t&, double, int)>();

    if (config.method == "pcm") {
        next = next_pcm<Geometry>;
    }
    if (config.method == "plm") {
        next = next_plm<Geometry>;
    }
    if (! next) {
        throw std::runtime_error(format("unrecognized method '%s'", config.method.data()));
    }

    // Wavespeed calculation is presently disabled, max wave speed of c is
    // assumed.
    //
    // auto wavespeed = [] HD (prim_t p) -> double
    // {
    //     auto a = outer_wavespeeds(p, 0);
    //     return max2(fabs(a[0]), fabs(a[1]));
    // };
    // update_prim(state, p);

    auto dx = config.dx;
    auto dt = dx * config.cfl;
    auto s0 = state;

    switch (config.rk)
    {
        case 1: {
            state = next(s0, geom, config, p, dt, 1);
            break;
        }
        case 2: {
            auto s1 = average(s0, next(s0, geom, config, p, dt, 1), 1./1);
            auto s2 = average(s0, next(s1, geom, config, p, dt, 1), 1./2);
            state = s2;
            break;
        }
        case 3: {
            auto s1 = average(s0, next(s0, geom, config, p, dt, 1), 1./1);
            auto s2 = average(s0, next(s1, geom, config, p, dt, 1), 1./4);
            auto s3 = average(s0, next(s2, geom, config, p, dt, 1), 2./3);
            state = s3;
            break;
        }
    }
}




static auto cell_coordinates(const Config& config)
{
    auto x0 = config.domain[0];
    auto x1 = config.domain[1];
    auto dx = config.dx;
    auto ni = int((x1 - x0) / dx);
    auto ic = range(ni);
    auto xc = (ic + 0.5) * dx + x0;
    return xc;
}

static auto face_coordinates(const Config& config)
{
    auto x0 = config.domain[0];
    auto x1 = config.domain[1];
    auto dx = config.dx;
    auto ni = int((x1 - x0) / dx);
    auto ic = range(ni + 1);
    auto xc = ic * dx + x0;
    return xc;
}




/**
 * 
 */
class Blast : public Simulation<Config, State, Product>
{
public:
    const char* name() const override
    {
        return "blast";
    }
    const char* author() const override
    {
        return "Jonathan Zrake (Clemson)";
    }
    const char* description() const override
    {
        return "GPU-accelerated hydrodynamics code for relativistic explosions";
    }
    const char* output_directory() const override
    {
        return config.outdir.data();
    }
    bool use_persistent_session() const override
    {
        return false;
    }
    double get_time(const State& state) const override
    {
        return state.time;
    }
    uint get_iteration(const State& state) const override
    {
        return round(state.iter);
    }
    void initial_state(State& state) const override
    {
        auto setup = setup_from_string(config.setup);
        auto x0 = config.domain[0];
        auto x1 = config.domain[1];
        auto bmk_gamma_shock = config.bmk_gamma_shock;
        auto initial_primitive = [=] HD (double x) -> prim_t
        {
            switch (setup)
            {
            case Setup::uniform: {
                    // Uniform gas (tests spherical geometry source terms)
                    //
                    return vec(1.0, 0.0, 1.0);
                }
            case Setup::sod: {
                    // Standard Sod shocktube (easy problem)
                    //
                    if (x < x0 + 0.5 * (x1 - x0)) {
                        return vec(1.0, 0.0, 1.0);
                    } else {
                        return vec(0.1, 0.0, 0.125);
                    }
                }
            case Setup::mm96p1: {
                    // Problem 1 from Marti & Muller 1996
                    //
                    if (x < x0 + 0.5 * (x1 - x0)) {
                        return vec(10.0, 0.0, 13.33);
                    } else {
                        return vec(1.0, 0.0, 1e-8);
                    }
                }
            case Setup::wind: {
                    // Steady-state cold wind, sub-relativistic velocity
                    //
                    auto f = 1.0; // mass outflow rate, per steradian, r^2 rho u
                    auto u = 0.1; // wind gamma-beta
                    auto rho = f / (x * x * u);
                    auto pre = 1e-10 * pow(rho, gamma); // uniform specfic entropy
                    return vec(rho, 0.1, pre);
                }
            case Setup::bmk: {
                    // Blandford-Mckee ultra-relativistic blast wave
                    //
                    auto shock_radius = 1.0;
                    auto eta = x / shock_radius;
                    if (eta < 1.0) {
                        auto Gamma = bmk_gamma_shock;
                        auto xi = Gamma * Gamma * (eta - 1.0);
                        auto f = 0.5 - 8.0 * xi;
                        auto g = 2.0 * sqrt(2.0) * pow((1.0 - 8.0 * xi), (-5.0 / 4.0));
                        auto h = (2.0 / 3.0) * pow((1.0 - 8.0 * xi), (-17.0 / 12.0));
                        auto rho = Gamma * g;
                        auto vsh = sqrt(1.0 - 1.0 / Gamma / Gamma);
                        auto vel = vsh * (1.0 - (1.0 / Gamma / Gamma) * f);
                        auto gb = vel / sqrt(1.0 - vel * vel);
                        auto pre = Gamma * Gamma * h;
                        return vec(rho, gb, pre);
                    } else {
                        return vec(1.0, 0.0, 1e-8);
                    }
                }
            case Setup::thermal_bomb: {
                auto energy = 1.0;
                auto r_in = 0.1 * (x1 - x0);
                auto rho_in = 1.0;
                auto rho_out = 1.0;
                auto p_in = energy / (4.0 / 3.0 * 3.14159 * pow(r_in, 3.0)) * (gamma - 1.0);
                auto p_out = rho_out * 1.0e-6;
                if (x < r_in) {
                    return vec(rho_in, 0.0, p_in);
                } else {
                    return vec(rho_out, 0.0, p_out);
                }
            }
            case Setup::fast_shell: {
                auto r_in = 0.1 * (x1 - x0);
                auto rho_in = x;
                auto rho_out = 1.0;
                if (x < r_in) {
                    auto bg = x * 30.0;
                    auto p = rho_in * 1.0E-1;
                    return vec(rho_in, bg, p);
                } else {
                    return vec(rho_out, 0.0, rho_out * 1.0e-6);
                }
            }
            default: return {};
            }
        };
        auto initial_conserved = [=] HD (double x) { return prim_to_cons(initial_primitive(x)); };
        state.time = 0.0;
        state.iter = 0.0;
        state.cons = cell_coordinates(config).map(initial_conserved).cache();
    }
    void update(State& state) const override
    {
        switch (coords_from_string(config.coords))
        {
        case CoordinateSystem::planar: return update_planar(state);
        case CoordinateSystem::spherical: return update_spherical(state);
        }
    }
    void update_planar(State& state) const
    {
        auto xc = cell_coordinates(config);
        auto xf = face_coordinates(config);
        auto ic = range(xc.space());
        auto da = xf.map([] HD (double) -> double { return 1.0; });
        auto dv = ic.map([xf] HD (int i) -> double { return xf[i + 1] - xf[i]; });
        auto st = [] HD (prim_t p, double rm, double rp) -> cons_t { return cons_t{}; };
        auto geom = grid_geometry(xf, da, dv, st);
        update_state(state, geom, config);
    }
    void update_spherical(State& state) const
    {
        auto rc = cell_coordinates(config);
        auto rf = face_coordinates(config);
        auto ic = range(rc.space());
        auto da = rf.map([] HD (double r) -> double { return r * r; });
        auto dv = ic.map([rf] HD (int i) -> double
        {
            auto x0 = rf[i + 0];
            auto r1 = rf[i + 1];
            return (r1 * r1 * r1 - x0 * x0 * x0) / 3.0;
        });
        auto st = [] HD (prim_t p, double rm, double rp) -> cons_t
        {
            return spherical_geometry_source_terms(p, rm, rp);
        };
        auto geom = grid_geometry(rf, da, dv, st);
        update_state(state, geom, config);
    }
    bool should_continue(const State& state) const override
    {
        return state.time < config.tfinal;
    }
    uint updates_per_batch() const override
    {
        return config.fold;
    }
    double checkpoint_interval() const override
    { 
        return config.cpi;
    }
    double product_interval() const override
    {
        return config.spi;
    }
    double timeseries_interval() const override
    {
        return config.tsi;
    }
    std::set<uint> get_product_cols() const override
    {
        return std::set(config.sp.begin(), config.sp.end());
    }
    const char* get_product_name(uint column) const override
    {
        switch (column) {
        case 0: return "comoving_mass_density";
        case 1: return "gamma_beta";
        case 2: return "gas_pressure";
        case 3: return "cell_coordinate";
        }
        return nullptr;
    }
    Product compute_product(const State& state, uint column) const override
    {
        auto cons_field = [] (uint n) {
            return [n] HD (cons_t u) {
                return cons_to_prim(u).get()[n];
            };
        };
        switch (column) {
        case 0: return state.cons.map(cons_field(0)).cache();
        case 1: return state.cons.map(cons_field(1)).cache();
        case 2: return state.cons.map(cons_field(2)).cache();
        case 3: return cell_coordinates(config).cache();
        }
        return {};
    }
    std::set<uint> get_timeseries_cols() const override
    {
        return std::set(config.ts.begin(), config.ts.end());
    }
    const char* get_timeseries_name(uint column) const override
    {
        switch (column) {
        case 0: return "time";
        }
        return nullptr;
    }
    double compute_timeseries_sample(const State& state, uint column) const override
    {
        return state.time;
    }
    vec_t<char, 256> status_message(const State& state, double secs_per_update) const override
    {
        return format("[%04d] t=%lf Mzps=%.2lf",
            get_iteration(state),
            get_time(state),
            1e-6 * state.cons.size() / secs_per_update);
    }
};




int main(int argc, const char **argv)
{
    try {
        return Blast().run(argc, argv);
    }
    catch (const std::exception& e) {
        vapor::print(e.what(), "\n");
    }
    return 0;
}
