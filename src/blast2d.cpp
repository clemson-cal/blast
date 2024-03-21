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
#include "envelope.hpp"




/**
 * 
 */
using namespace vapor;
using cons_t = vec_t<double, 4>;
using prim_t = vec_t<double, 4>;
using cons_array_t = memory_backed_array_t<2, cons_t, ref_counted_ptr_t>;
using prim_array_t = memory_backed_array_t<2, prim_t, ref_counted_ptr_t>;
using Product = memory_backed_array_t<2, double, ref_counted_ptr_t>;

#define gamma_law (4.0 / 3.0)
#define min2(a, b) ((a) < (b) ? (a) : (b))
#define max2(a, b) ((a) > (b) ? (a) : (b))
#define min3(a, b, c) min2(a, min2(b, c))
#define max3(a, b, c) max2(a, max2(b, c))
#define sign(x) copysign(1.0, x)
#define minabs(a, b, c) min3(fabs(a), fabs(b), fabs(c))

// HD static inline double plm_minmod(
//     double yl,
//     double yc,
//     double yr,
//     double plm_theta)
// {
//     double a = (yc - yl) * plm_theta;
//     double b = (yr - yl) * 0.5;
//     double c = (yr - yc) * plm_theta;
//     return 0.25 * fabs(sign(a) + sign(b)) * (sign(a) + sign(c)) * minabs(a, b, c);
// }




/**
 * 
 */
HD static auto gamma_beta_squared(prim_t p) -> double
{
    return p[1] * p[1] + p[2] + p[2];
}

HD static auto momentum_squared(cons_t u) -> double
{
    return u[1] * u[1] + u[2] + u[2];
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
    auto pre = p[3];
    return rho + pre * (1.0 + 1.0 / (gamma_law - 1.0));
}

HD static auto prim_to_cons(prim_t p) -> cons_t
{
    auto rho = p[0];
    auto pre = p[3];
    auto w = lorentz_factor(p);
    auto h = enthalpy_density(p) / rho;
    auto m = rho * w;
    auto u = cons_t{};
    u[0] = m;
    u[1] = m * (h * p[1]);
    u[2] = m * (h * p[2]);
    u[3] = m * (h * w - 1.0) - pre;
    return u;
}

HD static auto cons_to_prim(cons_t cons, double p=0.0) -> optional_t<prim_t>
{
    auto newton_iter_max = 50;
    auto error_tolerance = 1e-12 * (cons[0] + cons[3]);
    auto gm = gamma_law;
    auto m = cons[0];
    auto tau = cons[3];
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
    return some(prim_t{
        m / w0,
        w0 * cons[1] / (tau + m + p),
        w0 * cons[2] / (tau + m + p),
        p
    });
}

HD static auto prim_and_cons_to_flux(prim_t p, cons_t u, int axis) -> cons_t
{
    auto pre = p[3];
    auto vn = beta_component(p, axis);
    auto f = cons_t{};
    f[0] = vn * u[0];
    f[1] = vn * u[1] + pre * (axis == 0);
    f[2] = vn * u[2] + pre * (axis == 1);
    f[3] = vn * u[3] + pre * vn;
    return f;
}

HD static auto sound_speed_squared(prim_t p) -> double
{
    auto pre = p[3];
    auto rho_h = enthalpy_density(p);
    return gamma_law * pre / rho_h;
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

HD static auto riemann_hlle(prim_t pl, prim_t pr, cons_t ul, cons_t ur, int axis, double v_face) -> cons_t
{
    auto fl = prim_and_cons_to_flux(pl, ul, axis);
    auto fr = prim_and_cons_to_flux(pr, ur, axis);
    auto al = outer_wavespeeds(pl, axis);
    auto ar = outer_wavespeeds(pr, axis);
    auto am = min3(0.0, al[0], ar[0]);
    auto ap = max3(0.0, al[1], ar[1]);
    if (v_face < am) {
        return fl - ul * v_face;
    }
    if (v_face > ap) {
        return fr - ur * v_face;
    }
    auto u_hll = (ur * ap - ul * am + (fl - fr)) / (ap - am);
    auto f_hll = (fl * ap - fr * am - (ul - ur) * ap * am) / (ap - am);
    return f_hll - u_hll * v_face;
};

HD static auto spherical_geometry_source_terms(prim_t p, double r0, double r1, double q0, double q1)
{
    // Eqn A8 in Zhang & MacFadyen (2006), integrated over the spherical shell
    // between r0 and r1.
    // 
    auto dcosq = cos(q1) - cos(q0);
    auto dsinq = sin(q1) - sin(q0);
    auto df = 2.0 * M_PI;
    auto ur = p[1];
    auto uq = p[2];
    auto uf = 0.0;
    auto hg = enthalpy_density(p);
    auto pg = p[3];
    auto dr2 = r1 * r1 - r0 * r0;
    auto srdot = -0.5 * df * dr2 * dcosq * (hg * (uq * uq + uf * uf) + 2.0 * pg);
    auto sqdot = +0.5 * df * dr2 * (dcosq * hg * ur * uq + dsinq * (pg + hg * uf * uf));
    return cons_t{
        0.0,
        srdot,
        sqdot,
        0.0,
    };
}




enum class Setup
{
    uniform,
    wind,
    envelope,
};
static Setup setup_from_string(const std::string& name)
{
    if (name == "uniform") return Setup::uniform;
    if (name == "wind") return Setup::wind;
    if (name == "envelope") return Setup::envelope;
    throw std::runtime_error("unknown setup " + name);
}




/**
 * 
 */
struct Config
{
    int fold = 50;
    int rk = 2;
    double theta = 1.5;
    double tstart = 0.0;
    double tfinal = 0.0;
    double cfl = 0.4;
    double cpi = 0.0;
    double spi = 0.0;
    double tsi = 0.0;
    double dx = 1e-2;
    dvec_t<2> domain = {0.1, 0.99};
    bool expand = false;
    std::vector<uint> sp;
    std::vector<uint> ts;
    std::string outdir = ".";
    std::string method = "plm";
    std::string setup = "uniform";
};
VISITABLE_STRUCT(Config,
    fold,
    rk,
    theta,
    tstart,
    tfinal,
    cfl,
    cpi,
    spi,
    tsi,
    dx,
    domain,
    expand,
    sp,
    ts,
    outdir,
    method,
    setup
);




struct initial_model_t
{
    initial_model_t(const Config& config)
    : setup(setup_from_string(config.setup))
    {
    }
    prim_t initial_primitive(double r, double q, double t) const
    {
        switch (setup)
        {
        case Setup::uniform: {
            // Uniform gas (tests spherical geometry source terms)
            //
            return vec(1.0, 0.0, 0.0, 1.0);
        }
        case Setup::wind: {
            // Steady-state cold wind, sub-relativistic velocity
            //
            auto f = 1.0; // mass outflow rate, per steradian, r^2 rho u
            auto u = 0.1; // wind gamma-beta
            auto rho = f / (r * r * u);
            auto pre = 1e-10 * pow(rho, gamma_law); // uniform specfic entropy
            return vec(rho, 0.1, 0.0, pre);
        }
        case Setup::envelope: {
            // A relativistic envelope based on the BNS merger scenario
            // 
            auto envelope = envelope_t();
            auto m = envelope.shell_mass_rt(r, t);
            auto d = envelope.shell_density_mt(m, t);
            auto u = envelope.shell_gamma_beta_m(m);

            // narrow shell (parameters hard-coded for now)
            // 
            if (r < 0.2) {
                d *= 10.0 * exp(-pow(r - 0.2, 2.0) / 0.001) * exp(-q * q / 0.001) + 1.0;
                u += 10.0 * exp(-pow(r - 0.2, 2.0) / 0.001) * exp(-q * q / 0.001);
            }
            auto p = 1e-6 * d;

            return vec(d, u, 0.0, p);
        }
        default: return {};
        }
    }
    Setup setup;
};




struct log_spherical_geometry_t
{
    log_spherical_geometry_t(const Config& config)
    {
        y0 = config.domain[0];
        y1 = config.domain[1];
        q0 = 0.0;
        q1 = M_PI;
        ni = int((log(y1) - log(y0)) / config.dx);
        nj = int((q1 - q0) / config.dx);
        dlogy = (log(y1) - log(y0)) / ni;
        dq = (q1 - q0) / nj;
    }
    HD double face_position_i(int i, double t) const
    {
        return face_velocity_i(i) * t;
    }
    HD double face_position_j(int j) const
    {
        return q0 + dq * j;
    }
    HD double face_velocity_i(int i) const
    {
        return y0 * exp(dlogy * i);
    }
    HD double face_area_i(ivec_t<2> index, double t) const
    {
        auto rp = face_position_i(index[0], t);
        auto rm = face_position_i(index[0], t);
        auto qm = face_position_j(index[1]);
        auto qp = face_position_j(index[1] + 1);
        return area(vec(rm, qm), vec(rp, qp));
    }
    HD double face_area_j(ivec_t<2> index, double t) const
    {
        auto rm = face_position_i(index[0], t);
        auto rp = face_position_i(index[0] + 1, t);
        auto qm = face_position_j(index[1]);
        auto qp = face_position_j(index[1]);
        return area(vec(rm, qm), vec(rp, qp));
    }
    HD dvec_t<2> cell_position(ivec_t<2> index, double t) const
    {
        auto r = 0.5 * (face_position_i(index[0], t) + face_position_i(index[0] + 1, t));
        auto q = 0.5 * (face_position_j(index[1]) + face_position_j(index[1] + 1));
        return {r, q};
    }
    HD double cell_volume(ivec_t<2> index, double t) const
    {
        auto rm = face_position_i(index[0], t);
        auto rp = face_position_i(index[0] + 1, t);
        auto qm = face_position_j(index[1]);
        auto qp = face_position_j(index[1] + 1);
        auto dcost = -(cos(qp) - cos(qm));
        return 2.0 * M_PI * (pow(rp, 3) - pow(rm, 3)) / 3.0 * dcost;
    }
    HD cons_t source_terms(prim_t p, double xm_i, double xp_i, double xm_j, double xp_j) const
    {
        return spherical_geometry_source_terms(p, xm_i, xp_i, xm_j, xp_j);
    }
    index_space_t<2> cells_space() const
    {
        return index_space(ivec(0, 0), uvec(ni, nj));
    }
    static HD double area(dvec_t<2> c0, dvec_t<2> c1)
    {
        auto s0 = c0[0] * sin(c0[1]);
        auto s1 = c1[0] * sin(c1[1]);
        auto z0 = c0[0] * cos(c0[1]);
        auto z1 = c1[0] * cos(c1[1]);
        auto ds = s1 - s0;
        auto dz = z1 - z0;
        return M_PI * (s0 + s1) * sqrt(ds * ds + dz * dz);
    }
    double y0;
    double y1;
    double q0;
    double q1;
    double dlogy;
    double dq;
    int ni;
    int nj;
};




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




template<typename G>
void update_prim(const State& state, G g, prim_array_t& p, double t)
{
    auto u = state.cons;

    if (p.space() != u.space())
    {
        p = zeros<prim_t>(u.space()).cache();
    }
    p = indices(u.space()).map([=] HD (ivec_t<2> i)
    {
        auto dv = g.cell_volume(i, t);
        auto ui = u[i] / dv;
        auto pi = p[i];
        return cons_to_prim(ui, pi[3]).get();
    }).cache();
}




template<class G>
static State next_pcm(const State& state, const Config& config, prim_array_t& p, double dt, int prim_dirty)
{
    auto t = state.time;
    auto u = state.cons;
    auto g = G(config);
    auto cells_space = g.cells_space();
    auto faces_space_i = cells_space.nudge(vec(0, 0), vec(1, 0));
    auto faces_space_j = cells_space.nudge(vec(0, 0), vec(0, 1));
    auto vf_i = indices(faces_space_i).map([=] HD (ivec_t<2> i) { return g.face_velocity_i(i[0]); });
    auto xf_i = indices(faces_space_i).map([=] HD (ivec_t<2> i) { return g.face_position_i(i[0], t); });
    auto xf_j = indices(faces_space_j).map([=] HD (ivec_t<2> i) { return g.face_position_j(i[1]); });
    auto da_i = indices(faces_space_i).map([=] HD (ivec_t<2> i) { return g.face_area_i(i, t); });
    auto da_j = indices(faces_space_j).map([=] HD (ivec_t<2> i) { return g.face_area_j(i, t); });
    auto dv = indices(cells_space).map([=] HD (ivec_t<2> i) { return g.cell_volume(i, t); });

    if (prim_dirty) {
        update_prim(state, g, p, t);
    }

    auto model = initial_model_t(config);
    auto model_radial_fluxes = indices(faces_space_i).map([=] HD (ivec_t<2> i)
    {
        auto rf = xf_i[i];
        auto p = model.initial_primitive(rf, 0.0, t); // assumes no polar dependance of the initial state
        auto u = prim_to_cons(p);
        return prim_and_cons_to_flux(p, u, 0) - u * vf_i[i];
    });

    auto fhat_i = model_radial_fluxes.insert(indices(faces_space_i.contract(uvec(1, 0))).map([=] HD (ivec_t<2> i)
    {
        auto ul = u[i - vec(1, 0)] / dv[i - vec(1, 0)];
        auto ur = u[i - vec(0, 0)] / dv[i - vec(0, 0)];
        auto pl = p[i - vec(1, 0)];
        auto pr = p[i - vec(0, 0)];
        return riemann_hlle(pl, pr, ul, ur, 0, vf_i[i]);
    })).cache();

    auto fhat_j = zeros<cons_t>(faces_space_j).insert(indices(faces_space_j.contract(uvec(0, 1))).map([=] HD (ivec_t<2> i)
    {
        auto ul = u[i - vec(0, 1)] / dv[i - vec(0, 1)];
        auto ur = u[i - vec(0, 0)] / dv[i - vec(0, 0)];
        auto pl = p[i - vec(0, 1)];
        auto pr = p[i - vec(0, 0)];
        return riemann_hlle(pl, pr, ul, ur, 1, 0.0);
    })).cache();

    auto du = indices(cells_space).map([=] HD (ivec_t<2> i)
    {
        auto rm_i = xf_i[i + vec(0, 0)];
        auto rp_i = xf_i[i + vec(1, 0)];
        auto rm_j = xf_j[i + vec(0, 0)];
        auto rp_j = xf_j[i + vec(0, 1)];
        auto am_i = da_i[i + vec(0, 0)];
        auto ap_i = da_i[i + vec(1, 0)];
        auto am_j = da_j[i + vec(0, 0)];
        auto ap_j = da_j[i + vec(0, 1)];
        auto fm_i = fhat_i[i + vec(0, 0)];
        auto fp_i = fhat_i[i + vec(1, 0)];
        auto fm_j = fhat_j[i + vec(0, 0)];
        auto fp_j = fhat_j[i + vec(0, 1)];
        auto udot = g.source_terms(p[i], rm_i, rp_i, rm_j, rp_j);
        return fm_i * am_i - fp_i * ap_i + fm_j * am_j - fp_j * ap_j + udot;
    }) * dt;

    return State{
        state.time + dt,
        state.iter + 1.0,
        (u + du).cache(),
    };
}




template<class G>
static State next_plm(const State& state, const Config& config, prim_array_t& p, double dt, int prim_dirty)
{
    return state;
    // auto u = state.cons;
    // auto g = G(config, state.time);
    // auto ic = range(u.space());
    // auto iv = range(u.space().nudge(vec(0), vec(1)));
    // auto rf = iv.map([g] HD (int i) { return g.face_position(i); });
    // auto vf = iv.map([g] HD (int i) { return g.face_velocity(i); });
    // auto da = iv.map([g] HD (int i) { return g.face_area(i); });
    // auto gradient_cells = ic.space().contract(1);
    // auto interior_cells = ic.space().contract(2);
    // auto interior_faces = iv.space().contract(2);
    // auto plm_theta = config.theta;

    // if (prim_dirty) {
    //     update_prim(state, g, p);
    // }

    // auto grad = ic[gradient_cells].map([p, plm_theta] HD (int i)
    // {
    //     auto pl = p[i - 1];
    //     auto pc = p[i + 0];
    //     auto pr = p[i + 1];
    //     auto gc = prim_t{};

    //     for (int n = 0; n < 3; ++n)
    //     {
    //         gc[n] = plm_minmod(pl[n], pc[n], pr[n], plm_theta);
    //     }
    //     return gc;
    // }).cache();

    // auto fhat = iv[interior_faces].map([p, vf, grad] HD (int i)
    // {
    //     auto pl = p[i - 1] + grad[i - 1] * 0.5;
    //     auto pr = p[i + 0] - grad[i + 0] * 0.5;
    //     auto ul = prim_to_cons(pl);
    //     auto ur = prim_to_cons(pr);
    //     return riemann_hlle(pl, pr, ul, ur, vf[i]);
    // }).cache();

    // auto du = ic[interior_cells].map([rf, da, g, fhat, p] HD (int i)
    // {
    //     auto rm = rf[i + 0];
    //     auto rp = rf[i + 1];
    //     auto am = da[i + 0];
    //     auto ap = da[i + 1];
    //     auto fm = fhat[i + 0];
    //     auto fp = fhat[i + 1];
    //     auto udot = g.source_terms(p[i], rm, rp);
    //     return fm * am - fp * ap + udot;
    // }) * dt;

    // return State{
    //     state.time + dt,
    //     state.iter + 1.0,
    //     set_bc(u.add(du), config, 2),
    // };
}




template<class Geometry>
static void update_state(State& state, const Config& config)
{
    static prim_array_t p;
    auto next = std::function<State(State&, const Config&, prim_array_t&, double, int)>();

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

    auto t = state.time;
    auto g = Geometry(config);
    auto dx = g.face_position_i(1, t) - g.face_position_i(0, t);
    auto dt = dx * config.cfl;
    auto s0 = state;

    switch (config.rk)
    {
        case 1: {
            state = next(s0, config, p, dt, 1);
            break;
        }
        case 2: {
            auto s1 = average(s0, next(s0, config, p, dt, 1), 1./1);
            auto s2 = average(s0, next(s1, config, p, dt, 1), 1./2);
            state = s2;
            break;
        }
        case 3: {
            auto s1 = average(s0, next(s0, config, p, dt, 1), 1./1);
            auto s2 = average(s0, next(s1, config, p, dt, 1), 1./4);
            auto s3 = average(s0, next(s2, config, p, dt, 1), 2./3);
            state = s3;
            break;
        }
    }
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
        auto model = initial_model_t(config);
        auto t = config.tstart;
        auto g = log_spherical_geometry_t(config);
        auto u = indices(g.cells_space()).map([=] HD (ivec_t<2> i) {
            auto dv = g.cell_volume(i, t);
            auto xc = g.cell_position(i, t);
            auto p = model.initial_primitive(xc[0], xc[1], t);
            auto u = prim_to_cons(p);
            return u * dv;
        });
        state.time = t;
        state.iter = 0.0;
        state.cons = u.cache();
    }
    void update(State& state) const override
    {
        return update_state<log_spherical_geometry_t>(state, config);
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
        case 1: return "radial_gamma_beta";
        case 2: return "gas_pressure";
        case 3: return "face_positions_i";
        case 4: return "face_positions_j";
        }
        return nullptr;
    }
    Product compute_product(const State& state, uint column) const override
    {
        auto t = state.time;
        auto g = log_spherical_geometry_t(config);
        auto ic = indices(g.cells_space());
        auto xf_i = indices(index_space(vec(0, 0), uvec(g.ni + 1, 1))).map([=] HD (ivec_t<2> i) { return g.face_position_i(i[0], t); });
        auto xf_j = indices(index_space(vec(0, 0), uvec(1, g.nj + 1))).map([=] HD (ivec_t<2> i) { return g.face_position_j(i[1]); });
        auto dv = ic.map([g, t] HD (ivec_t<2> i) { return g.cell_volume(i, t); });
        auto cons_field = [] (uint n) {
            return [n] HD (cons_t u) {
                return cons_to_prim(u).get()[n];
            };
        };
        switch (column) {
        case 0: return (state.cons / dv).map(cons_field(0)).cache();
        case 1: return (state.cons / dv).map(cons_field(1)).cache();
        case 2: return (state.cons / dv).map(cons_field(3)).cache();
        case 3: return xf_i.cache();
        case 4: return xf_j.cache();
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
