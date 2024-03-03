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

#define gamma_law (4.0 / 3.0)
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
    return rho + pre * (1.0 + 1.0 / (gamma_law - 1.0));
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
    auto gm = gamma_law;
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

HD static auto riemann_hlle(prim_t pl, prim_t pr, cons_t ul, cons_t ur, double v_face) -> cons_t
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
    auto am = min2(alm, arm);
    auto ap = max2(alp, arp);
    (void)outer_wavespeeds; // silence unused function warning
    // These wave speed estimates are from Schneider et al. (1993):
    //
    // auto vb = 0.5 * (vl + vr);
    // auto cs = 0.5 * (cl + cr);
    // auto am = min2(0.0, (vb - cs) / (1.0 - vb * cs));
    // auto ap = max2(0.0, (vb + cs) / (1.0 + vb * cs));
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

HD static auto spherical_geometry_source_terms(prim_t p, double r0, double r1)
{
    // Eqn A8 in Zhang & MacFadyen (2006), integrated over the spherical shell
    // between r0 and r1, and specializing to radial velocity only.
    auto pg = p[2];
    auto dr2 = pow(r1, 2) - pow(r0, 2);
    auto srdot = pg * dr2;
    return cons_t{0.0, srdot, 0.0};
}




enum class Setup
{
    uniform,
    sod,
    mm96p1,
    wind,
    bmk,
    bomb,
    shell,
};
static Setup setup_from_string(const std::string& name)
{
    if (name == "uniform") return Setup::uniform;
    if (name == "sod") return Setup::sod;
    if (name == "mm96p1") return Setup::mm96p1;
    if (name == "wind") return Setup::wind;
    if (name == "bmk") return Setup::bmk;
    if (name == "bomb") return Setup::bomb;
    if (name == "shell") return Setup::shell;
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
    double bmk_gamma_shock = 5.0;
    double bomb_energy = 1.0E6;
    double bomb_rho_out = 1.0;
    double shell_u = 25.0;  // ejecta gamma-beta
    double shell_e = 1e-1;  // ejecta specific internal energy
    double shell_f = 100.0; // ejecta-to-ambient-medium comoving density ratio
    double shell_delta = 0.1;
    float dx = 1e-2;
    vec_t<float, 3> sod_l = {{1.0, 0.0, 1.000}};
    vec_t<float, 3> sod_r = {{0.1, 0.0, 0.125}};
    vec_t<float, 2> domain = {{1.0, 10.0}}; // x0, x1
    vec_t<float, 2> move = {{0.0, 0.0}}; // speed of mesh endpoints
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
    bomb_energy,
    bomb_rho_out,
    shell_u,
    shell_e,
    shell_f,
    shell_delta,
    sod_l,
    sod_r,
    domain,
    move,
    bc,
    sp,
    ts,
    outdir,
    method,
    setup,
    coords
);




struct planar_geometry_t
{
    planar_geometry_t(const Config& config, double t) : t(t)
    {
        x0 = config.domain[0];
        x1 = config.domain[1];
        v0 = config.move[0];
        v1 = config.move[1];
        ni = int((x1 - x0) / config.dx); // config.dx is essentially a hint
        dx = (x1 - x0) / ni;
    }
    HD double face_position(int i) const
    {
        return x0 + dx * i + face_velocity(i) * t;
    }
    HD double face_area(int i) const
    {
        return 1.0;
    }
    HD double face_velocity(int i) const
    {
        auto y = double(i) / ni;
        return y * v1 + (1.0 - y) * v0;
    }
    HD double cell_position(int i) const
    {
        return 0.5 * (face_position(i) + face_position(i + 1));
    }
    HD double cell_volume(int i) const
    {
        return face_position(i + 1) - face_position(i);
    }
    HD cons_t source_terms(prim_t p, double rm, double rp) const
    {
        return {};
    }
    HD index_space_t<1> cells_space() const
    {
        return range(ni).space();
    }
    double x0;
    double x1;
    double dx;
    double v0;
    double v1;
    double t;
    int ni;
};

struct spherical_geometry_t
{
    spherical_geometry_t(const Config& config, double t) : t(t)
    {
        r0 = config.domain[0];
        r1 = config.domain[1];
        v0 = config.move[0];
        v1 = config.move[1];
        ni = int((r1 - r0) / config.dx); // config.dx is essentially a hint
        dr = (r1 - r0) / ni;
    }
    HD double face_position(int i) const
    {
        return r0 + dr * i + face_velocity(i) * t;
    }
    HD double face_area(int i) const
    {
        auto r = face_position(i);
        return r * r;
    }
    HD double face_velocity(int i) const
    {
        auto y = double(i) / ni;
        return y * v1 + (1.0 - y) * v0;
    }
    HD double cell_position(int i) const
    {
        return 0.5 * (face_position(i) + face_position(i + 1));
    }
    HD double cell_volume(int i) const
    {
        auto r0 = face_position(i + 0);
        auto r1 = face_position(i + 1);
        return (r1 * r1 * r1 - r0 * r0 * r0) / 3.0;
    }
    HD cons_t source_terms(prim_t p, double rm, double rp) const
    {
        return spherical_geometry_source_terms(p, rm, rp);
    }
    HD index_space_t<1> cells_space() const
    {
        return range(ni).space();
    }
    double r0;
    double r1;
    double dr;
    double v0;
    double v1;
    double t;
    int ni;
};

struct grid_geometry_t
{
    grid_geometry_t(const Config& config, double t)
    : coords(coords_from_string(config.coords))
    , planar(config, t)
    , spherical(config, t) {
    }
    HD double face_position(int i) const
    {
        switch (coords) {
        case CoordinateSystem::planar: return planar.face_position(i);
        case CoordinateSystem::spherical: return spherical.face_position(i);
        }
    }
    HD double face_area(int i) const
    {
        switch (coords) {
        case CoordinateSystem::planar: return planar.face_area(i);
        case CoordinateSystem::spherical: return spherical.face_area(i);
        }
    }
    HD double face_velocity(int i) const
    {
        switch (coords) {
        case CoordinateSystem::planar: return planar.face_velocity(i);
        case CoordinateSystem::spherical: return spherical.face_velocity(i);
        }
    }
    HD double cell_position(int i) const
    {
        switch (coords) {
        case CoordinateSystem::planar: return planar.cell_position(i);
        case CoordinateSystem::spherical: return spherical.cell_position(i);
        }
    }
    HD double cell_volume(int i) const
    {
        switch (coords) {
        case CoordinateSystem::planar: return planar.cell_volume(i);
        case CoordinateSystem::spherical: return spherical.cell_volume(i);
        }
    }
    HD cons_t source_terms(prim_t p, double rm, double rp) const
    {
        switch (coords) {
        case CoordinateSystem::planar: return planar.source_terms(p, rm, rp);
        case CoordinateSystem::spherical: return spherical.source_terms(p, rm, rp);
        }
    }
    HD index_space_t<1> cells_space() const
    {
        switch (coords) {
        case CoordinateSystem::planar: return planar.cells_space();
        case CoordinateSystem::spherical: return spherical.cells_space();
        }
    }
    CoordinateSystem coords;
    planar_geometry_t planar;
    spherical_geometry_t spherical;
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
void update_prim(const State& state, G g, prim_array_t& p)
{
    auto u = state.cons;

    if (p.space() != u.space())
    {
        p = zeros<prim_t>(u.space()).cache();
    }
    p = range(u.space()).map([p, u, g] HD (int i)
    {
        auto dv = g.cell_volume(i);
        auto ui = u[i] / dv;
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

    if (bcl == 'f' && bcr == 'f')
    {
        return u.cache();
    }
    else if (bcl == 'f' && (bcr == 'o' || bcr == 'r'))
    {
        auto ur = u[u.size() - ng - 1];
        if (bcr == 'r') ur[1] *= -1.0;
        return u.at(ir).set(ur).cache();
    }
    else if ((bcl == 'o' || bcl == 'r') && bcr == 'f')
    {
        auto ul = u[ng];
        if (bcl == 'r') ul[1] *= -1.0;
        return u.at(il).set(ul).cache();
    }
    else if ((bcl == 'o' || bcl == 'r') && (bcr == 'o' || bcr == 'r'))
    {
        auto ul = u[ng];
        auto ur = u[u.size() - ng - 1];
        if (bcl == 'r') ul[1] *= -1.0;
        if (bcr == 'r') ur[1] *= -1.0;
        return u.at(il).set(ul).at(ir).set(ur).cache();
    }
    else
    {
        return u.cache();        
    }
}




template<class G>
static State next_pcm(const State& state, const Config& config, prim_array_t& p, double dt, int prim_dirty)
{
    auto u = state.cons;
    auto g = G(config, state.time);
    auto ic = range(u.space());
    auto iv = range(u.space().nudge(vec(0), vec(1)));
    auto rf = iv.map([g] HD (int i) { return g.face_position(i); });
    auto vf = iv.map([g] HD (int i) { return g.face_velocity(i); });
    auto da = iv.map([g] HD (int i) { return g.face_area(i); });
    auto dv = iv.map([g] HD (int i) { return g.cell_volume(i); });
    auto interior_cells = ic.space().contract(1);
    auto interior_faces = iv.space().contract(1);

    if (prim_dirty) {
        update_prim(state, g, p);
    }

    auto fhat = iv[interior_faces].map([p, u, dv, vf] HD (int i)
    {
        auto ul = u[i - 1] / dv[i];
        auto ur = u[i + 0] / dv[i];
        auto pl = p[i - 1];
        auto pr = p[i + 0];
        return riemann_hlle(pl, pr, ul, ur, vf[i]);
    }).cache();

    auto du = ic[interior_cells].map([rf, da, g, fhat, p] HD (int i)
    {
        auto rm = rf[i + 0];
        auto rp = rf[i + 1];
        auto am = da[i + 0];
        auto ap = da[i + 1];
        auto fm = fhat[i + 0];
        auto fp = fhat[i + 1];
        auto udot = g.source_terms(p[i], rm, rp);
        return fm * am - fp * ap + udot;
    }) * dt;

    return State{
        state.time + dt,
        state.iter + 1.0,
        set_bc(u.add(du), config, 1),
    };
}




template<class G>
static State next_plm(const State& state, const Config& config, prim_array_t& p, double dt, int prim_dirty)
{
    auto u = state.cons;
    auto g = G(config, state.time);
    auto ic = range(u.space());
    auto iv = range(u.space().nudge(vec(0), vec(1)));
    auto rf = iv.map([g] HD (int i) { return g.face_position(i); });
    auto vf = iv.map([g] HD (int i) { return g.face_velocity(i); });
    auto da = iv.map([g] HD (int i) { return g.face_area(i); });
    auto gradient_cells = ic.space().contract(1);
    auto interior_cells = ic.space().contract(2);
    auto interior_faces = iv.space().contract(2);
    auto plm_theta = config.theta;

    if (prim_dirty) {
        update_prim(state, g, p);
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

    auto fhat = iv[interior_faces].map([p, vf, grad] HD (int i)
    {
        auto pl = p[i - 1] + grad[i - 1] * 0.5;
        auto pr = p[i + 0] - grad[i + 0] * 0.5;
        auto ul = prim_to_cons(pl);
        auto ur = prim_to_cons(pr);
        return riemann_hlle(pl, pr, ul, ur, vf[i]);
    }).cache();

    auto du = ic[interior_cells].map([rf, da, g, fhat, p] HD (int i)
    {
        auto rm = rf[i + 0];
        auto rp = rf[i + 1];
        auto am = da[i + 0];
        auto ap = da[i + 1];
        auto fm = fhat[i + 0];
        auto fp = fhat[i + 1];
        auto udot = g.source_terms(p[i], rm, rp);
        return fm * am - fp * ap + udot;
    }) * dt;

    return State{
        state.time + dt,
        state.iter + 1.0,
        set_bc(u.add(du), config, 2),
    };
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

    auto dx = config.dx;
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
        auto setup = setup_from_string(config.setup);
        auto bmk_gamma_shock = config.bmk_gamma_shock;
        auto bomb_energy = config.bomb_energy;
        auto sod_l = cast<double>(config.sod_l);
        auto sod_r = cast<double>(config.sod_r);
        auto bomb_rho_out = config.bomb_rho_out;
        auto shell_u = config.shell_u;
        auto shell_e = config.shell_e;
        auto shell_f = config.shell_f;
        auto shell_delta = config.shell_delta;
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
                    // Sod shocktube -- or any Riemann problem using sod_l / sod_r
                    //
                    if (x < 0.5) { // x0 + 0.5 * (x1 - x0)) {
                        return sod_l;
                    } else {
                        return sod_r;
                    }
                }
            case Setup::mm96p1: {
                    // Problem 1 from Marti & Muller 1996
                    //
                    if (x < 0.5) {
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
                    auto pre = 1e-10 * pow(rho, gamma_law); // uniform specfic entropy
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
            case Setup::bomb: {
                auto r_in = 1.0;
                auto bomb_mass = 1.0;
                auto bomb_volume = 4.0 / 3.0 * M_PI * pow(r_in, 3.0);
                auto bomb_rho_in = bomb_mass / bomb_volume;
                auto p_in = bomb_energy / bomb_volume * (gamma_law - 1.0);
                auto p_out = bomb_rho_out * 1e-6;
                if (x < r_in) {
                    return vec(bomb_rho_in, 0.0, p_in);
                } else {
                    return vec(bomb_rho_out, 0.0, p_out);
                }
            }
            case Setup::shell: {
                auto rho_out = 1.0;
                auto p_out = rho_out * 1e-6;
                auto r_in = 1.0;
                auto rho_in = rho_out * shell_f * exp(-pow(x - r_in, 2.) / pow(shell_delta, 2.) / 2.);
                auto u_in = shell_u * exp(-pow(x - r_in, 2.) / pow(shell_delta, 2.) / 2.);
                auto p_in = rho_in * shell_e * (gamma_law - 1.);
                if (x < r_in) {
                    return vec(rho_in, u_in, p_in);
                } else {
                    return vec(rho_out, 0.0, p_out);
                }
            }
            default: return {};
            }
        };
        auto g = grid_geometry_t(config, 0.0);
        auto ic = range(g.cells_space());
        auto u = ic.map([=] (int i) {
            auto xc = g.cell_position(i);
            auto dv = g.cell_volume(i);
            auto p = initial_primitive(xc);
            auto u = prim_to_cons(p);
            return u * dv;
        });
        state.time = 0.0;
        state.iter = 0.0;
        state.cons = u.cache();
    }
    void update(State& state) const override
    {
        switch (coords_from_string(config.coords))
        {
        case CoordinateSystem::planar:
            return update_state<planar_geometry_t>(state, config);
        case CoordinateSystem::spherical:
            return update_state<spherical_geometry_t>(state, config);
        }
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
        auto g = grid_geometry_t(config, state.time);
        auto ic = range(g.cells_space());
        auto xc = ic.map([g] HD (int i) { return g.cell_position(i); });
        auto dv = ic.map([g] HD (int i) { return g.cell_volume(i); });
        auto cons_field = [] (uint n) {
            return [n] HD (cons_t u) {
                return cons_to_prim(u).get()[n];
            };
        };
        switch (column) {
        case 0: return (state.cons / dv).map(cons_field(0)).cache();
        case 1: return (state.cons / dv).map(cons_field(1)).cache();
        case 2: return (state.cons / dv).map(cons_field(2)).cache();
        case 3: return xc.cache();
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




#ifndef VAPOR_GUI

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

#else

/**
================================================================================

Approach to the UI
------------------

Note: the implementation below assumes MacOS + ImGui + SDL2 + Metal. Other
platforms and ImGui backends could easily be supported.

The UI loop passes commands to the simulation process. The commands include:

- restart the simulation
- start the simulation
- step the simulation
- new config instances

The simulation passes back to the UI a sequence of status objects. The
status objects include the simulation state, an iteration message, the
current configuration struct, and a map of science products.
================================================================================
*/

#import <Metal/Metal.h>
#import <QuartzCore/QuartzCore.h>
#include <tuple>
#include <variant>
#include <SDL.h>
#include "imgui.h"
#include "imgui_impl_metal.h"
#include "imgui_impl_sdl2.h"
#include "implot.h"
#include "imfilebrowser.h"

enum struct Action
{
    nothing,
    restart,
    step,
    quit,
};
struct Status {
    Config config;
    State state;
    vec_t<char, 256> message;
    std::map<std::string, Product> products;
};
using Command = std::variant<Action, Config>;




struct SimulationProcess
{
    SimulationProcess()
    {
        sim.initial_state(state);
    }
    bool operator()(Action action)
    {
        switch (action) {
        case Action::nothing:
            return true;
        case Action::restart:
            sim.initial_state(state);
            return true;
        case Action::step:
            secs = time_call(sim.updates_per_batch(), [&] { sim.update(state); });
            return true;
        case Action::quit:
            return false;
        }
    }
    bool operator()(Config config)
    {
        if (sim.get_config().dx != config.dx ||
            sim.get_config().domain != config.domain ||
            sim.get_config().setup != config.setup) {
            sim.get_config() = config;
            sim.initial_state(state);
        }
        else {
            sim.get_config() = config;            
        }
        return true;
    }
    bool operator()(Command command)
    {
        return std::visit(*this, command);
    }
    Status status() const
    {
        auto status = Status();
        status.config = sim.get_config();
        status.state = state;
        status.message = sim.status_message(state, secs);
        int col = 0;
        while (auto name = sim.get_product_name(col)) {
            status.products[name] = sim.compute_product(state, col);
            ++col;
        }
        return status;
    }
    Blast sim;
    State state;
    double secs = 1.0;
};

// ===> If putting the simulation on a background thread <===
// auto responder(std::function<Command(Status)> next)
// {
//     return [next] {
//         auto p = SimulationProcess();
//         while (p(next(p.status()))) {
//         }
//     };
// }




class App
{
public:
    App()
    {
        if (startup()) {
            printf("startup error\n");
            exit(1);
        }
        browser = ImGui::FileBrowser(ImGuiFileBrowserFlags_CloseOnEsc);
        browser.SetTitle("File Browser");
        browser.SetTypeFilters({".cfg"});
    }
    ~App()
    {
        shutdown();
    }
    Command draw(const Status& status)
    {
        Command command = Action::nothing;
        static bool show_config = false;
        static bool show_style = false;
        static bool auto_step = false;
        static bool draw_markers = true;

        ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f));
        ImGui::SetNextWindowSize(ImGui::GetIO().DisplaySize);
        ImGui::Begin("Window", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoBringToFrontOnFocus);

        if (ImGui::Button("Quit")) {
            auto_step = false;
            command = Action::quit;
        }
        ImGui::SameLine();
        if (ImGui::Button("Config")) {
            show_config = true;
        }
        ImGui::SameLine();
        if (ImGui::Button("Load")) {
            browser.Open();
        }
        ImGui::SameLine();
        if (ImGui::Button("Restart")) {
            auto_step = false;
            command = Action::restart;
        }
        ImGui::SameLine();
        ImGui::BeginDisabled(auto_step);
        if (ImGui::Button("Step")) {
            command = Action::step;
        }
        ImGui::EndDisabled();
        ImGui::SameLine();
        ImGui::Checkbox("Auto", &auto_step);
        if (auto_step) {
            command = Action::step;
        }
        ImGui::SameLine();
        ImGui::Text("%s", status.message.data);

        if (ImGui::Button("Style")) {
            show_style = true;
        }

        if (ImPlot::BeginPlot("##blast", ImVec2(-1.0, -1.0)))
        {
            auto x = status.products.at("cell_coordinate");
            if (true) {
                auto x0 = min(x);
                auto x1 = max(x);
                x = ((x - x0) / (x1 - x0)).cache();
            }
            for (const auto& [name, y] : status.products) {
                if (name != "cell_coordinate") {
                    if (draw_markers) {
                        ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
                        ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.25f);
                    }
                    ImPlot::PlotLine(name.data(), x.data(), y.data(), x.size());
                }
            }
            ImPlot::EndPlot();
        }

        ImGui::End();
        browser.Display();

        if (browser.HasSelected())
        {
            auto config = status.config;
            auto fname = browser.GetSelected().string();
            set_from_key_vals(config, readfile(fname.data()).data());
            browser.ClearSelected();
            command = config;
        }

        // const char* rk_options[] = {"None", "RK1", "RK2", "RK3"};
        // const char* method_options[] = {"pcm", "plm"};
        // static int method_index = 0;

        if (show_style) {
            ImGui::SetNextWindowSize(ImVec2(300.0, 0.0));
            ImGui::Begin("Style", &show_style, ImGuiWindowFlags_NoResize);
            ImGui::Checkbox("Markers", &draw_markers);
            ImGui::End();
        }
        if (show_config) {
            auto config = status.config;

            ImGui::SetNextWindowSize(ImVec2(300.0, 0.0));
            ImGui::Begin("Config", &show_config, ImGuiWindowFlags_NoResize);

            if (ImGui::SliderFloat("dx", &config.dx, 1e-1f, 1e-5f, "%.6g", ImGuiSliderFlags_Logarithmic)) {
                command = config;
            }
            if (ImGui::DragFloatRange2("domain", &config.domain[0], &config.domain[1], 0.05f)) {
                command = config;
            }
            // if (ImGui::Combo("rk", &config.rk, rk_options, IM_ARRAYSIZE(rk_options))) {
            // }
            // if (ImGui::Combo("method", &method_index, method_options, IM_ARRAYSIZE(method_options))) {
            //     config.method = method_options[method_index];
            // }
            ImGui::End();
        }
        return command;
    }
    int run(int argc, const char *argv[])
    {
        bool done = false;
        // ===> If putting the simulation on a background thread <===
        //
        // auto m = std::mutex();
        // auto queue = std::queue<Command>();
        // auto last_status = Status();
        // auto sim_thread = std::thread(responder([&m, &queue, &last_status] (Status status) -> Command
        // {
        //     {
        //         auto lock = std::lock_guard<std::mutex>(m);
        //         last_status = status;
        //     }
        //     while (true) {
        //         auto lock = std::lock_guard<std::mutex>(m);
        //         if (! queue.empty()) {
        //             auto command = queue.front();
        //             queue.pop();
        //             return command;
        //         }
        //     }
        //     assert(false);
        // }));
        //
        // ===> Else <===
        auto sim_process = SimulationProcess();
        {
            auto config = Config();
            for (int n = 1; n < argc; ++n)
            {
                if (argv[n][0] == '-') {
                    // do nothing
                } else if (strstr(argv[n], ".cfg")) {
                    set_from_key_vals(config, readfile(argv[n]).data());
                } else {
                    set_from_key_vals(config, argv[n]);
                }
            }
            // ===> If putting the simulation on a background thread <===
            // auto lock = std::lock_guard<std::mutex>(m);
            // queue.push(config);
            //
            // ===> Else <===
            sim_process(config);
        }
        while (! done)
        {
            SDL_Event event;
            while (SDL_PollEvent(&event))
            {
                ImGui_ImplSDL2_ProcessEvent(&event);
                if ((event.type == SDL_QUIT) ||
                    (event.type == SDL_WINDOWEVENT &&
                     event.window.event == SDL_WINDOWEVENT_CLOSE &&
                     event.window.windowID == SDL_GetWindowID(window)))
                {
                    // ===> If putting the simulation on a background thread <===
                    // auto guard = std::lock_guard<std::mutex>(m);
                    // queue.push(Action::quit);
                    // ===> Else <===
                    sim_process(Action::quit);
                    done = true;
                }
            }
            @autoreleasepool {
                auto d = new_frame();
                // ===> If putting the simulation on a background thread <===
                // auto guard = std::lock_guard<std::mutex>(m);
                // auto command = draw(last_status);
                // if (! queue.empty() &&
                //     std::holds_alternative<Action>(command) &&
                //     std::holds_alternative<Action>(queue.back()) &&
                //     std::get<Action>(command) == std::get<Action>(queue.back())) {
                //     // skip this command, it's already in the queue
                // } else {
                //     queue.push(command);
                // }
                // if (std::holds_alternative<Action>(command) && std::get<Action>(command) == Action::quit) {
                //     done = true;
                // }
                if (! sim_process(draw(sim_process.status()))) {
                    done = true;
                }
                // ImGui::ShowDemoWindow();
                // ImPlot::ShowDemoWindow();
                end_frame(d);
            }
        }
        // ===> If putting the simulation on a background thread <===
        // sim_thread.join();
        return 0;
    }
private:
    std::tuple<id<CAMetalDrawable>, id<MTLCommandBuffer>, id <MTLRenderCommandEncoder>> new_frame()
    {
        float clear_color[4] = {0.f, 0.f, 0.f, 1.f};
        int width, height;
        SDL_GetRendererOutputSize(renderer, &width, &height);
        layer.drawableSize = CGSizeMake(width, height);
        id<CAMetalDrawable> drawable = [layer nextDrawable];
        id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
        renderPassDescriptor.colorAttachments[0].clearColor = MTLClearColorMake(clear_color[0] * clear_color[3], clear_color[1] * clear_color[3], clear_color[2] * clear_color[3], clear_color[3]);
        renderPassDescriptor.colorAttachments[0].texture = drawable.texture;
        renderPassDescriptor.colorAttachments[0].loadAction = MTLLoadActionClear;
        renderPassDescriptor.colorAttachments[0].storeAction = MTLStoreActionStore;
        id <MTLRenderCommandEncoder> renderEncoder = [commandBuffer renderCommandEncoderWithDescriptor:renderPassDescriptor];
        [renderEncoder pushDebugGroup:@"ImGui demo"];
        ImGui_ImplMetal_NewFrame(renderPassDescriptor);
        ImGui_ImplSDL2_NewFrame();
        ImGui::NewFrame();
        return std::make_tuple(drawable, commandBuffer, renderEncoder);
    }
    void end_frame(std::tuple<id<CAMetalDrawable>, id<MTLCommandBuffer>, id <MTLRenderCommandEncoder>> t)
    {
        auto [drawable, commandBuffer, renderEncoder] = t;
        ImGui::Render();
        ImGui_ImplMetal_RenderDrawData(ImGui::GetDrawData(), commandBuffer, renderEncoder);
        [renderEncoder popDebugGroup];
        [renderEncoder endEncoding];
        [commandBuffer presentDrawable:drawable];
        [commandBuffer commit];
    }
    int startup()
    {
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImPlot::CreateContext();
        ImGuiIO& io = ImGui::GetIO();
        io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
        io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;
        ImGui::StyleColorsDark();
        // WARNING: font loading assumes pwd is the project root; should be fixed
        // io.Fonts->AddFontFromFileTTF("imgui/misc/fonts/DroidSans.ttf", 16.0f);
        // io.Fonts->AddFontFromFileTTF("imgui/misc/fonts/Roboto-Medium.ttf", 16.0f);
        // io.Fonts->AddFontFromFileTTF("imgui/misc/fonts/Cousine-Regular.ttf", 15.0f);

        if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) != 0)
        {
            printf("Error: %s\n", SDL_GetError());
            return -1;
        }
        SDL_SetHint(SDL_HINT_RENDER_DRIVER, "metal");
        SDL_SetHint(SDL_HINT_IME_SHOW_UI, "1");

        window = SDL_CreateWindow("Vapor Viewer",
            SDL_WINDOWPOS_CENTERED,
            SDL_WINDOWPOS_CENTERED,
            1280, 720,
            SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);

        if (window == nullptr)
        {
            printf("Error creating window: %s\n", SDL_GetError());
            return -2;
        }

        renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
        if (renderer == nullptr)
        {
            printf("Error creating renderer: %s\n", SDL_GetError());
            return -3;
        }

        layer = (__bridge CAMetalLayer*)SDL_RenderGetMetalLayer(renderer);
        layer.pixelFormat = MTLPixelFormatBGRA8Unorm;
        ImGui_ImplMetal_Init(layer.device);
        ImGui_ImplSDL2_InitForMetal(window);

        commandQueue = [layer.device newCommandQueue];
        renderPassDescriptor = [MTLRenderPassDescriptor new];

        return 0;
    }
    void shutdown()
    {
        ImGui_ImplMetal_Shutdown();
        ImGui_ImplSDL2_Shutdown();
        ImPlot::DestroyContext();
        ImGui::DestroyContext();
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
    }

    SDL_Window *window;
    SDL_Renderer* renderer;
    CAMetalLayer *layer;
    id<MTLCommandQueue> commandQueue;
    MTLRenderPassDescriptor* renderPassDescriptor;
    ImGui::FileBrowser browser;
};




int main(int argc, const char **argv)
{
    return App().run(argc, argv);
}

#endif
