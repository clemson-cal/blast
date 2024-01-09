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




// #define VAPOR_USE_SHARED_PTR_ALLOCATOR
#include "vapor/vapor.hpp"
using namespace vapor;




/**
 * 
 */
using cons_t = vec_t<double, 3>;
using prim_t = vec_t<double, 3>;
using cons_array_t = memory_backed_array_t<1, cons_t, ref_counted_ptr_t>;
using prim_array_t = memory_backed_array_t<1, prim_t, ref_counted_ptr_t>;
using Product = memory_backed_array_t<1, double, ref_counted_ptr_t>;




/**
 * 
 */
#define index_density 0
#define index_pressure 2
#define index_energy 2
#define gamma (4.0 / 3.0)
#define min2(a, b) ((a) < (b) ? (a) : (b))
#define max2(a, b) ((a) > (b) ? (a) : (b))
#define min3(a, b, c) min2(a, min2(b, c))
#define max3(a, b, c) max2(a, max2(b, c))




/**
 * 
 */
HD auto gamma_beta_squared(prim_t p) -> double
{
    return p[1] * p[1];
}

HD auto momentum_squared(cons_t u) -> double
{
    return u[1] * u[1];
}

HD auto lorentz_factor(prim_t p) -> double
{
    return sqrt(1.0 + gamma_beta_squared(p));
}

HD auto beta_component(prim_t p, int axis) -> double
{
    return p[axis + 1] / lorentz_factor(p);
}

HD auto enthalpy_density(prim_t p) -> double
{
    auto rho = p[index_density];
    auto pre = p[index_pressure];
    return rho + pre * (1.0 + 1.0 / (gamma - 1.0));
}

HD auto prim_to_cons(prim_t p) -> cons_t
{
    auto rho = p[index_density];
    auto pre = p[index_pressure];
    auto w = lorentz_factor(p);
    auto h = enthalpy_density(p) / rho;
    auto m = rho * w;
    auto u = cons_t{};
    u[0] = m;
    u[1] = m * (h * p[1]);
    u[2] = m * (h * w - 1.0) - pre;
    return u;
}

HD auto cons_to_prim(cons_t cons, double p=0.0) -> prim_t
{
    auto newton_iter_max = 50;
    auto error_tolerance = 1e-12 * (cons[index_density] + cons[index_energy]);
    auto gm = gamma;
    auto m = cons[index_density];
    auto tau = cons[index_energy];
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
            printf("c2p failed; D=%f tau=%f\n", cons[index_density], cons[index_energy]);
            exit(1);
        }
        if (fabs(f) < error_tolerance) {
            w0 = w;
            break;
        }        
        p -= f / g;
        n += 1;
    }
    return {m / w0, w0 * cons[1] / (tau + m + p), p};
}

HD auto prim_and_cons_to_flux(prim_t p, cons_t u, int axis) -> cons_t
{
    double pre = p[index_pressure];
    double vn = beta_component(p, axis);
    auto f = cons_t{};
    f[0] = vn * u[0];
    f[1] = vn * u[1] + pre * (axis == 0);
    f[2] = vn * u[2] + pre * vn;
    return f;
}

HD auto sound_speed_squared(prim_t p) -> double
{
    const double pre = p[index_pressure];
    const double rho_h = enthalpy_density(p);
    return gamma * pre / rho_h;
}

HD auto outer_wavespeeds(prim_t p, int axis) -> dvec_t<2>
{
    double a2 = sound_speed_squared(p);
    double uu = gamma_beta_squared(p);
    double vn = beta_component(p, axis);
    double vv = uu / (1.0 + uu);
    double v2 = vn * vn;
    double k0 = sqrt(a2 * (1.0 - vv) * (1.0 - vv * a2 - v2 * (1.0 - a2)));
    return vec(
        (vn * (1.0 - a2) - k0) / (1.0 - vv * a2),
        (vn * (1.0 - a2) + k0) / (1.0 - vv * a2)
    );
}

HD auto riemann_hlle(prim_t pl, prim_t pr, cons_t ul, cons_t ur) -> cons_t
{
    auto fl = prim_and_cons_to_flux(pl, ul, 0);
    auto fr = prim_and_cons_to_flux(pr, ur, 0);
    auto al = outer_wavespeeds(pl, 0);
    auto ar = outer_wavespeeds(pr, 0);
    auto alm = al[0];
    auto alp = al[1];
    auto arm = ar[0];
    auto arp = ar[1];
    auto am = min3(alm, arm, 0.0);
    auto ap = max3(alp, arp, 0.0);
    return (fl * ap - fr * am - (ul - ur) * ap * am) / (ap - am);
};




/**
 * 
 */
struct Config
{
    int num_zones = 100;
    int fold = 50;
    int rk = 1;
    double tfinal = 0.0;
    double cpi = 0.0;
    double spi = 0.0;
    double tsi = 0.0;
    std::vector<uint> sp;
    std::vector<uint> ts;
    std::string outdir = ".";
};
VISITABLE_STRUCT(Config, num_zones, fold, rk, tfinal, cpi, spi, tsi, sp, ts, outdir);




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




State average(const State& a, const State& b, double x)
{
    return x == 1.0 ? a : State{
        (a.time * (1.0 - x) + b.time * x),
        (a.iter * (1.0 - x) + b.iter * x),
        (a.cons * (1.0 - x) + b.cons * x).cache()
    };
}




/**
 * 
 */
using DiagnosticData = prim_array_t;




State next(const State& state, const Config& config, double dt)
{
    static prim_array_t primitive;

    auto ni = config.num_zones;
    auto dx = 1.0 / ni;
    auto iv = range(ni + 1);
    auto ic = range(ni);
    auto u = state.cons;
    auto interior_faces = iv.space().contract(1);
    auto interior_cells = ic.space().contract(1);

    if (primitive.space() != u.space())
    {
        primitive = zeros<prim_t>(u.space()).cache();
    }

    primitive = ic.map([p=primitive, u] HD (int i)
    {
        auto ui = u[i];
        auto pi = p[i];
        return cons_to_prim(ui, pi[2]);
    }).cache();

    auto fhat = iv[interior_faces].map([p=primitive, u] HD (int i)
    {
        auto ul = u[i - 1];
        auto ur = u[i];
        auto pl = p[i - 1];
        auto pr = p[i];
        return riemann_hlle(pl, pr, ul, ur);
    }).cache();

    auto du = ic[interior_cells].map([fhat, dt, dx] HD (int i)
    {
        auto fm = fhat[i];
        auto fp = fhat[i + 1];
        return (fp - fm) * (-dt / dx);
    });

    return State{
        state.time + dt,
        state.iter + 1.0,
        (u.at(interior_cells) + du).cache(),
    };
}




auto cell_coordinates(const Config& config)
{
    auto ni = config.num_zones;
    auto dx = 1.0 / ni;
    auto ic = range(ni);
    auto xc = (ic + 0.5) * dx;
    return xc;
}




void update_state(State& state, const Config& config)
{
    auto ni = config.num_zones;
    auto dx = 1.0 / ni;
    auto dt = dx * 0.3;
    auto s0 = state;

    switch (config.rk)
    {
        case 1: {
            state = next(s0, config, dt);
            break;
        }
        case 2: {
            auto s1 = average(s0, next(s0, config, dt), 1./1);
            auto s2 = average(s0, next(s1, config, dt), 1./2);
            state = s2;
            break;
        }
        case 3: {
            auto s1 = average(s0, next(s0, config, dt), 1./1);
            auto s2 = average(s0, next(s1, config, dt), 1./4);
            auto s3 = average(s0, next(s2, config, dt), 2./3);
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
        auto initial_conserved = [] HD (double x)
        {
            if (x < 0.5)
                return prim_to_cons(vec(1.0, 0.0, 1.0));
            else
                return prim_to_cons(vec(0.1, 0.0, 0.125));
        };
        state.time = 0.0;
        state.iter = 0.0;
        state.cons = cell_coordinates(config).map(initial_conserved).cache();
    }
    void update(State& state) const override
    {
        update_state(state, config);
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
        switch (column) {
        case 0: return state.cons.map([] HD (cons_t u) { return cons_to_prim(u)[0]; }).cache();
        case 1: return state.cons.map([] HD (cons_t u) { return cons_to_prim(u)[1]; }).cache();
        case 2: return state.cons.map([] HD (cons_t u) { return cons_to_prim(u)[2]; }).cache();
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
            1e-6 * config.num_zones / secs_per_update);
    }
};




int main(int argc, const char **argv)
{
    try {
        return Blast().run(argc, argv);
    }
    catch (const std::exception& e) {
        vapor::print(e.what());
        vapor::print("\n");
    }
    return 0;
}
