#include "vapor/vapor.hpp"
using namespace vapor;

struct cons_t: public vec_t<double, 3> {};
struct prim_t: public vec_t<double, 3> {};

template<> struct hdf5_repr<cons_t> : public hdf5_repr<vec_t<double, 3>> {};
template<> struct hdf5_repr<prim_t> : public hdf5_repr<vec_t<double, 3>> {};




struct Config
{
    int num_zones = 100;
    double tfinal = 1.0;
    double cpi = 0.0;
    std::vector<int> ts;
};
VISITABLE_STRUCT(Config, num_zones, tfinal, cpi, ts);




struct State
{
    double time;
    int iteration;
    memory_backed_array_t<1, cons_t, ref_counted_ptr_t> conserved;
};
VISITABLE_STRUCT(State, time, iteration, conserved);




using DiagnosticData = memory_backed_array_t<1, cons_t, ref_counted_ptr_t>;




class ShockSimulation : public Simulation<Config, State, DiagnosticData>
{
public:
    double get_time(const State& state) const override
    {
        return state.time;
    }
    virtual uint get_iteration(const State& state) const override
    {
        return state.iteration;
    }
    void initial_state(State& state) const override
    {
        state.time = 0.0;
        state.iteration = 0;
        state.conserved = zeros<cons_t>(uvec(config.num_zones)).cache();
    }
    void update(State& state) const override
    {
        state.time += 0.0085216;
        state.iteration += 1;
    }
    bool should_continue(const State& state) const override
    {
        return state.time < config.tfinal;
    }
    double checkpoint_interval() const override
    { 
        return config.cpi;
    }
};




int main(int argc, const char **argv)
{
    try {
        auto sim = ShockSimulation();
        return run(argc, argv, sim);
    }
    catch (const std::exception& e) {
        printf("[error]: %s\n", e.what());
    }
    return 0;
}
