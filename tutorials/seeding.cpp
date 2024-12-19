/** TRACCC tutorial for beginners
 *
 * (c) 2024 UC Berkeley and LBNL
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/edm/spacepoint.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"

// Detray include(s).
#include "detray/navigation/detail/helix.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

using namespace traccc;

int main()
{
    // 2T magnetic field in the z-axis
    const vector3 B{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                    2.f * unit<scalar>::T};

    /***********************************************
     * Generate spacepoint data with helices
     ***********************************************/

    // Make two helices of negative charged particle
    const scalar q{-1.f * unit<scalar>::e};
    const point3 pos{0.f, 0.f, 0.f};
    const scalar time{0.f};

    // The momentum of the 1st helix [GeV/c]: (1, 0, 1)
    const vector3 mom_1st{1.f * unit<scalar>::GeV, 0.f, 1.f * unit<scalar>::GeV};
    detray::detail::helix<traccc::default_algebra> hlx_1st(
        pos, time, vector::normalize(mom_1st), q / vector::norm(mom_1st), &B);

    // The momentum of the 2nd helix [GeV/c]: (0, 3, 4)
    const vector3 mom_2nd{0.f, 3.f * unit<scalar>::GeV, 4.f * unit<scalar>::GeV};
    detray::detail::helix<traccc::default_algebra> hlx_2nd(
        pos, time, vector::normalize(mom_2nd), q / vector::norm(mom_2nd), &B);

    // Make three spacepoints from each helix
    // Spacepoint constructor inputs:
    // 1. 3D position
    // 2. corresponding measurement object (optional)
    spacepoint_collection_types::host spacepoints;
    spacepoints.push_back({hlx_1st(50 * unit<scalar>::mm), {}});
    spacepoints.push_back({hlx_1st(100 * unit<scalar>::mm), {}});
    spacepoints.push_back({hlx_1st(150 * unit<scalar>::mm), {}});
    spacepoints.push_back({hlx_2nd(50 * unit<scalar>::mm), {}});
    spacepoints.push_back({hlx_2nd(100 * unit<scalar>::mm), {}});
    spacepoints.push_back({hlx_2nd(150 * unit<scalar>::mm), {}});

    /***********************************************
     * Run Seeding
     ***********************************************/

    // Config objects
    traccc::seedfinder_config finder_config;
    finder_config.bFieldInZ = B[2]; // Set the B field assumption for seeding

    traccc::spacepoint_grid_config grid_config(finder_config);
    traccc::seedfilter_config filter_config;

    // Seeding algorithm object
    vecmem::host_memory_resource host_mr;
    traccc::seeding_algorithm sa(finder_config, grid_config, filter_config,
                                 host_mr);

    // Run seeding
    auto seeds = sa(spacepoints);

    std::cout << "Number of seeds found: " << seeds.size() << std::endl;

    /***********************************************
     * Run Track Parameter Estimation
     ***********************************************/

    // Track parameter estimation algorithm object
    traccc::track_params_estimation tp(host_mr);

    // Run track parameter estimation
    auto bound_params = tp(spacepoints, seeds, B);

    std::cout << "Momentum of the 1st seed [GeV/c]: " << bound_params[0].p(q) << std::endl;
    std::cout << "Momentum of the 2nd seed [GeV/c]: " << bound_params[1].p(q) << std::endl;

    /*************************************************************
     * Tasks
     * 
     * 1. Create the 3rd helix with a positive charge 
     * 2. Generate new three spacepoints from the 3rd helix
     * 3. Run the seeding and track parameter estimation again 
     * 4. Check the results (momentum and position)
     *************************************************************/

    return 1;
}