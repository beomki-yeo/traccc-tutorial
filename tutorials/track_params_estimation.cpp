/** TRACCC tutorial
 *
 * (c) 2024 UC Berkeley and LBNL
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/edm/spacepoint.hpp"
#include "traccc/seeding/track_params_estimation.hpp"

// Detray include(s).
#include "detray/navigation/detail/helix.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

using namespace traccc;

int main()
{
    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;

    // 2 T B field in z-axis
    const vector3 B{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                    2.f * unit<scalar>::T};

    // Make a helix of negative charged particle with a momentum of (1 GeV/c, 0, 1 GeV/c)
    const scalar q{-1.f * unit<scalar>::e};
    const point3 pos{0.f, 0.f, 0.f};
    const scalar time{0.f};
    const vector3 mom{1.f * unit<scalar>::GeV, 0.f, 1.f * unit<scalar>::GeV};
    detray::detail::helix<traccc::default_algebra> hlx(
        pos, time, vector::normalize(mom), q / vector::norm(mom), &B);

    // Make three spacepoints with the helix
    spacepoint_collection_types::host spacepoints;
    spacepoints.push_back({hlx(50 * unit<scalar>::mm), {}});
    spacepoints.push_back({hlx(100 * unit<scalar>::mm), {}});
    spacepoints.push_back({hlx(150 * unit<scalar>::mm), {}});

    // Make a seed from the three spacepoints
    seed_collection_types::host seeds;
    seeds.push_back({0u, 1u, 2u, 0.f, 0.f});

    // Track parameter estimation algorithm object
    traccc::track_params_estimation tp(host_mr);

    // Run track parameter estimation
    auto bound_params = tp(spacepoints, seeds, B);

    // Print the momentum of esitmated track parameter
    const vector3 p = bound_params[0].mom(q);

    std::cout << "Truth momentum [GeV/c]    : (" << mom[0] << "," << mom[1]
              << "," << mom[2] << ")" << std::endl;

    std::cout << "Estimated momentum [GeV/c]: (" << p[0] << "," << p[1]
              << "," << p[2] << ")" << std::endl;
}