/** TRACCC tutorial for beginners
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Traccc include(s).
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/geometry/silicon_detector_description.hpp"
#include "traccc/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"

// Detray include(s).
#include "detray/test/utils/detectors/build_telescope_detector.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

using namespace traccc;

int main()
{
    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;

    // Prepare input cell data (Two clusters)
    /*
     * [ ] [ ] [ ] [ ] [ ] [ ] [ ] [ ]
     * [ ] [ ] [ ] [ ] [ ] [ ] [ ] [ ]
     * [ ] [X] [X] [X] [X] [ ] [ ] [ ]
     * [ ] [ ] [ ] [ ] [ ] [ ] [ ] [ ]
     * [ ] [ ] [ ] [ ] [ ] [ ] [X] [ ]
     * [ ] [ ] [ ] [ ] [ ] [X] [X] [X]
     * [ ] [ ] [ ] [ ] [ ] [ ] [X] [ ]
     * [ ] [ ] [ ] [ ] [ ] [ ] [ ] [ ]
     */
    traccc::edm::silicon_cell_collection::host cells{host_mr};
    cells.resize(9u);
    cells.channel0() = {1u, 2u, 3u, 4u, 6u, 5u, 6u, 7u, 6u};
    cells.channel1() = {2u, 2u, 2u, 2u, 4u, 5u, 5u, 5u, 6u};
    cells.activation() = {1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f};
    cells.time() = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    cells.module_index() = {0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u};

    // Make a tracking geometry with a single pixel module
    // 1. The normal axis of the module plane is the x-axis
    // 2. Its origin is at (10 mm, 0 mm, 0 mm) 
    const vector3 align_axis{1.f, 0.f, 0.f}; // normal axis in x-axis
    detray::detail::ray<traccc::default_algebra> pilot_track{
        {0, 0, 0}, 0, align_axis, -1};
    detray::mask<detray::rectangle2D> rect{0u,
                                           100.f * traccc::unit<scalar>::mm,
                                           100.f * traccc::unit<scalar>::mm};
    detray::tel_det_config tel_cfg{rect, pilot_track};
    tel_cfg.positions({10 * traccc::unit<scalar>::mm}); // The origin is at (10 mm, 0, 0)

    const auto [det, name_map] = build_telescope_detector(host_mr, tel_cfg);

    // Define the digitization information of the module
    traccc::silicon_detector_description::host dd{host_mr};
    dd.resize(1u);
    dd.reference_x()[0] = 0.f;
    dd.reference_y()[0] = 0.f;
    dd.pitch_x()[0] = 1.f; // Pitch of cells in the local x-axis
    dd.pitch_y()[0] = 1.f; // Pitch of cells in the local y-axis
    dd.dimensions()[0] = 2; // 2 dimensions (i.e. pixel)
    dd.geometry_id()[0] = detray::geometry::barcode{0u};

    /**************************************************
     * Run the clusterization and create measurements
     **************************************************/

    traccc::host::clusterization_algorithm ca(host_mr);
    auto measurements = ca(vecmem::get_data(cells), vecmem::get_data(dd));

    std::cout << "Number of measurements: " << measurements.size() << std::endl;
    std::cout << "Local position of the 1st measurement: ("
              << measurements[0].local[0] << "," << measurements[0].local[1] << ")"
              << std::endl;
    std::cout << "Local position of the 2nd measurement: ("
              << measurements[1].local[0] << "," << measurements[1].local[1] << ")"
              << std::endl;

    /*****************************************************************
     * Run the spacepoint formation (local-to-global transformation)
     *****************************************************************/

    traccc::host::silicon_pixel_spacepoint_formation_algorithm sf(host_mr);
    auto spacepoints = sf(det, vecmem::get_data(measurements));

    std::cout << "Number of spacepoints: " << spacepoints.size() << std::endl;
    std::cout << "Global position of the 1st spacepoint: ("
              << spacepoints[0].global[0] << "," 
              << spacepoints[0].global[1] << "," 
              << spacepoints[0].global[2] << ")" << std::endl;
    std::cout << "Global position of the 2nd spacepoint: ("
              << spacepoints[1].global[0] << "," 
              << spacepoints[1].global[1] << "," 
              << spacepoints[1].global[2] << ")" << std::endl;              

    /*************************************************************
     * Practice
     * 
     * 1. Make the 3rd cluster by adding cells
     * 2. Run the seeding and track parameter estimation again 
     * 4. Check the results (local and global position)
     *************************************************************/

    return 1;
}