/** TRACCC tutorial for beginners
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// traccc include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/finding/combinatorial_kalman_filter_algorithm.hpp"

// detray include(s).
#include "detray/core/detector.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/io/frontend/detector_reader.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

using namespace traccc;

int main()
{
    /// Type declarations
    using host_detector_type = detray::detector<detray::default_metadata,
                                                detray::host_container_types>;

    /*******************************
     * Read the telescope geometry
     *******************************/

    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;

    detray::io::detector_reader_config reader_cfg{};
    std::string file{__FILE__};
    std::string dir{file.substr(0, file.rfind("/"))};
    reader_cfg.add_file(dir + "/../geometry/telescope_detector_geometry.json");
    reader_cfg.add_file(dir + "/../geometry/telescope_detector_homogeneous_material.json");

    const auto [host_det, names] =
        detray::io::read_detector<host_detector_type>(host_mr, reader_cfg);

    /********************************************************************
     * Prepare measurements and track parameters estimated from seeding
     ********************************************************************/

    measurement_collection_types::host measurements;

    // Measurements from the 1st particle
    measurements.push_back({{0.0050959279760718346f, -0.03403015062212944f},
                            {0.0024999999441206455, 0.0024999999441206455},
                            detray::geometry::barcode{281474976710783},
                            0u});
    measurements.push_back({{3.1412324905395508f, -38.085170745849609f},
                            {0.0024999999441206455, 0.0024999999441206455},
                            detray::geometry::barcode{281475110928575},
                            1u});
    measurements.push_back({{12.875617980957031f, -76.4844970703125f},
                            {0.0024999999441206455, 0.0024999999441206455},
                            detray::geometry::barcode{281475245146367},
                            2u});
    measurements.push_back({{29.173341751098633f, -115.179931640625f},
                            {0.0024999999441206455, 0.0024999999441206455},
                            detray::geometry::barcode{281475379364159},
                            3u});
    measurements.push_back({{52.272651672363281f, -154.30535888671875f},
                            {0.0024999999441206455, 0.0024999999441206455},
                            detray::geometry::barcode{281475513581951},
                            4u});

    // Measurements from the 2nd particle
    measurements.push_back({{0.093167312443256378, 0.033621422946453094},
                            {0.0024999999441206455, 0.0024999999441206455},
                            detray::geometry::barcode{281474976710783},
                            5u});
    measurements.push_back({{3.5960044860839844, 67.298873901367188},
                            {0.0024999999441206455, 0.0024999999441206455},
                            detray::geometry::barcode{281475110928575},
                            6u});
    measurements.push_back({{14.534117698669434, 135.09417724609375},
                            {0.0024999999441206455, 0.0024999999441206455},
                            detray::geometry::barcode{281475245146367},
                            7u});
    measurements.push_back({{32.869495391845703, 203.50968933105469},
                            {0.0024999999441206455, 0.0024999999441206455},
                            detray::geometry::barcode{281475379364159},
                            8u});
    measurements.push_back({{59.011554718017578, 273.01675415039062},
                            {0.0024999999441206455, 0.0024999999441206455},
                            detray::geometry::barcode{281475513581951},
                            9u});

    // @@@@IMPORTANT@@@@
    // Measurements need to be sorted w.r.t. geometry barcode
    std::sort(measurements.begin(), measurements.end(),
              measurement_sort_comp());

    bound_track_parameters_collection_types::host params;

    // Initial estimation for the 1st particle's track parameters
    params.push_back(bound_track_parameters(
        detray::geometry::barcode{281474976710783},
        detray::bound_parameters_vector<traccc::default_algebra>({0.f, 0.f}, 0.f, 1.935, -1.f, 0.f),
        matrix::identity<detray::bound_matrix<traccc::default_algebra>>()));

    // Initial estimation for the 2nd particle's track parameters
    params.push_back(bound_track_parameters(
        detray::geometry::barcode{281474976710783},
        detray::bound_parameters_vector<traccc::default_algebra>({0.f, 0.f}, 0.f, 0.978, -1.f, 0.f),
        matrix::identity<detray::bound_matrix<traccc::default_algebra>>()));

    /******************************
     * Run Finding
     ******************************/

    using finding_algorithm =
        traccc::host::combinatorial_kalman_filter_algorithm;
    finding_algorithm::config_type finding_cfg;
    finding_cfg.propagation.stepping.rk_error_tol = 1e-8f * unit<float>::mm;
    finding_cfg.min_track_candidates_per_track = 5;

    finding_algorithm finding_alg(finding_cfg);

    const traccc::vector3 B{0, 0, 2 * detray::unit<traccc::scalar>::T};
    auto field = detray::bfield::create_const_field(B);

    auto track_candidates =
        finding_alg(host_det, field,
                    vecmem::get_data(measurements),
                    vecmem::get_data(params));

    std::cout << std::endl;
    std::cout << "Number of tracks found: " << track_candidates.size() << std::endl;
    std::cout << std::endl;

    for (std::size_t i = 0; i < track_candidates.size(); i++)
    {
        std::cout << i << "-th track's candidates:" << std::endl;
        auto &cands = track_candidates.at(i).items;

        for (std::size_t j = 0; j < cands.size(); j++)
        {
            std::cout << "meas id = " << cands.at(j).measurement_id << " | "
                      << cands.at(j).surface_link << std::endl;
        }
        std::cout << std::endl;
    }

    return 1;
}