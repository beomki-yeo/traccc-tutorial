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
#include "traccc/fitting/kalman_fitting_algorithm.hpp"

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
        detray::io::read_detector<traccc::default_detector::host>(host_mr, reader_cfg);

    /***************************
     * Prepare track candidate
     ***************************/

    track_candidate_container_types::host track_candidates;

    // There are two tracks to fit
    track_candidates.resize(2u);

    // Initial esitmation for the 1st particle's track parameters
    track_candidates.at(0u).header = bound_track_parameters(
        detray::geometry::barcode{281474976710783},
        detray::bound_parameters_vector<traccc::default_algebra>({10.f, -10.f}, 0.f, constant<scalar>::pi_2, -1.1f, 0.f),
        matrix::identity<detray::bound_matrix<traccc::default_algebra>>());

    // Initial esitmation for the 2nd particle's track parameters
    track_candidates.at(1u).header = bound_track_parameters(
        detray::geometry::barcode{281474976710783},
        detray::bound_parameters_vector<traccc::default_algebra>({10.f, -10.f}, 0.f, constant<scalar>::pi_2, -0.9f, 0.f),
        matrix::identity<detray::bound_matrix<traccc::default_algebra>>());

    // Measurements from the 1st particle
    track_candidates.at(0u).items.push_back({{0.0050959279760718346f, -0.03403015062212944f},
                                             {0.0024999999441206455, 0.0024999999441206455},
                                             detray::geometry::barcode{281474976710783}});
    track_candidates.at(0u).items.push_back({{3.1412324905395508f, -38.085170745849609f},
                                             {0.0024999999441206455, 0.0024999999441206455},
                                             detray::geometry::barcode{281475110928575}});
    track_candidates.at(0u).items.push_back({{12.875617980957031f, -76.4844970703125f},
                                             {0.0024999999441206455, 0.0024999999441206455},
                                             detray::geometry::barcode{281475245146367}});
    track_candidates.at(0u).items.push_back({{29.173341751098633f, -115.179931640625f},
                                             {0.0024999999441206455, 0.0024999999441206455},
                                             detray::geometry::barcode{281475379364159}});
    track_candidates.at(0u).items.push_back({{52.272651672363281f, -154.30535888671875f},
                                             {0.0024999999441206455, 0.0024999999441206455},
                                             detray::geometry::barcode{281475513581951}});

    // Measurements from the 2nd particle
    track_candidates.at(1u).items.push_back({{0.093167312443256378, 0.033621422946453094},
                                             {0.0024999999441206455, 0.0024999999441206455},
                                             detray::geometry::barcode{281474976710783}});
    track_candidates.at(1u).items.push_back({{3.5960044860839844, 67.298873901367188},
                                             {0.0024999999441206455, 0.0024999999441206455},
                                             detray::geometry::barcode{281475110928575}});
    track_candidates.at(1u).items.push_back({{14.534117698669434, 135.09417724609375},
                                             {0.0024999999441206455, 0.0024999999441206455},
                                             detray::geometry::barcode{281475245146367}});
    track_candidates.at(1u).items.push_back({{32.869495391845703, 203.50968933105469},
                                             {0.0024999999441206455, 0.0024999999441206455},
                                             detray::geometry::barcode{281475379364159}});
    track_candidates.at(1u).items.push_back({{59.011554718017578, 273.01675415039062},
                                             {0.0024999999441206455, 0.0024999999441206455},
                                             detray::geometry::barcode{281475513581951}});

    /******************************
     * Run Fitting
     ******************************/

    // Fitting algorithm object
    traccc::fitting_config fit_cfg;
    fit_cfg.propagation.stepping.rk_error_tol = 1e-8f * unit<float>::mm;
    //@TIP: Kalman fitter can be repeated for more precise result
    //fit_cfg.n_iterations = 2;
    fit_cfg.use_backward_filter = true;
    traccc::host::kalman_fitting_algorithm fitting(fit_cfg, host_mr);

    const traccc::vector3 B{0, 0, 2 * detray::unit<traccc::scalar>::T};
    auto field = detray::bfield::create_const_field(B);

    // Run fitting
    auto track_states =
        fitting(host_det, field, traccc::get_data(track_candidates));

    const scalar q = -1.f;

    std::cout << std::endl;
    std::cout << "---- 1st track fitting result ----" << std::endl;
    std::cout << "NDF: " << track_states.at(0u).header.ndf << std::endl;
    std::cout << "Chi2: " << track_states.at(0u).header.chi2 << std::endl;
    std::cout << "Fitted momentum [GeV/c]: " << track_states.at(0u).header.fit_params.p(q) << std::endl;
    std::cout << std::endl;

    std::cout << "---- 2nd track fitting result ----" << std::endl;
    std::cout << "NDF: " << track_states.at(1u).header.ndf << std::endl;
    std::cout << "Chi2: " << track_states.at(1u).header.chi2 << std::endl;
    std::cout << "Fitted momentum [GeV/c]: " << track_states.at(1u).header.fit_params.p(q) << std::endl;
    std::cout << std::endl;

    return 1;
}