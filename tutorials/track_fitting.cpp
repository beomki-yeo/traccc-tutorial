/** TRACCC tutorial for beginners
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// traccc include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/fitting/kalman_filter/kalman_fitter.hpp"
#include "traccc/geometry/detector.hpp"

// detray include(s).
#include "detray/core/detector.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/core/detector_metadata.hpp"
#include "detray/test/utils/detectors/build_telescope_detector.hpp"
#include "detray/io/frontend/detector_reader.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

using namespace traccc;

int main()
{
    /// Type declarations
    using host_detector_type = detray::detector<detray::default_metadata,
                                                detray::host_container_types>;

    using b_field_t = covfie::field<detray::bfield::const_bknd_t>;
    using rk_stepper_type =
        detray::rk_stepper<b_field_t::view_t, traccc::default_algebra,
                           detray::constrained_step<>>;

    using host_navigator_type = detray::navigator<const host_detector_type>;
    using host_fitter_type =
        traccc::kalman_fitter<rk_stepper_type, host_navigator_type>;

    /*******************************
     * Read the telescope geometry
     *******************************/

    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;

    detray::io::detector_reader_config reader_cfg{};
    std::string file {__FILE__};
    std::string dir {file.substr(0, file.rfind("/"))};
    reader_cfg.add_file(dir + "/../geometry/telescope_detector_geometry.json");
    reader_cfg.add_file(dir + "/../geometry/telescope_detector_homogeneous_material.json");

    const auto [host_det, names] =
        detray::io::read_detector<host_detector_type>(host_mr, reader_cfg);

    /*
    traccc::io::read_detector(
        host_detector, host_mr, "telescope_detector_geometry.json",
        "telescope_detector_homogeneous_material.json");
    */
    /******************************
     * Prepare input measurements
     ******************************/

    /******************************
     * Run Fitting
     ******************************/
    
    // Fitting algorithm object
    traccc::fitting_config fit_cfg;
    fit_cfg.use_backward_filter = true;

    return 1;
}