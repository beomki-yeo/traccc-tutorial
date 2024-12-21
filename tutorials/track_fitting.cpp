/** TRACCC tutorial for beginners
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// traccc include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/fitting/kalman_filter/kalman_fitter.hpp"

// detray include(s).
#include "detray/core/detector.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/core/detector_metadata.hpp"
#include "detray/test/utils/detectors/build_telescope_detector.hpp"

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

    /*****************************
     * Build a geometry
     *****************************/

    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;

    // 2T magnetic field in z-axis
    const traccc::vector3 B{0, 0, 2 * detray::unit<traccc::scalar>::T};
    auto field = detray::bfield::create_const_field(B);

    // Senstive module dimension
    detray::mask<detray::rectangle2D> sensitive_rect{0u,
                                                     100.f * traccc::unit<scalar>::mm,
                                                     100.f * traccc::unit<scalar>::mm};
    // Module alignment directio in x-axis
    const vector3 align_axis{1.f, 0.f, 0.f};
    detray::detail::ray<traccc::default_algebra> pilot_track{
        {0, 0, 0}, 0, align_axis, -1};

    // Make a telescope geometry with 5 pixel modules
    detray::tel_det_config tel_cfg{sensitive_rect, pilot_track};
    tel_cfg.positions({100, 200, 300, 400, 500}); // unit in mm
    tel_cfg.module_material(detray::silicon<scalar>());
    tel_cfg.mat_thickness(80.f * unit<scalar>::um);
    tel_cfg.envelope(200);

    const auto [det, name_map] = build_telescope_detector(host_mr, tel_cfg);

    /******************************
     * Prepare input measurements
     ******************************/

    /******************************
     * Run Fitting
     ******************************/

    // Fitting algorithm object
    traccc::fitting_config fit_cfg;
    fit_cfg.use_backward_filter = true;
    fit_cfg.covariance_inflation_factor = 1e3f;


    return 1;
}