/** TRACCC tutorial for beginners
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// traccc include(s).
#include "traccc/definitions/common.hpp"

// detray include(s).
#include "detray/core/detector.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/test/utils/detectors/build_telescope_detector.hpp"
#include "detray/io/frontend/detector_writer.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

using namespace traccc;

int main()
{

    /*****************************
     * Build and write a geometry
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
    tel_cfg.positions({0, 100, 200, 300, 400}); // unit in mm
    tel_cfg.module_material(detray::silicon<scalar>());
    tel_cfg.mat_thickness(80.f * unit<scalar>::um);

    const auto [det, name_map] = build_telescope_detector(host_mr, tel_cfg);

    // Create detector file
    auto writer_cfg = detray::io::detector_writer_config{}
                          .format(detray::io::format::json)
                          .replace_files(true);
    detray::io::write_detector(det, name_map, writer_cfg);

    return 1;
}