#ifndef SIM_DATA_HPP
#define SIM_DATA_HPP 1

#include "initial_data.hpp"
#include "units.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "fermi_table.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/source_terms/CenterGravity.hpp"
#include "monopole_self_gravity.hpp"
#include "source/newtonian/two_dimensional/source_terms/cylindrical_complementary.hpp"
#include "source/newtonian/two_dimensional/source_terms/SeveralSources.hpp"
#include "inner_bc.hpp"
#include "lazy_extensive_updater.hpp"
#include "lazy_cell_updater.hpp"
#include "create_grid.hpp"
#include "generate_atomic_properties.hpp"
#include "source/misc/vector_initialiser.hpp"
#include "circular_section.hpp"
#include "calc_init_cond.hpp"

class SimData
{
public:

  SimData(const InitialData& id,
	  const Units& u,
	  const CircularSection& domain);

  hdsim& getSim(void);

  const FermiTable& getEOS(void) const;

private:
  const CylindricalSymmetry pg_;
  const SquareBox outer_;
  VoronoiMesh tess_;
  const FermiTable eos_;
  const Hllc rs_;
  Eulerian point_motion_;
  CenterGravity gravity_acc_;
  ConservativeForce gravity_force_;
  MonopoleSelfGravity msg_;
  CylindricalComplementary geom_force_;
  SeveralSources force_;
  const SimpleCFL tsf_;
  const InnerBC fc_;
  const LazyExtensiveUpdater eu_;
  const LazyCellUpdater cu_;
  hdsim sim_;
};

#endif // SIM_DATA_HPP
