#include "sim_data.hpp"
#include "calc_bottom_area.hpp"

SimData::SimData(const InitialData& id,
		 const Units& u,
		 const CircularSection& domain):
  pg_(Vector2D(0,0), Vector2D(1,0)),
  outer_(Vector2D(-0.5*id.radius_mid.front(),0.9*id.radius_mid.front()),
	 Vector2D(0.5*id.radius_mid.front(),1.2*id.radius_mid.back())),
  tess_(create_grid(outer_.getBoundary(),1e-2,0.9*id.radius_list.front()),
	outer_),
  eos_("eos_tab.coded",1,1,0,generate_atomic_properties()),
  rs_(),
  point_motion_(),
  cag_
  (u.core_mass,
   linspace(id.radius_list.front(),id.radius_list.back(),100),
   u.gravitation_constant,
   domain.getAngles()),
  geom_force_(pg_.getAxis()),
  force_(VectorInitialiser<SourceTerm*>
	 (&cag_)
	 (&geom_force_)
	 ()),
  tsf_(0.3),
  fc_(),
  eu_(),
  cu_(),
  sim_(tess_,
       outer_,
       pg_,
       calc_init_cond(tess_,eos_,id,domain),
       eos_,
       point_motion_,
       force_,
       tsf_,
       fc_,
       eu_,
       cu_) {}

hdsim& SimData::getSim(void)
{
  return sim_;
}

const FermiTable& SimData::getEOS(void) const
{
  return eos_;
}
