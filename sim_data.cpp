#include "sim_data.hpp"

SimData::SimData(const InitialData& id, const Units& u):
  pg_(Vector2D(0,0), Vector2D(1,0)),
  outer_(Vector2D(-0.5*id.radius_mid.front(),0.9*id.radius_mid.front()),
	 Vector2D(0.5*id.radius_mid.front(),1.2*id.radius_mid.back())),
  tess_(create_grid(outer_.getBoundary(),1e-2,0.9*id.radius_list.front()),
	outer_),
  eos_("eos_tab.coded",1,1,0,generate_atomic_properties()),
  rs_(),
  point_motion_(),
  gravity_acc_(u.gravitation_constant*u.core_mass,
	       0, Vector2D(0,0)),
  gravity_force_(gravity_acc_),
  msg_(linspace(id.radius_list.front(),id.radius_list.back(),100),
       u.gravitation_constant,
       pair<double,double>(M_PI*0.46,M_PI*0.54)),
  geom_force_(pg_.getAxis()),
  force_(VectorInitialiser<SourceTerm*>(&gravity_force_)
	 (&geom_force_)(&msg_)()),
  tsf_(0.3),
  fc_(rs_,string("ghost"),id.radius_mid.back(),
      u.gravitation_constant*u.core_mass/
      pow(id.radius_list.back(),2)),
  eu_(),
  cu_(),
  sim_(tess_,
       outer_,
       pg_,
       calc_init_cond(tess_,eos_,id,CircularSection(id.radius_mid.front(),
						    id.radius_mid.back(),
						    0.46*M_PI,
						    0.54*M_PI)),
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
