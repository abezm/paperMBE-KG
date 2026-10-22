#include "stdafx.h"
#include "FEKinematicGrowthZ.h"
#include "FEUncoupledMaterial.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/DumpStream.h>
#include <FECore/log.h>
#define _USE_MATH_DEFINES				// M_PI
#include <cmath>
#include "FECore/log.h"				// to print to log file and/or screen
#include <limits>
#include <algorithm>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEKinematicGrowthZ, FEElasticMaterial)

	// material properties
	ADD_PROPERTY(m_Base, "elastic");
	
	// ##############  Hyperelastic model ##############
	ADD_PARAMETER(m_mu    , FE_RANGE_GREATER_OR_EQUAL(0.0), "c_e"); 
	ADD_PARAMETER(m_cm    , FE_RANGE_GREATER_OR_EQUAL(0.0), "c_m");
	ADD_PARAMETER(m_dm    , FE_RANGE_GREATER(0.0), "d_m");
	ADD_PARAMETER(m_cc    , FE_RANGE_GREATER_OR_EQUAL(0.0), "c_c");
	ADD_PARAMETER(m_dc    , FE_RANGE_GREATER(0.0), "d_c");
	ADD_PARAMETER(m_alpha , "alpha");
	
	// ############## mass fractions ###################
	
	// constituent mass fractions (elastin, collagen)
	ADD_PARAMETER(m_phieo, FE_RANGE_CLOSED(0.0,1.0), "f_e"); 
	ADD_PARAMETER(m_phico, FE_RANGE_CLOSED(0.0,1.0), "f_c"); 
	
	// collagen relative fractions 
	ADD_PARAMETER(m_betaz, FE_RANGE_CLOSED(0.0,1.0), "beta_z"); 
	ADD_PARAMETER(m_betat, FE_RANGE_CLOSED(0.0,1.0), "beta_t"); 

	// ##############  Kinematic Growth model ##############
	ADD_PARAMETER(m_kPar, FE_RANGE_GREATER(0.0), "k_parameter");			// rate parameter K_sigma k (rate-dependent growth)
	ADD_PARAMETER(m_sgrSta, FE_RANGE_GREATER_OR_EQUAL(0.0), "gr_start");		// end of Stage I, start of G&R
	ADD_PARAMETER(m_sgrEnd, FE_RANGE_GREATER_OR_EQUAL(0.0), "gr_end");		// end of the degradation ramp
	
	ADD_PARAMETER(m_MBE_opt,  "MBE");						// 1: equilibrated growth (MBE), 0: rate-dependent growth
	ADD_PARAMETER(m_asym,     "asymm_aneurysm");					// degradation profile, 0: axisymmetric, 1: asymmetric
	ADD_PARAMETER(m_lo_opt,   "lo_opt");						// 1: F_e = F_o after gr_start
	ADD_PARAMETER(m_writeRes, "pr_opt");						// 1: debug output to the log
	
	ADD_PARAMETER(m_kJ, "k_J"); // volumetric penalty of Stage I

	// ############## deposition stretches ##############
	ADD_PARAMETER(m_Get, FE_RANGE_GREATER_OR_EQUAL(0.0), "Ge_t");
	ADD_PARAMETER(m_Gez, FE_RANGE_GREATER_OR_EQUAL(0.0), "Ge_z");
	ADD_PARAMETER(m_Gm,  FE_RANGE_GREATER_OR_EQUAL(0.0), "Gm");
	ADD_PARAMETER(m_Gc,  FE_RANGE_GREATER_OR_EQUAL(0.0), "Gc");
	
	ADD_PARAMETER(m_muout, FE_RANGE_CLOSED(0.0,1.0), "mu_out"); 
	ADD_PARAMETER(m_muin,  FE_RANGE_CLOSED(0.0,1.0), "mu_in"); 
	ADD_PARAMETER(m_ccout, FE_RANGE_CLOSED(0.0,1.0), "cc_out"); 
	ADD_PARAMETER(m_ccin,  FE_RANGE_CLOSED(0.0,1.0), "cc_in");

	ADD_PARAMETER(m_hypSel, "hyperels_mod"); // 1: neo-Hooke, 2: St. Venant-Kirchhoff, 3: elastin + SMC + four collagen families
	ADD_PARAMETER(m_degr,   "ec_degrad"); // 1: elastin and collagen degradation
	ADD_PARAMETER(m_alpm,   "alp_update"); // 1: update of the diagonal fibre angle
	
	// ############## secant tangent option ##############
	ADD_PARAMETER(m_secant_tangent, "secant_tangent");
	
	// ############## degradation profile ##############
	ADD_PARAMETER(m_z_center, "z_center");
	ADD_PARAMETER(m_z_spread, "z_spread");
	ADD_PARAMETER(m_nuz, "z_nu");
	ADD_PARAMETER(m_nut, "t_nu");
	ADD_PARAMETER(m_pif, "t_pif");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEKinematicMaterialPointZ::FEKinematicMaterialPointZ(FEMaterialPointData* mp) : FEMaterialPointData(mp)
{
	m_sed = 0.0;
	m_sedp = 0.0;
	
	m_Jp = 1.0;
	m_J = 1.0;
	m_Fg.unit();
	m_Fgp.unit();
	
	m_svo = 0.0;
	
	m_Fo.unit();
	
	m_Wt = 0.0;
	
}

//-----------------------------------------------------------------------------
//! Create a shallow copy of the material point data
FEMaterialPointData* FEKinematicMaterialPointZ::Copy()
{
	FEKinematicMaterialPointZ* pt = new FEKinematicMaterialPointZ(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
//! Initializes material point data.
void FEKinematicMaterialPointZ::Init()
{
	// intialize data to zero
	m_Se.zero();
	m_Sep.zero();
	m_sed = 0.0;
    m_sedp = 0.0;
    m_Jp = 1.0;
    m_J = 1.0;
    m_Fg.unit();
    m_Fgp.unit();
    
    m_Fo.unit();
    
    m_svo = 0.0;	
    
    m_elas_dmg = 1.0;
    m_coll_dmg = 1.0;	
    
    m_Wt = 0.0;

    // don't forget to initialize the base class
    FEMaterialPointData::Init();
}

//-----------------------------------------------------------------------------
//! Update material point data.
void FEKinematicMaterialPointZ::Update(const FETimeInfo& timeInfo)
{

	// the elastic stress stored in pt is the Cauchy stress.
	// however, we need to store the 2nd PK stress
	m_Sep = m_Se;
    m_sedp = m_sed;
    m_Jp = m_J;	// updates J_n at end of current timestep
    m_Fgp = m_Fg;	// updates Fg_n at end of current timestep
    
    m_Wt = m_sed/m_J;
    
    // don't forget to call the base class
    FEMaterialPointData::Update(timeInfo);
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FEKinematicMaterialPointZ::Serialize(DumpStream& ar)
{
    FEMaterialPointData::Serialize(ar);
	ar & m_Se;
	ar & m_Sep;
    ar & m_sed & m_sedp;
    ar & m_J & m_Jp;
    ar & m_Fg & m_Fgp;
    ar & m_svo;
    ar & m_Fo;	
    
    ar & m_elas_dmg & m_coll_dmg;
    
    ar & m_Wt;
}

//-----------------------------------------------------------------------------
//! constructor
FEKinematicGrowthZ::FEKinematicGrowthZ(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_Base = 0;
}

//-----------------------------------------------------------------------------
//! get the elastic base material \todo I want to call this GetElasticMaterial, but this name is being used
FEElasticMaterial* FEKinematicGrowthZ::GetBaseMaterial()
{ 
	return m_Base; 
}

//-----------------------------------------------------------------------------
//! Set the base material
void FEKinematicGrowthZ::SetBaseMaterial(FEElasticMaterial* pbase)
{ 
	m_Base = pbase; 
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPointData* FEKinematicGrowthZ::CreateMaterialPointData()
{
	return new FEKinematicMaterialPointZ(m_Base->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! Stress function
mat3ds FEKinematicGrowthZ::Stress(FEMaterialPoint& mp)
{
	// get machine epsilon 
	double eps = std::numeric_limits<double>::epsilon();
	
	// Get current time, end time and time-step
    	double t = GetFEModel()->GetTime().currentTime;

	double dt = GetTimeInfo().timeIncrement;

	double partialtime = m_sgrEnd;			// end of the degradation ramp
	double sgr = std::min(t,partialtime);		// time of the degradation ramp

	if (dt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
    
	// get the elastic part
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the KinematicGrowth point data
	FEKinematicMaterialPointZ& pt = *mp.ExtractData<FEKinematicMaterialPointZ>();
	
	// Get Jacobian of deformation gradient, J_n (previous time-step)
	double Jp = pt.m_Jp;	// J_n
	
	// Get total deformation gradient, F_n+1
	mat3d F = ep.m_F;		// F_n+1
	double J = ep.m_J;		// J_n+1 
	pt.m_J = J;			// saves J_n+1 to later update Jp at the end of current timestep
	
	// Get growth deformation gradient, Fg_n+1
	mat3dd I(1.0);				// identity (2nd order)
	mat3d  Fg = pow(J,1.0/3.0)*I; if (t <= m_sgrSta + eps) Fg = mat3d(1,0,0, 0,1,0, 0,0,1);	// isotropic-volumetric growth ONLY; F_g = 1 in Stage I
	mat3d  Fgi = Fg.inverse();		// Fg_n+1 inverse
	double Jgi = Fgi.det();		// det Fg^-1 (not used)
	pt.m_Fg = Fg;				// saves Fg_n+1 to later update Jp at the end of current timestep
	
	// elastic deformation gradient, Fe_n+1 (for elastic stress calculation)
	mat3d Fe;
	
	mat3d Fo = pt.m_Fo;
	
	Fe = F*Fgi;
	if (m_lo_opt && t > m_sgrSta + eps){
		Fe = Fo;	// lo_opt = 1: F_e = F_o after gr_start
	}

	// MATERIAL PARAMETERS ----------------------------

	// mass/volume fractions
	double phieo = m_phieo;
	double phico = m_phico;
	double phimo = 1.0 - phieo - phico;
	// f_e + f_c > 1: the three mass fractions are set to 1
	if (phieo + phico > 1.0){
		phieo = 1.0;
		phico = 1.0;
		phimo = 1.0;
	}

	double mu  = m_mu; 
	
	double cm = m_cm;
	double dm = m_dm;
	double cc = m_cc;
	double dc = m_dc;

	double Gm = m_Gm;
	double Gc = m_Gc;
	double Get = m_Get;
	double Gez = m_Gez;
	
	double betat = m_betat;
	double betaz = m_betaz;
	
	double betad = 0.5*(1.0 - betat - betaz);
	// beta_t + beta_z > 1: beta_t = beta_z = 1, beta_d = 0.5
	if (betat + betaz > 1.0){
		betat = 1.0; 
		betaz = 1.0; 
		betad = 0.5; 
	}
	
	double alpha_d = m_alpha;

	double aexp = 1.0;									// fibre-angle update: tan(alpha) = tan(alpha_0) (lt/lz)^aexp
	double alpha = alpha_d*M_PI/180;	// rad
	
	// retrieve material position
	vec3d X = mp.m_r0;
	vec3d Xcl = {0.0, 0.0, X.z};				// center line
	vec3d NX = {X.x-Xcl.x,X.y-Xcl.y,X.z-Xcl.z};		// radial vector
	double ro = sqrt(NX*NX);
	NX /= ro;

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// local basis: Nz axial (local system of the mesh), Nt circumferential, Nr radial
	vec3d Nz(Q(0,2),Q(1,2),Q(2,2)); 	// axial
	vec3d Nt = vec3d(-NX.y, NX.x, 0.0);  // circ
	vec3d Nr = Nz ^ Nt; 			// radial
	double lenz = Nz.norm(); double lent = Nt.norm(); double lenr = Nr.norm();
	if (lenz > 1e-12) Nz /= lenz;
	if (lent > 1e-12) Nt /= lent;
	if (lenr > 1e-12) Nr /= lenr;
	vec3d N[3] = {Nr, Nt, Nz};
	
    	double lt = 1.;
    	double lz = 1.;
    	if (t > m_sgrSta + eps && m_alpm==1) {
    	
    		// update alpha
		// stretches measured from the original homeostatic state (F Fo^-1), not from the mesh (F)
		mat3d Fr = F * Fo.inverse();
		lt = (Fr * N[1]).norm();
		lz = (Fr * N[2]).norm();
		alpha = atan(tan(alpha) * pow(lt / lz, aexp));
    	
	    	}
    	
    	vec3d Np = Nt * sin(alpha) + Nz * cos(alpha);
    	vec3d Nn = Nt * sin(alpha) - Nz * cos(alpha);

	double azimuth = acos(-NX.y);			// azimuth wrt axis -Y
	
	mat3ds Ge = I;	// Ge (unit matrix)
	// Ge from spectral decomposition
	Ge = 1.0/Get/Gez*dyad(N[0]) + Get*dyad(N[1]) + Gez*dyad(N[2]);

	double muout = mu;
	double muin  = mu;
	
	double ccout  = cc;
	double ccin   = cc;

	if (t > m_sgrSta + eps && m_degr==1) {
		//	ELASTIN DAMAGE
		muout = m_muout*mu;
		muin  = m_muin*muout;
		
		if (m_asym == 0) {
			mu = muout+(muin-muout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz));	// axisymmetric aneurysm
		}
		if (m_asym == 1) {
			mu = muout+(muin-muout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz))*exp(-pow(std::abs((azimuth-M_PI)/(M_PI/m_pif)),m_nut)); 	// asymmetric aneurysm	
		}
	
		//	REDUCED COLLAGEN CROSS-LINKING
		ccout  = m_ccout*cc;
		ccin   = m_ccin*ccout;
		
		if (m_asym == 0) {
			cc = ccout+(ccin-ccout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz));	// axisymmetric aneurysm
		}
		if (m_asym == 1) {
			cc = ccout+(ccin-ccout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz))*exp(-pow(std::abs((azimuth-M_PI)/(M_PI/m_pif)),m_nut)); 	// asymmetric aneurysm	
		}
	}
	
	pt.m_elas_dmg = mu/m_mu; 
	pt.m_coll_dmg = cc/m_cc;
	
	mat3ds Se(0.0);
	
	double lto = 1.;
	double lzo = 1.;
	double lpo = 1.;
	double lno = 1.;
	
	double lmt2 = 1.; 
	double lct2 = 1.;
	double lcz2 = 1.; 
	double lcp2 = 1.; 
	double lcn2 = 1.;
	
	if (m_hypSel == 1) Se = mu*I;	// Neo-hooke
	if (m_hypSel == 2){		// St.Venant-Kirchhoff
		mat3ds Ae = ((Fe.transpose()*Fe-I).sym())*0.5; // Green-Lagrange elastic strain
		Se = Ae*2.0*mu;	
	}
	if (m_hypSel == 3){		// 4-FIBER FAMILY MOD. 
		
		lto = (Fe*N[1]).norm();
		lzo = (Fe*N[2]).norm();
		lpo = (Fe*Np).norm();
		lno = (Fe*Nn).norm();

		lmt2 = (Gm*lto)*(Gm*lto);
		lct2 = (Gc*lto)*(Gc*lto);
		lcz2 = (Gc*lzo)*(Gc*lzo);
		lcp2 = (Gc*lpo)*(Gc*lpo);						// deposition stretch calculation (computational purposes)
		lcn2 = (Gc*lno)*(Gc*lno);						// idem for symmetric
		
		mat3ds Sel = mu*Ge.sqr();	
		
		double Shm  = cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm);
		
		double Shct = cc*(lct2-1.0)*exp(dc*(lct2-1.0)*(lct2-1.0))*(Gc*Gc);
		double Shcz = cc*(lcz2-1.0)*exp(dc*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc);
		double Shcp = cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc);
		double Shcn = cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc);
		
		mat3ds Sm = Shm*dyad(N[1]);
		mat3ds Sc = (  Shct*dyad(N[1])*betat + Shcz*dyad(N[2])*betaz +
				Shcp*dyad( Np )*betad + Shcn*dyad( Nn )*betad );
				
		Se = phieo * Sel + phimo * Sm + phico * Sc;	// elastin + SMCs + collagen
		
	}

	// extra PK-2 stress, Sx_n+1
	// isotropic growth ONLY
	mat3ds Sx = pow(J,1.0/3.0) * Se;
	
	// Cauchy-Green strain
	mat3ds C = ep.RightCauchyGreen();
	mat3ds Ci = C.inverse();
	
	// volumetric extra Cauchy stress, svx_n+1
	double svx = 1.0/3.0/J*Sx.dotdot(C);
	
	mat3ds S(0.0);
	double svo = 0.0;
	
	double Jdep = 0.9999;	// the penalty term k_J ln(Jdep J) is zero at J = 1/Jdep
	double lm = m_kJ;
	
	if (t <= m_sgrSta + eps) {
    	
		// total Cauchy stress, s_n+1
	    	S = J*Se + Ci*lm*log(Jdep*J);	// Stage I: S^x = J S_e, since F_g = 1
	    	
	    	// store volumetric Cauchy stress  (original-homeostatic)
	    	svo   = 1.0/3.0/J*S.dotdot(C);
	    	pt.m_svo = svo;
	    	
	    	// store Fo
	    	pt.m_Fo = F;
    		
    	}
    	else {
	    	// retrieve volumetric Cauchy stress  (original-homeostatic)
	    	svo = pt.m_svo;
	    	
	    	// Lagrange multiplier, p_n+1
	    	
	    	double p = 0.0;
	    	if (m_MBE_opt == 0){ // (RATE DEPENDENT)
	    		p = svx - svo*(1.0+log(J/Jp)/(m_kPar*dt));	
		}
		if (m_MBE_opt == 1) p = svx - svo;	// (MBE)
		
	    	// total PK-2 stress, S_n+1
	    	S = Sx - J*p*Ci;   
	    	
	    	// debug output (pr_opt = 1)
		if (m_writeRes) {
		    	char *msg = "F = %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"
				"Fo = %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"
    				"S = %f\t%f\t%f\t%f\t%f\t%f\n"
    				"--------------------------------------------\n";
			feLog(msg, 
    			F(0,0),F(1,1),F(2,2),F(0,1),F(1,0),F(0,2),F(2,0),F(1,2),F(2,1),
			Fo(0,0),Fo(1,1),Fo(2,2),Fo(0,1),Fo(1,0),Fo(0,2),Fo(2,0),Fo(1,2),Fo(2,1), 
    			S.xx(),S.yy(),S.zz(),S.xy(),S.xz(),S.yz());
    		}
	    	
	}
	
	// total cauchy stress (push-forward of S)
	mat3ds sigma = 1.0/J*(F*(S*F.transpose())).sym();
	
	// debug output (pr_opt = 1)
	if (m_writeRes) {
	    	char *msg = "F = %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"
    			"sig = %f\t%f\t%f\t%f\t%f\t%f\n"
    			"lo = %f\t%f\t%f\t%f\n"
    			"Np = %f\t%f\t%f\n"
    			"Nn = %f\t%f\t%f\n"
    			"Nz = %f\t%f\t%f\n"
    			"Nc = %f\t%f\t%f\n"
    			"--------------------------------------------\n";
		feLog(msg,
    		F(0,0),F(1,1),F(2,2),F(0,1),F(1,0),F(0,2),F(2,0),F(1,2),F(2,1),
    		sigma.xx(),sigma.yy(),sigma.zz(),sigma.xy(),sigma.xz(),sigma.yz(),
    		lto, lzo, lpo, lno,
    		Np(0), Np(1), Np(2),
    		Nn(0), Nn(1), Nn(2),
    		N[2].x, N[2].y, N[2].z,
    		N[1].x, N[1].y, N[1].z
    		);
    	}

	return sigma;
}

//-----------------------------------------------------------------------------
//! Material tangent
tens4dmm FEKinematicGrowthZ::TangentMm(FEMaterialPoint& mp)	
{

	// get machine epsilon 
	double eps = std::numeric_limits<double>::epsilon();
	
	// Get current time, end time and time-step
    	double t = GetFEModel()->GetTime().currentTime;

	double dt = GetTimeInfo().timeIncrement;

	double partialtime = m_sgrEnd;			// end of the degradation ramp
	double sgr = std::min(t,partialtime);		// time of the degradation ramp
 
 	// Get elastic material point data (def grad, F_n+1)
 	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
 	
 	// get the KinematicGrowth point data 
	FEKinematicMaterialPointZ& pt = *mp.ExtractData<FEKinematicMaterialPointZ>();
 	
 	// Get Jacobian of deformation gradient, J_n (previous time-step)
 	double Jp = pt.m_Jp;	// J_n
 	
 	// Get total deformation gradient, F_n+1
	mat3d &F = ep.m_F;		// F_n+1
	double J = ep.m_J;		// J_n+1 
	
	// Dyadic products of identity
	mat3dd I(1.0);				// identity (2nd order)
	tens4dmm IoI = dyad4mm(I,I);
	tens4dmm IxI = dyad1mm(I,I);
	
	// Get growth deformation gradient, Fg_n+1
	mat3d  Fg = pow(J,1.0/3.0)*I; if (t <= m_sgrSta + eps) Fg = mat3d(1,0,0, 0,1,0, 0,0,1);	// isotropic-volumetric growth ONLY; F_g = 1 in Stage I
	mat3d  Fgi = Fg.inverse();		// Fg_n+1 inverse
	double Jgi = Fgi.det();		// det Fg^-1 (not used)
	
	// elastic deformation gradient, Fe_n+1 (for elastic stress calculation)
	mat3d Fe;
	
	mat3d Fo = pt.m_Fo;
	
	Fe = F*Fgi;
	if (m_lo_opt && t > m_sgrSta + eps){
		Fe = Fo;	// lo_opt = 1: F_e = F_o after gr_start
	}

	// MATERIAL PARAMETERS ----------------------------

	// mass/volume fractions
	double phieo = m_phieo;
	double phico = m_phico;
	double phimo = 1.0 - phieo - phico;
	// f_e + f_c > 1: the three mass fractions are set to 1
	if (phieo + phico > 1.0){
		phieo = 1.0;
		phico = 1.0;
		phimo = 1.0;
	}

	double mu  = m_mu; 
	
	double cm = m_cm;
	double dm = m_dm;
	double cc = m_cc;
	double dc = m_dc;

	double Gm = m_Gm;
	double Gc = m_Gc;
	double Get = m_Get;
	double Gez = m_Gez;
	
	double betat = m_betat;
	double betaz = m_betaz;
	
	double betad = 0.5*(1.0 - betat - betaz);
	// beta_t + beta_z > 1: beta_t = beta_z = 1, beta_d = 0.5
	if (betat + betaz > 1.0){
		betat = 1.0; 
		betaz = 1.0; 
		betad = 0.5; 
	}
	
	double alpha_d = m_alpha;

	double aexp = 1.0;			// fibre-angle update: tan(alpha) = tan(alpha_0) (lt/lz)^aexp
	double alpha = alpha_d*M_PI/180;	// rad

	// retrieve material position
	vec3d X = mp.m_r0;
	vec3d Xcl = {0.0, 0.0, X.z};				// center line
	vec3d NX = {X.x-Xcl.x,X.y-Xcl.y,X.z-Xcl.z};		// radial vector
	double ro = sqrt(NX*NX);
	NX /= ro;

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);
	
	// local basis: Nz axial (local system of the mesh), Nt circumferential, Nr radial
	vec3d Nz(Q(0,2),Q(1,2),Q(2,2)); 	// axial
	vec3d Nt = vec3d(-NX.y, NX.x, 0.0);  // circ
	vec3d Nr = Nz ^ Nt; 			// radial
	double lenz = Nz.norm(); double lent = Nt.norm(); double lenr = Nr.norm();
	if (lenz > 1e-12) Nz /= lenz;
	if (lent > 1e-12) Nt /= lent;
	if (lenr > 1e-12) Nr /= lenr;
	vec3d N[3] = {Nr, Nt, Nz};

	mat3ds tent = dyad(Nt);
	mat3ds tenz = dyad(Nz);

    	mat3ds dtanadC(0.0);
    	vec3d dNpdta(0.0);
    	vec3d dNndta(0.0);
    	
    	if (t > m_sgrSta + eps && m_alpm==1) {
    		
    		// update alpha
		// stretches measured from the original homeostatic state (F Fo^-1), not from the mesh (F)
		mat3d Foi = Fo.inverse();
		vec3d Mt = Foi * Nt;
		vec3d Mz = Foi * Nz;
		double lt = (F * Mt).norm();
		double lz = (F * Mz).norm();
		alpha = atan(tan(alpha) * pow(lt / lz, aexp));

		// derivative of the fiber direction vectors Np and Nn with respect to tan(alpha)
		dNpdta = (Nt - Nz*tan(alpha))*pow(1 + pow(tan(alpha),2),-1.5);	// d(Np)/d(tan(alpha))
		dNndta = (Nt + Nz*tan(alpha))*pow(1 + pow(tan(alpha),2),-1.5);
		
		// derivative 2 * d(tan(alpha)) / dC
		// lt^2 = Mt.C.Mt and lz^2 = Mz.C.Mz, so the derivative uses dyad(Mt), dyad(Mz) in place of dyad(Nt), dyad(Nz)
		dtanadC = (dyad(Mt) * (1.0 / (lt*lt)) - dyad(Mz) * (1.0 / (lz*lz))) * (tan(alpha)*aexp);	// 2*d(tan(alpha))/d(C) : (Ft)o(Ft)
    	}
    	
    	vec3d Np = Nt * sin(alpha) + Nz * cos(alpha);
    	vec3d Nn = Nt * sin(alpha) - Nz * cos(alpha);
    	
    	mat3ds tenp = dyad(Np);
	mat3ds tenn = dyad(Nn);

	double azimuth = acos(-NX.y);			// azimuth wrt axis -Y

	mat3ds Ge = I;	// Ge (unit matrix)
	// Ge from spectral decomposition
	N[0] = N[2]^N[1];
	Ge = 1.0/Get/Gez*dyad(N[0]) + Get*tent + Gez*tenz;

	double muout = mu;
	double muin  = mu;
	
	double ccout  = cc;
	double ccin   = cc;
	
	if (t > m_sgrSta + eps && m_degr==1) {
		//	ELASTIN DAMAGE
		muout = m_muout*mu;
		muin  = m_muin*muout;
		
		if (m_asym == 0) {
			mu = muout+(muin-muout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz));	// axisymmetric aneurysm
		}
		if (m_asym == 1) {
			mu = muout+(muin-muout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz))*exp(-pow(std::abs((azimuth-M_PI)/(M_PI/m_pif)),m_nut)); 	// asymmetric aneurysm	
		}

		//	REDUCED COLLAGEN CROSS-LINKING
		ccout  = m_ccout*cc;
		ccin   = m_ccin*ccout;
		
		if (m_asym == 0) {
			cc = ccout+(ccin-ccout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz));	// axisymmetric aneurysm
		}
		if (m_asym == 1) {
			cc = ccout+(ccin-ccout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz))*exp(-pow(std::abs((azimuth-M_PI)/(M_PI/m_pif)),m_nut)); 	// asymmetric aneurysm	
		}
	}
	
	pt.m_elas_dmg = mu/m_mu; 
	pt.m_coll_dmg = cc/m_cc;

	mat3ds Se(0.0);
	mat3ds dSedJ(0.0);
	tens4dmm Ce(0.0); 
	tens4ds Cm(0.0); 
	tens4ds Cc(0.0); 
	
	mat3ds dSpdta(0.0);
	mat3ds dSndta(0.0);
	tens4dmm Calph(0.0);
	tens4dmm calph(0.0);
	
	if (m_hypSel == 1) Se = mu*I;	// Neo-hooke
	if (m_hypSel == 2){		// St.Venant-Kirchhoff
	
		mat3ds CGe = (Fe.transpose()*Fe).sym();
		mat3ds Ae = (CGe-I)*0.5; 			   // Green-Lagrange elastic strain
		Se = Ae*2.0*mu;	
		
		Ce = dyad4s(I,I)*mu; // St.Venant-Kirchhoff
		
		dSedJ = (-2.0*mu/3.0)*pow(J,-5.0/3.0)*(F.transpose()*F).sym();

	}
	if (m_hypSel == 3){		// 4-FIBER FAMILY MOD. 
		
		double lto = (Fe*N[1]).norm();
		double lzo = (Fe*N[2]).norm();
		double lpo = (Fe*Np).norm();
		double lno = (Fe*Nn).norm();

		double lmt2 = (Gm*lto)*(Gm*lto);
		double lct2 = (Gc*lto)*(Gc*lto);
		double lcz2 = (Gc*lzo)*(Gc*lzo);
		double lcp2 = (Gc*lpo)*(Gc*lpo);						// deposition stretch calculation (computational purposes)
		double lcn2 = (Gc*lno)*(Gc*lno);						// idem for symmetric

		mat3ds Sel = mu*Ge.sqr();	
		
		double Shm  = cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm);
		
		double Shct = cc*(lct2-1.0)*exp(dc*(lct2-1.0)*(lct2-1.0))*(Gc*Gc);
		double Shcz = cc*(lcz2-1.0)*exp(dc*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc);
		double Shcp = cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc);
		double Shcn = cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc);
		
		mat3ds Sm = Shm*dyad(N[1]);
		mat3ds Sc = (  Shct*dyad(N[1])*betat + Shcz*dyad(N[2])*betaz +
				Shcp*dyad( Np )*betad + Shcn*dyad( Nn )*betad );
				
		Se = phieo * Sel + phimo * Sm + phico * Sc;	// elastin + SMCs + collagen
		
		// (dS.hat_a / dI4_a)
		double dShm  = cm*(1.0+2.0*dm*(lmt2-1.0)*(lmt2-1.0))*exp(dm*(lmt2-1.0)*(lmt2-1.0))*pow(Gm,4);
		
		double dShct = cc*(1.0+2.0*dc*(lct2-1.0)*(lct2-1.0))*exp(dc*(lct2-1.0)*(lct2-1.0))*pow(Gc,4);
		double dShcz = cc*(1.0+2.0*dc*(lcz2-1.0)*(lcz2-1.0))*exp(dc*(lcz2-1.0)*(lcz2-1.0))*pow(Gc,4);
		double dShcp = cc*(1.0+2.0*dc*(lcp2-1.0)*(lcp2-1.0))*exp(dc*(lcp2-1.0)*(lcp2-1.0))*pow(Gc,4);
		double dShcn = cc*(1.0+2.0*dc*(lcn2-1.0)*(lcn2-1.0))*exp(dc*(lcn2-1.0)*(lcn2-1.0))*pow(Gc,4);

		Cm = 2.0 * dShm*dyad1s(tent);
		Cc = 2.0 * (dShct*dyad1s(tent)*betat + dShcz*dyad1s(tenz)*betaz +
			  dShcp*dyad1s(tenp)*betad + dShcn*dyad1s(tenn)*betad);
		
		// Compute the variation of the 2nd Piola-Kirchhoff stress wrt tan(alpha)
		dSpdta = 2.0 * Shcp * pow(J,1./3.) * dyad1s(Np, dNpdta);
		dSndta = 2.0 * Shcn * pow(J,1./3.) * dyad1s(Nn, dNndta);
		// tangent: variation of the stress magnitude with tan(alpha), lcp2 = Gc^2 (Fe Np).(Fe Np)
		dSpdta += 2.0 * dShcp * ((Fe*Np) * (Fe*dNpdta)) * pow(J,1./3.) * dyad(Np);
		dSndta += 2.0 * dShcn * ((Fe*Nn) * (Fe*dNndta)) * pow(J,1./3.) * dyad(Nn);
		
		// Tangent stiffness (var. of alpha)
		tens4dmm Calphp = dyad1mm(dSpdta, dtanadC);
		tens4dmm Calphn = dyad1mm(dSndta, dtanadC);
		Calph = phico * betad * (Calphp + Calphn);
		
		Ce = phimo * Cm + phico * Cc;	// SMCs + collagen
		
		dSedJ = (-2.0/(3.0*J))*(phimo * dShm*lmt2*dyad(N[1]) + 
			 phico * dShct*lct2*dyad(N[1])*betat + phico * dShcz*lcz2*dyad(N[2])*betaz + 
			 phico * dShcp*lcp2*dyad( Np )*betad + phico * dShcn*lcn2*dyad( Nn )*betad);
		
		if(m_degr==1 && t > m_sgrSta + eps) {	// elastic terms dropped with ec_degrad = 1; exact only with lo_opt = 1 (F_e = F_o)
			Ce = Ce*0.0;
			dSedJ = dSedJ*0.0;	// zero with constant F_o
		}  
		
	}

	// spatial elastic tangent
	tens4dmm ce = Ce.pp(F)/J;
	
	calph = Calph.pp(F)/J;
	
	if (dt == 0.0) return ce;

	// EXTRA PK-2 stress, Sx_n+1
	mat3ds Sx = pow(J,1.0/3.0)*Se; // isotropic growth ONLY
	
	// Cauchy-Green strain
	mat3ds C = ep.RightCauchyGreen();
	mat3ds Ci = C.inverse();

	tens4dmm CixCi = dyad1mm(Ci,Ci);
	tens4dmm CioCi = dyad4mm(Ci,Ci);	// (Ci dyad4 Ci)_ijkl = 1/2 (Ci_ik Ci_jl + Ci_il Ci_jk) (Holzapfel, p. 255)
	
	// volumetric extra Cauchy stress, svx_n+1
	double svx = 1.0/3.0/J*Sx.dotdot(C);
	
	tens4dmm c(0.0);
	tens4dmm Cp(0.0);
	tens4dmm cp(0.0);
	mat3ds S(0.0);
	
	// retrieves volumetric Cauchy stress  (original-homeostatic)
	double svo = pt.m_svo;
	if (t <= m_sgrSta + eps) {
    
	    	double Jdep = 0.9999;	// the penalty term k_J ln(Jdep J) is zero at J = 1/Jdep
		double lm = m_kJ;
	    	
	    	// material volumetric tangent
	    	Cp = lm*(CixCi-2.0*log(Jdep*J)*CioCi);	// p = dW_vol/dC = d/dC (ln J)**2
	    	
	    	// spatial volumetric tangent
	    	cp = Cp.pp(F)/J;
	    	
	    	// total spatial tangent, c_n+1
	    	// Stage I: material tangent J Ce + J Se x C^-1 + Cp; spatial: J ce + (F Se F^T) x 1 + cp
	    	mat3ds FSeFt = (F*Se*F.transpose()).sym();
	    	mat3ds Iden(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
	    	c = J*ce + dyad1mm(FSeFt, Iden) + cp + calph;
	    	
	    	// total PK-2 stress, S_n+1 / only to print
	    	S = J*Se + Ci*lm*log(Jdep*J);	// Stage I: S^x = J S_e, since F_g = 1	
    	}
    	
    	else {
	    	
	    	// update Lagrange multiplier
	    	
	    	double p = 0.0;
	    	tens4dmm cpAl(0.0);		// Spatial Algorithmic tangent
	    	if (m_MBE_opt == 0){ // (RATE DEPENDENT)
	    	
	    		p = svx - svo*(1.0+log(J/Jp)/(m_kPar*dt));
	    		
	    		// Algorithmic volumetric tangent, (dp_n+1_alg)/(dC)
	    		tens4dmm CpAl = (J*svo/(m_kPar*dt))*((1.0+log(J/Jp))*CixCi-2.0*log(J/Jp)*CioCi);
	    		cpAl = CpAl.pp(F)/J;
		}
		if (m_MBE_opt == 1) p = svx - svo;	// (MBE)

	    	// total PK-2 stress, S_n+1
	    	S = Sx - J*p*Ci;   

	    	// Extra material tangent (for isotropic growth ONLY)
	    	tens4dmm Cx11 = pow(J,1.0/3.0)*dyad1mm(Se,Ci)*(1.0/3.0);
	    	tens4dmm Cx12 = pow(J,4.0/3.0)*dyad1mm(dSedJ,Ci);
	    	tens4dmm Cx2  = pow(J,-1.0/3.0)*Ce;
		tens4dmm Cx = Cx11 + Cx12 + Cx2;	
	    	
	    	// spatial extra tangent (pushforward)
		tens4dmm cx = Cx.pp(F)/J;	

	    	// extra Cauchy stress
		mat3ds sx = 1.0/J*(F*(Sx*F.transpose())).sym();

		// volumetric spatial tangent (continuous)
		// tangent: the multiplier depends on the angle through svx, so the extra tangent in dp/dC includes calph
		tens4dmm cp = 1.0/3.0*(2.0*sx.tr()*IoI-2.0*dyad1mm(I,sx)-ddot(IxI,cx+calph))+ svo*(IxI-2.0*IoI);
			
		// total spatial tangent, c_n+1
		c = cx + cp + cpAl + calph;	// includes the variation of alpha
	}

	// debug output (pr_opt = 1)
	if (m_writeRes) {
	    	char *msg = 
	    		"F = %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"
    			"S = %f\t%f\t%f\t%f\t%f\t%f\n"
    			"c =\n"
    			"%f\t%f\t%f\t%f\t%f\t%f\n"
    			"%f\t%f\t%f\t%f\t%f\t%f\n"
    			"%f\t%f\t%f\t%f\t%f\t%f\n"
    			"%f\t%f\t%f\t%f\t%f\t%f\n"
    			"%f\t%f\t%f\t%f\t%f\t%f\n"
    			"%f\t%f\t%f\t%f\t%f\t%f\n"
    			"--------------------------------------------\n";
		feLog(msg, 
    		F(0,0),F(1,1),F(2,2),F(0,1),F(1,0),F(0,2),F(2,0),F(1,2),F(2,1),
    		S.xx(),S.yy(),S.zz(),S.xy(),S.xz(),S.yz(),
    		c(0,0,0,0),c(0,0,1,1),c(0,0,2,2),c(0,0,0,1),c(0,0,1,2),c(0,0,2,0),
    		c(1,1,0,0),c(1,1,1,1),c(1,1,2,2),c(1,1,0,1),c(1,1,1,2),c(1,1,2,0),
   		c(2,2,0,0),c(2,2,1,1),c(2,2,2,2),c(2,2,0,1),c(2,2,1,2),c(2,2,2,0),
   		c(0,1,0,0),c(0,1,1,1),c(0,1,2,2),c(0,1,0,1),c(0,1,1,2),c(0,1,2,0),
    		c(1,2,0,0),c(1,2,1,1),c(1,2,2,2),c(1,2,0,1),c(1,2,1,2),c(1,2,2,0),
    		c(2,0,0,0),c(2,0,1,1),c(2,0,2,2),c(2,0,0,1),c(2,0,1,2),c(2,0,2,0));
    	}
	return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density function
double FEKinematicGrowthZ::StrainEnergyDensity(FEMaterialPoint& mp)
{
	// get machine epsilon 
	double eps = std::numeric_limits<double>::epsilon();
	
	// Get current time, end time and time-step
    	double t = GetFEModel()->GetTime().currentTime;

	double partialtime = m_sgrEnd;			// end of the degradation ramp
	double sgr = std::min(t,partialtime);		// time of the degradation ramp
	
	// get the KinematicGrowth point data
	FEKinematicMaterialPointZ& pt = *mp.ExtractData<FEKinematicMaterialPointZ>();
    	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
    	mat3d &F = et.m_F; double J = et.m_J;

	mat3dd I(1.0);				// identity (2nd order)
	
	// Get growth deformation gradient, Fg_n+1
	mat3d  Fg = pow(J,1.0/3.0)*I; if (t <= m_sgrSta + eps) Fg = mat3d(1,0,0, 0,1,0, 0,0,1);	// isotropic-volumetric growth ONLY; F_g = 1 in Stage I
	mat3d  Fgi = Fg.inverse();		// Fg_n+1 inverse
	double Jgi = Fgi.det();		// det Fg^-1 (not used)
	
	// elastic deformation gradient, Fe_n+1 (for elastic stress calculation)
	mat3d Fe;
	
	mat3d Fo = pt.m_Fo;
	
	Fe = F*Fgi;
	if (m_lo_opt && t > m_sgrSta + eps){
		Fe = Fo;	// lo_opt = 1: F_e = F_o after gr_start
	}

	// MATERIAL PARAMETERS ----------------------------

	// mass/volume fractions
	double phieo = m_phieo;
	double phico = m_phico;
	double phimo = 1.0 - phieo - phico;
	// f_e + f_c > 1: the three mass fractions are set to 1
	if (phieo + phico > 1.0){
		phieo = 1.0;
		phico = 1.0;
		phimo = 1.0;
	}
	
	double mu  = m_mu; 
	
	double cm = m_cm;
	double dm = m_dm;
	double cc = m_cc;
	double dc = m_dc;

	double Gm = m_Gm;
	double Gc = m_Gc;
	double Get = m_Get;
	double Gez = m_Gez;
	
	double betat = m_betat;
	double betaz = m_betaz;
	
	double betad = 0.5*(1.0 - betat - betaz);
	// beta_t + beta_z > 1: beta_t = beta_z = 1, beta_d = 0.5
	if (betat + betaz > 1.0){
		betat = 1.0; 
		betaz = 1.0; 
		betad = 0.5; 
	}
	
	double alpha_d = m_alpha;	

	double aexp = 1.0;									// fibre-angle update: tan(alpha) = tan(alpha_0) (lt/lz)^aexp
	double alpha = alpha_d*M_PI/180;	// rad
	
	// retrieve material position
	vec3d X = mp.m_r0;
	vec3d Xcl = {0.0, 0.0, X.z};				// center line
	vec3d NX = {X.x-Xcl.x,X.y-Xcl.y,X.z-Xcl.z};		// radial vector
	double ro = sqrt(NX*NX);
	NX /= ro;

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// local basis: Nz axial (local system of the mesh), Nt circumferential, Nr radial
	vec3d Nz(Q(0,2),Q(1,2),Q(2,2)); 	// axial
	vec3d Nt = vec3d(-NX.y, NX.x, 0.0);  // circ
	vec3d Nr = Nz ^ Nt; 			// radial
	double lenz = Nz.norm(); double lent = Nt.norm(); double lenr = Nr.norm();
	if (lenz > 1e-12) Nz /= lenz;
	if (lent > 1e-12) Nt /= lent;
	if (lenr > 1e-12) Nr /= lenr;
	vec3d N[3] = {Nr, Nt, Nz};

    	if (t > m_sgrSta + eps && m_alpm==1) {
    	
    		// update alpha 
		// stretches measured from the original homeostatic state (F Fo^-1), not from the mesh (F)
		mat3d Fr = F * Fo.inverse();
		double lt = (Fr * N[1]).norm();
		double lz = (Fr * N[2]).norm();
		alpha = atan(tan(alpha) * pow(lt / lz, aexp));
    	   	
	    	}
    	
    	vec3d Np = Nt * sin(alpha) + Nz * cos(alpha);
    	vec3d Nn = Nt * sin(alpha) - Nz * cos(alpha);

	double azimuth = acos(-NX.y);			// azimuth wrt axis -Y

	mat3ds Ge = I;	// Ge (unit matrix)
	// Ge from spectral decomposition
	N[0] = N[2]^N[1];
	Ge = 1.0/Get/Gez*dyad(N[0]) + Get*dyad(N[1]) + Gez*dyad(N[2]);
	
	double muout = mu;
	double muin  = mu;
	
	double ccout  = cc;
	double ccin   = cc;
	
	if (t > m_sgrSta + eps && m_degr==1) {
		//	ELASTIN DAMAGE
		muout = m_muout*mu;
		muin  = m_muin*muout;
		
		if (m_asym == 0) {
			mu = muout+(muin-muout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz));	// axisymmetric aneurysm
		}
		if (m_asym == 1) {
			mu = muout+(muin-muout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz))*exp(-pow(std::abs((azimuth-M_PI)/(M_PI/m_pif)),m_nut)); 	// asymmetric aneurysm	
		}
	
		//	REDUCED COLLAGEN CROSS-LINKING
		ccout  = m_ccout*cc;
		ccin   = m_ccin*ccout;
		
		if (m_asym == 0) {
			cc = ccout+(ccin-ccout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz));	// axisymmetric aneurysm
		}
		if (m_asym == 1) {
			cc = ccout+(ccin-ccout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz))*exp(-pow(std::abs((azimuth-M_PI)/(M_PI/m_pif)),m_nut)); 	// asymmetric aneurysm	
		}
	}
	
	pt.m_elas_dmg = mu/m_mu; 
	pt.m_coll_dmg = cc/m_cc;
	
	// Calculate the new elastic strain energy density
	if (m_hypSel == 1){
		mat3ds Ce = (Fe.transpose()*Fe).sym();
		pt.m_sed = J*0.5*mu*(Ce.tr()-3.0);	// Neo-hooke: W_R = J mu/2 (tr Ce - 3), consistent with Se = mu I
	}
	if (m_hypSel == 2){				// St.Venant-Kirchhoff
		mat3ds Ae = ((Fe.transpose()*Fe-I).sym())*0.5; // Green-Lagrange elastic strain
		pt.m_sed = J*mu*Ae.dotdot(Ae);	// W_R = J mu Ae:Ae, consistent with Se = 2 mu Ae
	}
	if (m_hypSel == 3){
		
		double lto = (Fe*N[1]).norm();
		double lzo = (Fe*N[2]).norm();
		double lpo = (Fe*Np).norm();
		double lno = (Fe*Nn).norm();

		double lmt2 = (Gm*lto)*(Gm*lto);
		double lct2 = (Gc*lto)*(Gc*lto);
		double lcz2 = (Gc*lzo)*(Gc*lzo);
		double lcp2 = (Gc*lpo)*(Gc*lpo);						// deposition stretch calculation (computational purposes)
		double lcn2 = (Gc*lno)*(Gc*lno);						// idem for symmetric
	
		double sede = ( mu/2.0*((Fe*Ge).dotdot(Fe*Ge)-3.0) );

		double sedm = (cm/(4.0*dm)*(exp(dm*(lmt2-1.0)*(lmt2-1.0))-1.0) );
		double sedc =  (
				cc/(4.0*dc)*(exp(dc*(lct2-1.0)*(lct2-1.0))-1.0)*betat +
				cc/(4.0*dc)*(exp(dc*(lcz2-1.0)*(lcz2-1.0))-1.0)*betaz +
				cc/(4.0*dc)*(exp(dc*(lcp2-1.0)*(lcp2-1.0))-1.0)*betad +
				cc/(4.0*dc)*(exp(dc*(lcn2-1.0)*(lcn2-1.0))-1.0)*betad 
				);

		pt.m_sed = J*(phieo * sede + phimo * sedm + phico * sedc);
	}

    	double sed = pt.m_sed;
    
	double sedt = sed;

	// return the total strain energy density
	return sedt;
}
//-----------------------------------------------------------------------------
// Dummy symmetric tangent to satisfy C++ base class requirements
tens4ds FEKinematicGrowthZ::Tangent(FEMaterialPoint& mp)
{
	return tens4ds(0.0);
}
