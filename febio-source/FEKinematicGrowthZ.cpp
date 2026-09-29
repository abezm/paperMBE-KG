#include "stdafx.h"
#include "FEKinematicGrowthZ.h"
#include "FEUncoupledMaterial.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/DumpStream.h>
#include <FECore/log.h>
#define _USE_MATH_DEFINES				// to introduce pi constant (1/2)
//#include <math.h>					// to introduce pi constant (2/2)
#include <cmath>					// Replaced math.h with cmath for std::abs
#include "FECore/log.h"				// to print to log file and/or screen
#include <limits>
#include <algorithm>                    		// Added for std::min

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEKinematicGrowthZ, FEElasticMaterial)

	// 1. define the material properties
	ADD_PROPERTY(m_Base, "elastic");
	
	// ##############  Hyperelastic model ##############
	ADD_PARAMETER(m_mu    , FE_RANGE_GREATER_OR_EQUAL(0.0), "c_e"); 
	ADD_PARAMETER(m_cm    , FE_RANGE_GREATER_OR_EQUAL(0.0), "c_m");
	ADD_PARAMETER(m_dm    , FE_RANGE_GREATER(0.0), "d_m");
	ADD_PARAMETER(m_cc    , FE_RANGE_GREATER_OR_EQUAL(0.0), "c_c");
	ADD_PARAMETER(m_dc    , FE_RANGE_GREATER(0.0), "d_c");
	//ADD_PARAMETER(m_cct   , FE_RANGE_GREATER_OR_EQUAL(0.0), "c_ct");
	//ADD_PARAMETER(m_dct   , FE_RANGE_GREATER(0.0), "d_ct");
	//ADD_PARAMETER(m_ccz   , FE_RANGE_GREATER_OR_EQUAL(0.0), "c_cz");
	//ADD_PARAMETER(m_dcz   , FE_RANGE_GREATER(0.0), "d_cz");
	ADD_PARAMETER(m_alpha , "alpha");
	
	// ############## mass fractions ###################
	
	// constituent mass fractions (elastin, collagen)
	ADD_PARAMETER(m_phieo, FE_RANGE_CLOSED(0.0,1.0), "f_e"); 
	ADD_PARAMETER(m_phico, FE_RANGE_CLOSED(0.0,1.0), "f_c"); 
	
	// collagen relative fractions 
	ADD_PARAMETER(m_betaz, FE_RANGE_CLOSED(0.0,1.0), "beta_z"); 
	ADD_PARAMETER(m_betat, FE_RANGE_CLOSED(0.0,1.0), "beta_t"); 
	
	
	// ##############  Kinematic Growth model ##############
	ADD_PARAMETER(m_kPar, FE_RANGE_GREATER(0.0), "k_parameter");			// 2. GAIN*RATE parameter
	ADD_PARAMETER(m_sgrSta, FE_RANGE_GREATER_OR_EQUAL(0.0), "gr_start");		// 3. initial G&R time
	ADD_PARAMETER(m_sgrEnd, FE_RANGE_GREATER_OR_EQUAL(0.0), "gr_end");		// FINAL G&R time
	
	ADD_PARAMETER(m_MBE_opt,  "MBE");						// 4. option for MBE G&R
	ADD_PARAMETER(m_asym,     "asymm_aneurysm");					// 5. type of aneurysm
	ADD_PARAMETER(m_lo_opt,   "lo_opt");						// option for stretch calc
	ADD_PARAMETER(m_writeRes, "pr_opt");						// option for printing 
	
	ADD_PARAMETER(m_kJ, "k_J"); // penalizer (incompress)
	//ADD_PARAMETER(m_mu, "C_1");  // shear modulus

	// ############## deposition stretches ##############
	ADD_PARAMETER(m_Get, FE_RANGE_GREATER_OR_EQUAL(0.0), "Ge_t");
	ADD_PARAMETER(m_Gez, FE_RANGE_GREATER_OR_EQUAL(0.0), "Ge_z");
	ADD_PARAMETER(m_Gm,  FE_RANGE_GREATER_OR_EQUAL(0.0), "Gm");
	ADD_PARAMETER(m_Gc,  FE_RANGE_GREATER_OR_EQUAL(0.0), "Gc");
	
	ADD_PARAMETER(m_muout, FE_RANGE_CLOSED(0.0,1.0), "mu_out"); 
	ADD_PARAMETER(m_muin,  FE_RANGE_CLOSED(0.0,1.0), "mu_in"); 
	ADD_PARAMETER(m_ccout, FE_RANGE_CLOSED(0.0,1.0), "cc_out"); 
	ADD_PARAMETER(m_ccin,  FE_RANGE_CLOSED(0.0,1.0), "cc_in");

	ADD_PARAMETER(m_hypSel, "hyperels_mod"); // 6. hyperelastic model option 
	ADD_PARAMETER(m_degr,   "ec_degrad"); // 7. properties are degraded w this option 
	ADD_PARAMETER(m_alpm,   "alp_update"); // 8. alpha is updated w this option 
	
	// ############## secant tangent option ##############
	ADD_PARAMETER(m_secant_tangent, "secant_tangent");	// 8. secant tangent
	
	// ############## geometry of cylinder (for degradation) ##############
	ADD_PARAMETER(m_z_center, "z_center");
	ADD_PARAMETER(m_z_spread, "z_spread");
	ADD_PARAMETER(m_nuz, "z_nu");
	ADD_PARAMETER(m_nut, "t_nu");
	ADD_PARAMETER(m_pif, "t_pif");
	
	// ## ACTIVE STRESS ##
	//ADD_PARAMETER(m_Tmax, FE_RANGE_GREATER_OR_EQUAL(0.0), "T_max");
	//ADD_PARAMETER(m_lamM, FE_RANGE_GREATER(0.0), "lam_M"); // 
	//ADD_PARAMETER(m_lam0, FE_RANGE_GREATER(0.0), "lam_0"); // 
	//ADD_PARAMETER(m_cb,   FE_RANGE_GREATER(0.0), "c_b");   // 

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
	//ar & m_H & m_Hp;
    ar & m_sed & m_sedp;
    //ar & m_alpha & m_alphap;
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
	
	//double endtime = GetFEModel()->GetCurrentStep()->m_tend;

	//double endtime = 11.0;				// 11.0 | 31.0-32.0 (TEVG)
	double partialtime = m_sgrEnd; //endtime;			// partialtime <= endtime | 10.0 | 10.4 (for TI calculation)
	double sgr = std::min(t,partialtime);		// min(t,partialtime) | min(t,9.0)
	
	
	if (dt == 0) return mat3ds(0, 0, 0, 0, 0, 0);
    
	// get the elastic part
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the KinematicGrowth point data
	FEKinematicMaterialPointZ& pt = *mp.ExtractData<FEKinematicMaterialPointZ>();
	
	// Get Jacobian of deformation gradient, J_n (previous time-step)
	double Jp = pt.m_Jp;	// Fg_n
	
	// Get total deformation gradient, F_n+1
	mat3d F = ep.m_F;		// F_n+1
	double J = ep.m_J;		// J_n+1 
	pt.m_J = J;			// saves J_n+1 to later update Jp at the end of current timestep
	
	// Get growth deformation gradient, Fg_n+1
	mat3dd I(1.0);				// identity (2nd order)
	mat3d  Fg = pow(J,1.0/3.0)*I;		// isotropic-volumetric growth ONLY
	mat3d  Fgi = Fg.inverse();		// Fg_n+1 inverse
	double Jgi = Fgi.det();		// 1/J_n
	pt.m_Fg = Fg;				// saves Fg_n+1 to later update Jp at the end of current timestep
	
	// elastic deformation gradient, Fe_n+1 (for elastic stress calculation)
	mat3d Fe; // = et.m_F*Fgi;		// Note: px.m_Fe is Fe_n
	//double Je = et.m_J*Jgi;
	
	mat3d Fo = pt.m_Fo;
	
	Fe = F*Fgi;		// 
	if (m_lo_opt && t > m_sgrSta + eps){
		Fe = Fo;	// elastic strain does not evolve during G&R of aneurysms?
	}
	
	// =======================================
	
	// substitute elastic deformation in material point
	//ep.m_F = Fe;
	//ep.m_J = Je;

	// ------------------------------------------------
	// MATERIAL PARAMETERS ----------------------------
	//double mu = m_mu; //100.0;	// NEO-HOOKE / ST.VENANT-K (UNIAXIAL EXAMPLE)

	// mass/volume fractions
	///*
	double phieo = m_phieo;
	double phico = m_phico;
	double phimo = 1.0 - phieo - phico;
	if (phieo + phico > 1.0){
		phieo = 1.0; //0.34;
		phico = 1.0; //0.33;
		phimo = 1.0; //
	}
	//*/
	//double phieo = 1.0; //0.34;
	//double phico = 1.0; //0.33;
	//double phimo = 1.0; //
	//
	//
	
	// --- TO PRINT ---
	//feLog("Mass fractions: Elastin(phieo)=%f, Collagen(phico)=%f, SMC(phimo)=%f\n", phieo, phico, phimo);
	// -------------------------


	double mu  = m_mu; 
	
	double cm = m_cm;
	double dm = m_dm;
	double cc = m_cc;
	double dc = m_dc;
	
	//double ccz = m_ccz;
	//double dcz = m_dcz;
	//double cct = m_cct;
	//double dct = m_dct;
	
	double Gm = m_Gm;
	double Gc = m_Gc;
	double Get = m_Get;
	double Gez = m_Gez;
	
	double betat = m_betat;		// 0.045
	double betaz = m_betaz;		// 0.062
	
	double betad = 0.5*(1.0 - betat - betaz);
	if (betat + betaz > 1.0){
		betat = 1.0; 
		betaz = 1.0; 
		betad = 0.5; 
	}
	
	double alpha_d = m_alpha;
	
	// ======================================
	/*
	// KNOCKOUTS
	//double phieo = 0.34;
	//double phico = 0.33;
	//double phimo = 1.0 - phieo - phico;
	
	double mu = 89.71; 
	double Gm = m_Gm; //1.20;  // dep stretch smc
	double Gc = m_Gc; //1.25; // dep stretch col
	//double Get = 1.90; //m_Get;		
	//double Gez = 1.62; //m_Gez;
	
	// ---------------
	double cm = 261.4;	// kPa
	double dm = 0.24;
	double cc = 234.9;	// kPa
	double dc = 4.08;
	double alpha_d = 29.9;	// deg
	//
	double betat = 0.056;	// orientation fractions for collagen
	double betaz = 0.067;
	*/
	//======================================
	/*
	// TORTUOSITY
	//double phieo = 0.331;
	//double phimo = 0.265;
	//double phico = 0.404;

	double mu = 88.46;									// 88.46

	double alpha_d = 38.1; //0.665;							// 0.615 (110) / 0.665 (120)

	double cm = 419.1;									// 419.1
	double dm = 0.105;									// 0.105
	double Gm = m_Gm; //1.198;									// 1.198
	double cc = 1278.;									// 1278.
	double dc = 11.10;									// 11.10
	double Gc = m_Gc; //1.092;								// 1.092
	
	//double Get = 2.153;			// 1.993 (110) / 2.153 (120)
	//double Gez = 1.572;			// 1.572

	double betat = 0.045;		// 0.045
	double betaz = 0.062;		// 0.062
	*/
	
	// ======================================
	//
	/*	
	// ELASTIN/COLL. REDUCTION FACTORS
	double m_muin; double m_ccin;
	if (m_asym == 1) {	// A-symmetric aneurysms
		m_muin = 0.50; 
		m_ccin = 0.40;	
	}
	if (m_asym == 0) {	// AXI-symmetric aneurysms
		m_muin = 0.70; 
		m_ccin = 0.70;	
	
	}
	*/
	///*
	//double m_muout = 1.00;		// 1.00 (alone) | 1.00-0.30 = 0.70 (aging) | 1.00 (hypertension) | 1.00 (alone & contractility (50))
	//double m_muin  = 0.40;		// 0.38 (alone) | 0.68-0.30 = 0.38 (aging) | 0.63 (hypertension) | 0.30 (alone & contractility (50))
	
	//double m_ccout = 1.00;
	//double m_ccin  = 0.60;		// 0.77 (alone) | 0.91 (aging) | 0.91 (hypertension)
	//*/
	// ------------------------------------------------
	double aexp = 1.0;									// 1.0 (KNOCKOUTS) | 0.0 (CMAME and TORTUOSITY)
	double alpha = alpha_d*M_PI/180;	// rad
	
	// retrieve material position
	vec3d X = mp.m_r0;
	vec3d Xcl = {0.0, 0.0, X.z};				// center line
	vec3d NX = {X.x-Xcl.x,X.y-Xcl.y,X.z-Xcl.z};		// radial vector
	double ro = sqrt(NX*NX);
	NX /= ro;

	// retrieve local element basis directions
	//vec3d N[3];

	// pointwise, consistent with mesh generated with Matlab script
	//N[2] = {0.0, 0.0, 1.0}; N[2] = N[2]/sqrt(N[2]*N[2]);	// axial = d(Xcl)/d(z)
	//N[1] = {-NX.y, NX.x, NX.z};				// circu.
	//N[0] = N[2]^N[1];					// radial
	//N[2] = {1.0, 0.0, 0.0}; N[2] = N[2]/sqrt(N[2]*N[2]);	// axial = d(Xcl)/d(z)
	//N[1] = {0.0, 1.0, 0.0}; N[1] = N[1]/sqrt(N[1]*N[1]);	// circumferential
	
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// Unified fiber directions 
	vec3d Nz(Q(0,2),Q(1,2),Q(2,2)); 	// axial
	//vec3d Nt(Q(0,1),Q(1,1),Q(2,1)); 	// circ
	vec3d Nt = vec3d(-NX.y, NX.x, 0.0);  // circ
	vec3d Nr = Nz ^ Nt; 			// radial
	double lenz = Nz.norm(); double lent = Nt.norm(); double lenr = Nr.norm();
	if (lenz > 1e-12) Nz /= lenz;
	if (lent > 1e-12) Nt /= lent;
	if (lenr > 1e-12) Nr /= lenr;
	vec3d N[3] = {Nr, Nt, Nz};
	
	/*
	if (m_degr==0) {
    		// Copy the local element basis directions to n
		//vec3d n[2];
    		//n[0].x = Q[0][0]; n[0].y = Q[1][0]; n[0].z = Q[2][0];	// axial? 	
    		//n[1].x = Q[0][1]; n[1].y = Q[1][1]; n[1].z = Q[2][1];	// circ.?
    		N[2].x = Q[0][0]; N[2].y = Q[1][0]; N[2].z = Q[2][0];	// axial? 	
    		N[1].x = Q[0][1]; N[1].y = Q[1][1]; N[1].z = Q[2][1];	// circ.?
    		Np = N[2]*cos(alpha) + N[1]*sin(alpha); 
		Nn = N[2]*cos(alpha) - N[1]*sin(alpha); 
    	}
    	*/
    	double lt = 1.;
    	double lz = 1.;
    	if (t > m_sgrSta + eps && m_alpm==1) { // && m_muin < 1.0) {
    	
    		// update alpha
		lt = (F * N[1]).norm();
		lz = (F * N[2]).norm();
		alpha = atan(tan(alpha) * pow(lt / lz, aexp));
    	
		// pointwise, consistent with mesh generated with Matlab script
		/*
		N[2] = {0.0, 0.0, 1.0}; N[2] = N[2]/sqrt(N[2]*N[2]);	// axial = d(Xcl)/d(z)
		N[1] = {-NX.y, NX.x, NX.z};				// circu.
		Np = N[1]*sin(alpha)+N[2]*cos(alpha);	// original diagonal fiber dir.
		Nn = N[1]*sin(alpha)-N[2]*cos(alpha);	// idem for symmetric
		*/
    	}
    	
    	vec3d Np = Nt * sin(alpha) + Nz * cos(alpha);
    	vec3d Nn = Nt * sin(alpha) - Nz * cos(alpha);
    	
    	// --- CÓDIGO DE DEPURACIÓN PARA IMPRIMIR ---
	// Imprime los vectores base
	//feLog("N[1] (Nt) = [%.5f, %.5f, %.5f]\n", Nt.x, Nt.y, Nt.z);
	//feLog("N[2] (Nz) = [%.5f, %.5f, %.5f]\n", Nz.x, Nz.y, Nz.z);
	
	// Imprime el ángulo y estiramientos
	//feLog("lt = %.5f, lz = %.5f, alpha_act = %.5f rad\n", lt, lz, alpha_act);
	
	// Imprime las fibras resultantes
	//feLog("Np = [%.5f, %.5f, %.5f]\n", Np.x, Np.y, Np.z);
	//feLog("Nn = [%.5f, %.5f, %.5f]\n", Nn.x, Nn.y, Nn.z);
	//feLog("--------------------------------------------------\n");
    	
    	// Evaluate the structural direction in the reference configuration
	//vec3d  Np = N[1]*sin(alpha)+N[2]*cos(alpha);	// original diagonal fiber dir.
	//vec3d  Nn = N[1]*sin(alpha)-N[2]*cos(alpha);	// idem for symmetric
	//vec3d  Np = N[2]*cos(alpha) + N[1]*sin(alpha); 
	//vec3d  Nn = N[2]*cos(alpha) - N[1]*sin(alpha); 
	
	//feLog("axial(N2)=%f\t%f\t%f\n",N[2].x,N[2].y,N[2].z);
	//feLog("circ.(N1)=%f\t%f\t%f\n",N[1].x,N[1].y,N[1].z);
	
	
	double azimuth = acos(-NX.y);			// azimuth wrt axis -Y
	
	mat3ds Ge = I;	// Ge (unit matrix)
	///*
	//if (m_degr==1) {
	// Ge from spectral decomposition
	//N[0] = N[2]^N[1];
	Ge = 1.0/Get/Gez*dyad(N[0]) + Get*dyad(N[1]) + Gez*dyad(N[2]);
	//} 
	//*/


	double muout = mu;
	double muin  = mu;
	
	double ccout  = cc;
	double ccin   = cc;		// 0.77 (alone) | 0.91 (aging) | 0.91 (hypertension)
		
	
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
		ccin   = m_ccin*ccout;		// 0.77 (alone) | 0.91 (aging) | 0.91 (hypertension)
		//double cctout = m_ccout*cct;
		//double cctin  = m_ccin*cct;		
		//double cczout = m_ccout*ccz;
		//double cczin  = m_ccin*ccz;		// 
		
		if (m_asym == 0) {
			cc = ccout+(ccin-ccout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz));	// axisymmetric aneurysm
			//cct = cctout+(ccin-cctout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),5));	// axisymmetric aneurysm
			//ccz = cczout+(ccin-cczout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),5));	// axisymmetric aneurysm
		}
		if (m_asym == 1) {
			cc = ccout+(ccin-ccout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz))*exp(-pow(std::abs((azimuth-M_PI)/(M_PI/m_pif)),m_nut)); 	// asymmetric aneurysm	
			//cct = cctout+(ccin-cctout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),5))*exp(-pow(std::abs((azimuth-M_PI)/(M_PI/3.0)),5)); 
			//ccz = cczout+(ccin-cczout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),5))*exp(-pow(std::abs((azimuth-M_PI)/(M_PI/3.0)),5)); 
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
		
		//int x_opt = 1;
		/*
		if (m_lo_opt && t <= m_sgrSta + eps){
			lto = (F*N[1]).norm();
			lzo = (F*N[2]).norm();
			lpo = (F*Np).norm();
			lno = (F*Nn).norm();
			// hacer G = 1?? REV; MARCOS NO LO HACE
			// G > 1 uniaxial leads to SDB!?
			//Gm = 1.0;
			//Gc = 1.0;
		}
		if (m_lo_opt && t > m_sgrSta + eps){
			//mat3d Fo = pt.m_Fo;
			lto = (Fo*N[1]).norm();
			lzo = (Fo*N[2]).norm();
			lpo = (Fo*Np).norm();
			lno = (Fo*Nn).norm();
		}
		*/
		
		lmt2 = (Gm*lto)*(Gm*lto);
		lct2 = (Gc*lto)*(Gc*lto);
		lcz2 = (Gc*lzo)*(Gc*lzo);
		lcp2 = (Gc*lpo)*(Gc*lpo);						// deposition stretch calculation (computational purposes)
		lcn2 = (Gc*lno)*(Gc*lno);						// idem for symmetric
		
		mat3ds Sel = mu*Ge.sqr();	
		/*
		mat3ds Sm = (cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm)*dyad(N[1]));
		mat3ds Sc = (cc*(lct2-1.0)*exp(dc*(lct2-1.0)*(lct2-1.0))*(Gc*Gc)*dyad(N[1])*betat +
				cc*(lcz2-1.0)*exp(dc*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc)*dyad(N[2])*betaz +
				cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc)*dyad( Np )*betad +
				cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc)*dyad( Nn )*betad );
		*/
		/*
		double Shm  = cm*(lt2-1.0)*exp(dm*(lt2-1.0)*(lt2-1.0));
		
		double Shct = cc*(lt2-1.0)*exp(dc*(lt2-1.0)*(lt2-1.0));
		double Shcz = cc*(lz2-1.0)*exp(dc*(lz2-1.0)*(lz2-1.0));
		double Shcp = cc*(lp2-1.0)*exp(dc*(lp2-1.0)*(lp2-1.0));
		double Shcn = cc*(ln2-1.0)*exp(dc*(ln2-1.0)*(ln2-1.0));
		*/
		
		double Shm  = cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm);
		
		//double Shct = cct*(lct2-1.0)*exp(dct*(lct2-1.0)*(lct2-1.0))*(Gc*Gc);
		//double Shcz = ccz*(lcz2-1.0)*exp(dcz*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc);
		double Shct = cc*(lct2-1.0)*exp(dc*(lct2-1.0)*(lct2-1.0))*(Gc*Gc);
		double Shcz = cc*(lcz2-1.0)*exp(dc*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc);
		double Shcp = cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc);
		double Shcn = cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc);
		
		mat3ds Sm = Shm*dyad(N[1]);
		mat3ds Sc = (  Shct*dyad(N[1])*betat + Shcz*dyad(N[2])*betaz +
				Shcp*dyad( Np )*betad + Shcn*dyad( Nn )*betad );
		//
				
		Se = phieo * Sel + phimo * Sm + phico * Sc;	// elastin + SMCs + collagen
		
	}
	
	// restore safe copy
	//ep.m_F = F;
	//ep.m_J = J;
	
	//char *msg = "F = %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"		"--------------------------------------------\n";
	//feLog(msg,F(0,0),F(1,1),F(2,2),F(0,1),F(1,0),F(0,2),F(2,0),F(1,2),F(2,1));

	//char *msgx = "X = %f\t%f\t%f\n";
	//feLog(msgx,X.x,X.y,X.z);
	
	// =======================================
	
	// active ("relaxed" / MBE-only: see knockouts)
	//double Tmax = 250.0 * 0.0;			// 250.0 | 50.0 | 150.0 (for uniform cases, except for contractility -> 250)
	//double lamM = 1.1;
	//double lam0 = 0.4;
	//double CB = sqrt(log(2.0));			// such that (1-exp(-CB^2)) = 0.5 (0.8326)
	//double CS = 0.5*m_cb * 1.0;			// such that (1-exp( -C^2)) = 0.0 for lt = 1/(1+CB/CS)^(1/3) = 0.7 and (1-exp(-C^2)) = 0.75 for lt = 2.0
	//mat3ds Sao = phimo * m_Tmax*(1.0-exp(-m_cb*m_cb))*(1.0-pow((m_lamM-1.0)/(m_lamM-m_lam0),2))*(lto*lto)*dyad(Nt);
	
	// extra PK-2 stress, Sx_n+1
	// isotropic growth ONLY
	//mat3ds Sx = pow(J,1.0/3.0)*Se + phimo * Sa; 
	//mat3ds Sx = pow(J,1.0/3.0)*(Se + phimo * Sa); // if Sa is a function of elastic strain
	mat3ds Sx = pow(J,1.0/3.0) * Se; // TODO add active stress
	
	// Cauchy-Green strain
	mat3ds C = ep.RightCauchyGreen();
	mat3ds Ci = C.inverse();
	
	// volumetric extra Cauchy stress, svx_n+1
	double svx = 1.0/3.0/J*Sx.dotdot(C);
	
	mat3ds S(0.0);
	double svo = 0.0;
	
	double Jdep = 0.9999;
	double lm = m_kJ; //1.0e5;
	
	if (t <= m_sgrSta + eps) {
    	
		// total Cauchy stress, s_n+1
	    	S = Se + Ci*lm*log(Jdep*J);
	    	//S = Ci*lm*log(Jdep*J);	
	    	
	    	// store volumetric Cauchy stress  (original-homeostatic)
	    	svo   = 1.0/3.0/J*S.dotdot(C);	 //(1.0/3.0)*se.dotdot(I);
	    	pt.m_svo = svo;
	    	
	    	// store Fo
	    	pt.m_Fo = F;
    		
    	}
    	else {
	    	// retrieve volumetric Cauchy stress  (original-homeostatic)
	    	svo = pt.m_svo;
	    	
	    	// Lagrange multiplier, p_n+1
	    	// Note: stimulus func. for arteries ONLY; f(sigx_1)
	    	double p = 0.0;
	    	if (m_MBE_opt == 0){ // (RATE DEPENDENT)
	    		p = svx - svo*(1.0+log(J/Jp)/(m_kPar*dt));	
		}
		if (m_MBE_opt == 1) p = svx - svo;	// (MBE) // p = 1.0/3.0/J*Sx.dotdot(C) - svo;
			//p = svx - svo/(1.0-delta)*(1.0+KsKi*(EPS*pow(rIrIo,-3)-1.0)-KfKi*inflam);		// Ups = 1 -> p
		
	    	// total PK-2 stress, S_n+1
	    	S = Sx - J*p*Ci;   
	    	
	    	// prints stuff
		//int m_writeRes = 0;
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
	//mat3ds sigma = ep.push_forward(S);
	mat3ds sigma = 1.0/J*(F*(S*F.transpose())).sym();
	
	// prints stuff
	//int m_writeRes = 0;
	if (m_writeRes) {
	    	char *msg = "F = %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"
			//"Fo = %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"
    			//"S = %f\t%f\t%f\t%f\t%f\t%f\n"
    			"sig = %f\t%f\t%f\t%f\t%f\t%f\n"
    			//"lsq = %f\t%f\t%f\t%f\t%f\n"
    			"lo = %f\t%f\t%f\t%f\n"
    			"Np = %f\t%f\t%f\n"
    			"Nn = %f\t%f\t%f\n"
    			"Nz = %f\t%f\t%f\n"
    			"Nc = %f\t%f\t%f\n"
    			"--------------------------------------------\n";
		feLog(msg,
    		F(0,0),F(1,1),F(2,2),F(0,1),F(1,0),F(0,2),F(2,0),F(1,2),F(2,1),
		//Fo(0,0),Fo(1,1),Fo(2,2),Fo(0,1),Fo(1,0),Fo(0,2),Fo(2,0),Fo(1,2),Fo(2,1), 
    		//S.xx(),S.yy(),S.zz(),S.xy(),S.xz(),S.yz());
    		sigma.xx(),sigma.yy(),sigma.zz(),sigma.xy(),sigma.xz(),sigma.yz(),
    		lto, lzo, lpo, lno,
    		//lmt2, lct2, lcz2, lcp2, lcn2,
    		Np(0), Np(1), Np(2),
    		Nn(0), Nn(1), Nn(2),
    		N[2].x, N[2].y, N[2].z,
    		N[1].x, N[1].y, N[1].z
    		);
    		// 	print(lm2, lct2, lcz2, lcp2, lcn2)
    	}

	//feLog("sig = %f\t%f\t%f\t%f\t%f\t%f\n",sigma.xx(),sigma.yy(),sigma.zz(),sigma.xy(),sigma.xz(),sigma.yz());
	
	return sigma; //ep.push_forward(S);
	//return s;
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
	
	//double endtime = GetFEModel()->GetCurrentStep()->m_tend;
	//double endtime = GetFEModel()->GetCurrentStep()->m_tend;

	//double endtime = 11.0;				// 11.0 | 31.0-32.0 (TEVG)  // end of time-interval for stimulus application?
	double partialtime = m_sgrEnd;			// partialtime <= endtime | 10.0 | 10.4 (for TI calculation)
	double sgr = std::min(t,partialtime);		// min(t,partialtime) | min(t,9.0)
 
 	// Get elastic material point data (def grad, F_n+1)
 	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
 	
 	// get the KinematicGrowth point data 
	FEKinematicMaterialPointZ& pt = *mp.ExtractData<FEKinematicMaterialPointZ>();
 	
 	// Get Jacobian of deformation gradient, J_n (previous time-step)
 	double Jp = pt.m_Jp;	// Fg_n
 	
 	// Get total deformation gradient, F_n+1
	mat3d &F = ep.m_F;		// F_n+1
	double J = ep.m_J;		// J_n+1 
	//pt.m_J = J;			// already saved J_n+1; 
	
	// Dyadic products of identity
	mat3dd I(1.0);				// identity (2nd order)
	tens4dmm IoI = dyad4mm(I,I);
	tens4dmm IxI = dyad1mm(I,I);
	
	// Get growth deformation gradient, Fg_n+1
	mat3d  Fg = pow(J,1.0/3.0)*I;		// isotropic-volumetric growth ONLY
	mat3d  Fgi = Fg.inverse();		// Fg_n+1 inverse
	double Jgi = Fgi.det();		// 1/J_n
	
	// elastic deformation gradient, Fe_n+1 (for elastic stress calculation)
	mat3d Fe; // = et.m_F*Fgi;		// Note: px.m_Fe is Fe_n
	//double Je = et.m_J*Jgi;
	
	mat3d Fo = pt.m_Fo;
	
	Fe = F*Fgi;		// 
	if (m_lo_opt && t > m_sgrSta + eps){
		Fe = Fo;	// elastic strain does not evolve during G&R of aneurysms?
	}
	
	// =======================================
	
	// substitute elastic deformation in material point
	//ep.m_F = Fe;
	//ep.m_J = Je;


	// MATERIAL PARAMETERS ----------------------------
	//double mu = m_mu; //100.0;	// NEO-HOOKE / ST.VENANT-K (UNIAXIAL EXAMPLE)

	// mass/volume fractions
	///*
	double phieo = m_phieo;
	double phico = m_phico;
	double phimo = 1.0 - phieo - phico;
	if (phieo + phico > 1.0){
		phieo = 1.0; //0.34;
		phico = 1.0; //0.33;
		phimo = 1.0; //
	}
	//*/
	//double phieo = 1.0; //0.34;
	//double phico = 1.0; //0.33;
	//double phimo = 1.0; //
	//
	//


	double mu  = m_mu; 
	
	double cm = m_cm;
	double dm = m_dm;
	double cc = m_cc;
	double dc = m_dc;
	
	//double ccz = m_ccz;
	//double dcz = m_dcz;
	//double cct = m_cct;
	//double dct = m_dct;
	
	double Gm = m_Gm;
	double Gc = m_Gc;
	double Get = m_Get;
	double Gez = m_Gez;
	
	double betat = m_betat;		// 0.045
	double betaz = m_betaz;		// 0.062
	
	double betad = 0.5*(1.0 - betat - betaz);
	if (betat + betaz > 1.0){
		betat = 1.0; 
		betaz = 1.0; 
		betad = 0.5; 
	}
	
	double alpha_d = m_alpha;
	
	// ======================================
	/*
	// KNOCKOUTS
	//double phieo = 0.34;
	//double phico = 0.33;
	//double phimo = 1.0 - phieo - phico;
	
	double mu = 89.71; 
	double Gm = m_Gm; //1.20; //m_Gm; // dep stretch smc
	double Gc = m_Gc; //1.25; //m_Gc; // dep stretch col
	//double Get = 1.90; //m_Get;		
	//double Gez = 1.62; //m_Gez;
	
	// ---------------
	double cm = 261.4;	// kPa
	double dm = 0.24;
	double cc = 234.9;	// kPa
	double dc = 4.08;
	double alpha_d = 29.9;	// deg
	//
	double betat = 0.056;	// orientation fractions for collagen
	double betaz = 0.067;
	
	*/
	// ======================================
	/*
	// TORTUOSITY
	//double phieo = 0.331;
	//double phimo = 0.265;
	//double phico = 0.404;

	double mu = 88.46;									// 88.46

	double alpha_d = 38.1; //0.665;							// 0.615 (110) / 0.665 (120)

	double cm = 419.1;									// 419.1
	double dm = 0.105;									// 0.105
	double Gm = m_Gm; //1.198;									// 1.198
	double cc = 1278.;									// 1278.
	double dc = 11.10;									// 11.10
	double Gc = m_Gc; //1.092;									// 1.092
	
	//double Get = 2.153;			// 1.993 (110) / 2.153 (120)
	//double Gez = 1.572;			// 1.572

	double betat = 0.045;		// 0.045
	double betaz = 0.062;		// 0.062
	*/
	// ======================================
	//
	/*	
	// ELASTIN/COLL. REDUCTION FACTORS
	double m_muin; double m_ccin;
	if (m_asym == 1) {	// A-symmetric aneurysms
		m_muin = 0.50; 
		m_ccin = 0.40;	
	}
	if (m_asym == 0) {	// AXI-symmetric aneurysms
		m_muin = 0.70; 
		m_ccin = 0.70;	
	
	}
	*/
	///*
	//double m_muout = 1.00;		// 1.00 (alone) | 1.00-0.30 = 0.70 (aging) | 1.00 (hypertension) | 1.00 (alone & contractility (50))
	//double m_muin  = 0.40;		// 0.38 (alone) | 0.68-0.30 = 0.38 (aging) | 0.63 (hypertension) | 0.30 (alone & contractility (50))
	
	//double m_ccout = 1.00;
	//double m_ccin  = 0.60;		// 0.77 (alone) | 0.91 (aging) | 0.91 (hypertension)
	//*/
	//
	// ------------------------------------------------
	double aexp = 1.0;			// 1.0 (KNOCKOUTS) | 0.0 (CMAME and TORTUOSITY)
	double alpha = alpha_d*M_PI/180;	// rad
	
	
	// retrieve material position
	vec3d X = mp.m_r0;
	vec3d Xcl = {0.0, 0.0, X.z};				// center line
	vec3d NX = {X.x-Xcl.x,X.y-Xcl.y,X.z-Xcl.z};		// radial vector
	double ro = sqrt(NX*NX);
	NX /= ro;

	// retrieve local element basis directions
	//vec3d N[3];
	
	// pointwise, consistent with mesh generated with Matlab script
	//N[2] = {0.0, 0.0, 1.0}; N[2] = N[2]/sqrt(N[2]*N[2]);	// axial = d(Xcl)/d(z)
	//N[1] = {-NX.y, NX.x, NX.z};				// circu.
	//N[0] = N[2]^N[1];					// radial
	//N[2] = {1.0, 0.0, 0.0}; N[2] = N[2]/sqrt(N[2]*N[2]);	// axial = d(Xcl)/d(z)
	//N[1] = {0.0, 1.0, 0.0}; N[1] = N[1]/sqrt(N[1]*N[1]);	// circumferential
	
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);
	
	// Unified fiber directions 
	vec3d Nz(Q(0,2),Q(1,2),Q(2,2)); 	// axial
	//vec3d Nt(Q(0,1),Q(1,1),Q(2,1)); 	// circ
	vec3d Nt = vec3d(-NX.y, NX.x, 0.0);  // circunferencial (tangente al radio)
	vec3d Nr = Nz ^ Nt; 			// radial
	double lenz = Nz.norm(); double lent = Nt.norm(); double lenr = Nr.norm();
	if (lenz > 1e-12) Nz /= lenz;
	if (lent > 1e-12) Nt /= lent;
	if (lenr > 1e-12) Nr /= lenr;
	vec3d N[3] = {Nr, Nt, Nz};

	mat3ds tent = dyad(Nt);
	mat3ds tenz = dyad(Nz);


	/*
	if (m_degr==0) {
    		// Copy the local element basis directions to n
		//vec3d n[2];
    		//n[0].x = Q[0][0]; n[0].y = Q[1][0]; n[0].z = Q[2][0];	// axial? 	
    		//n[1].x = Q[0][1]; n[1].y = Q[1][1]; n[1].z = Q[2][1];	// circ.?
    		N[2].x = Q[0][0]; N[2].y = Q[1][0]; N[2].z = Q[2][0];	// axial? 	
    		N[1].x = Q[0][1]; N[1].y = Q[1][1]; N[1].z = Q[2][1];	// circ.?
    		Np = N[2]*cos(alpha) + N[1]*sin(alpha); 
		Nn = N[2]*cos(alpha) - N[1]*sin(alpha); 
    	}
    	*/
    	mat3ds dtanadC(0.0);
    	vec3d dNpdta(0.0);
    	vec3d dNndta(0.0);
    	
    	if (t > m_sgrSta + eps && m_alpm==1) {
    		
    		// update alpha
		double lt = (F * Nt).norm();
		double lz = (F * Nz).norm();
		alpha = atan(tan(alpha) * pow(lt / lz, aexp));
    	
		// pointwise, consistent with mesh generated with Matlab script
		/*
		N[2] = {0.0, 0.0, 1.0}; N[2] = N[2]/sqrt(N[2]*N[2]);	// axial = d(Xcl)/d(z)
		N[1] = {-NX.y, NX.x, NX.z};				// circu.
		Np = N[1]*sin(alpha)+N[2]*cos(alpha);	// original diagonal fiber dir.
		Nn = N[1]*sin(alpha)-N[2]*cos(alpha);	// idem for symmetric
		*/
		
		// derivative of the fiber direction vectors Np and Nn with respect to tan(alpha)
		dNpdta = (Nt - Nz*tan(alpha))*pow(1 + pow(tan(alpha),2),-1.5);	// d(Np)/d(tan(alpha))
		dNndta = (Nt + Nz*tan(alpha))*pow(1 + pow(tan(alpha),2),-1.5);
		
		// derivative 2 * d(tan(alpha)) / dC
		dtanadC = (tent * (1.0 / (lt*lt)) - tenz * (1.0 / (lz*lz))) * (tan(alpha)*aexp);	// 2*d(tan(alpha))/d(C) : (Ft)o(Ft)
    	}
    	
    	vec3d Np = Nt * sin(alpha) + Nz * cos(alpha);
    	vec3d Nn = Nt * sin(alpha) - Nz * cos(alpha);
    	
    	mat3ds tenp = dyad(Np);
	mat3ds tenn = dyad(Nn);
    	
    	// Evaluate the structural direction in the reference configuration
	//vec3d  Np = N[1]*sin(alpha)+N[2]*cos(alpha);	// original diagonal fiber dir.
	//vec3d  Nn = N[1]*sin(alpha)-N[2]*cos(alpha);	// idem for symmetric
	//vec3d  Np = N[2]*cos(alpha) + N[1]*sin(alpha); 
	//vec3d  Nn = N[2]*cos(alpha) - N[1]*sin(alpha); 
	
	double azimuth = acos(-NX.y);			// azimuth wrt axis -Y
	
	
	mat3ds Ge = I;	// Ge (unit matrix)
	///*
	//if (m_degr==1) {
	// Ge from spectral decomposition
	N[0] = N[2]^N[1];
	Ge = 1.0/Get/Gez*dyad(N[0]) + Get*tent + Gez*tenz;
	//} 
	//*/
	
		
	
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
		ccin   = m_ccin*ccout;		// 0.77 (alone) | 0.91 (aging) | 0.91 (hypertension)
		//double cctout = m_ccout*cct;
		//double cctin  = m_ccin*cct;		
		//double cczout = m_ccout*ccz;
		//double cczin  = m_ccin*ccz;		// 
		
		if (m_asym == 0) {
			cc = ccout+(ccin-ccout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz));	// axisymmetric aneurysm
			//cct = cctout+(ccin-cctout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),5));	// axisymmetric aneurysm
			//ccz = cczout+(ccin-cczout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),5));	// axisymmetric aneurysm
		}
		if (m_asym == 1) {
			cc = ccout+(ccin-ccout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz))*exp(-pow(std::abs((azimuth-M_PI)/(M_PI/m_pif)),m_nut)); 	// asymmetric aneurysm	
			//cct = cctout+(ccin-cctout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),5))*exp(-pow(std::abs((azimuth-M_PI)/(M_PI/3.0)),5)); 
			//ccz = cczout+(ccin-cczout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),5))*exp(-pow(std::abs((azimuth-M_PI)/(M_PI/3.0)),5)); 
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
		//mat3ds CGei = CGe.inverse();
		mat3ds Ae = (CGe-I)*0.5; 			   // Green-Lagrange elastic strain
		//mat3ds Ae = ((Fe.transpose()*Fe-I).sym())*0.5; // Green-Lagrange elastic strain
		Se = Ae*2.0*mu;	
		
		Ce = dyad4s(I,I)*mu; // St.Venant-Kirchhoff
		
		dSedJ = (-2.0*mu/3.0)*pow(J,-5.0/3.0)*(F.transpose()*F).sym();

	}
	if (m_hypSel == 3){		// 4-FIBER FAMILY MOD. 
		
		double lto = (Fe*N[1]).norm();
		double lzo = (Fe*N[2]).norm();
		double lpo = (Fe*Np).norm();
		double lno = (Fe*Nn).norm();
		
		//int x_opt = 1;
		/*
		if (m_lo_opt && t <= m_sgrSta + eps){
			lto = (F*N[1]).norm();
			lzo = (F*N[2]).norm();
			lpo = (F*Np).norm();
			lno = (F*Nn).norm();
			// hacer G = 1?? REV; MARCOS NO LO HACE
			// G > 1 uniaxial leads to SDB!?
			//Gm = 1.0;
			//Gc = 1.0;
		}
		if (m_lo_opt && t > m_sgrSta + eps){
			//mat3d Fo = pt.m_Fo;
			lto = (Fo*N[1]).norm();
			lzo = (Fo*N[2]).norm();
			lpo = (Fo*Np).norm();
			lno = (Fo*Nn).norm();
		}
		*/
		
		double lmt2 = (Gm*lto)*(Gm*lto);
		double lct2 = (Gc*lto)*(Gc*lto);
		double lcz2 = (Gc*lzo)*(Gc*lzo);
		double lcp2 = (Gc*lpo)*(Gc*lpo);						// deposition stretch calculation (computational purposes)
		double lcn2 = (Gc*lno)*(Gc*lno);						// idem for symmetric

		
		mat3ds Sel = mu*Ge.sqr();	
		
		double Shm  = cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm);
		
		//double Shct = cct*(lct2-1.0)*exp(dct*(lct2-1.0)*(lct2-1.0))*(Gc*Gc);
		//double Shcz = ccz*(lcz2-1.0)*exp(dcz*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc);
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
		
		//double dShct = cct*(1.0+2.0*dct*(lct2-1.0)*(lct2-1.0))*exp(dct*(lct2-1.0)*(lct2-1.0))*pow(Gc,4);
		//double dShcz = ccz*(1.0+2.0*dcz*(lcz2-1.0)*(lcz2-1.0))*exp(dcz*(lcz2-1.0)*(lcz2-1.0))*pow(Gc,4);
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
		
		// Tangent stiffness (var. of alpha)
		tens4dmm Calphp = dyad1mm(dSpdta, dtanadC);
		tens4dmm Calphn = dyad1mm(dSndta, dtanadC);
		Calph = phico * betad * (Calphp + Calphn);
		
		Ce = phimo * Cm + phico * Cc;	// SMCs + collagen
		
		dSedJ = (-2.0/(3.0*J))*(phimo * dShm*lmt2*dyad(N[1]) + 
			 phico * dShct*lct2*dyad(N[1])*betat + phico * dShcz*lcz2*dyad(N[2])*betaz + 
			 phico * dShcp*lcp2*dyad( Np )*betad + phico * dShcn*lcn2*dyad( Nn )*betad);
		
		///*	 
		if(m_degr==1 && t > m_sgrSta + eps) {	// constant Fo (aneurysms)
			Ce = Ce*0.0;
			dSedJ = dSedJ*0.0;	// [REV; commmented it does not converge]
		}  
		//*/
		
	}
	
	// restore safe copy
	//ep.m_F = F;
	//ep.m_J = J;
	
	// spatial elastic tangent
	//tens4dmm ce = m_Base->Tangent(mp);	// no funciona?
	tens4dmm ce = Ce.pp(F)/J;
	//tens4dmm Ce = ce.pp(F.inverse())*J;// pull-back to get material tangent
	
	calph = Calph.pp(F)/J;
	
	if (dt == 0.0) return ce;
	// restore safe copy
	//ep.m_F = F;
	//ep.m_J = J;
	
	
	// =======================================
	
	// EXTRA PK-2 stress, Sx_n+1
	mat3ds Sx = pow(J,1.0/3.0)*Se; // isotropic growth ONLY
	
	// Cauchy-Green strain
	mat3ds C = ep.RightCauchyGreen();
	mat3ds Ci = C.inverse();

	tens4dmm CixCi = dyad1mm(Ci,Ci);	// 	rev
	tens4dmm CioCi = dyad4mm(Ci,Ci);	// debe ser el simetrico: 1/2 ...etc  (Holzapfel, p.255)
	//tens4ds CixCi = dyad1s(Ci);
	//tens4ds CioCi = dyad4s(Ci);
	
	// volumetric extra Cauchy stress, svx_n+1
	double svx = 1.0/3.0/J*Sx.dotdot(C);
	
	tens4dmm c(0.0);
	tens4dmm Cp(0.0);
	tens4dmm cp(0.0);
	mat3ds S(0.0);
	
	// retrieves volumetric Cauchy stress  (original-homeostatic)
	double svo = pt.m_svo;
	//int m_writeRes = 0;
	if (t <= m_sgrSta + eps) {
		//m_writeRes = 1; // to write tangent
    
	    	double Jdep = 0.9999;
		double lm = m_kJ; //1.0e5;
	    	
	    	// material volum tangent
	    	Cp = lm*(CixCi-2.0*log(Jdep*J)*CioCi);	// p = dW_vol/dC = d/dC (ln J)**2
	    	
	    	// spatial volum tangent
	    	//cp = lm/J*(IxI-2.0*log(Jdep*J)*IoI);
	    	cp = Cp.pp(F)/J;
	    	
	    	// total spatial tangent, c_n+1
	    	c = ce + cp + calph;		// [XXX] added contrib from var. of alpha
	    	
	    	// total PK-2 stress, S_n+1 / only to print
	    	S = Se + Ci*lm*log(Jdep*J);	
    	}
    	
    	else {
   		//m_writeRes = 1; // to write tangent
	    	
	    	// update Lagrange multiplier
	    	// Note: stimulus func. for arteries ONLY; f(sigx_1)
	    	double p = 0.0;
	    	tens4dmm cpAl(0.0);		// Spatial Algorithmic tangent
	    	if (m_MBE_opt == 0){ // (RATE DEPENDENT)
	    	
	    		p = svx - svo*(1.0+log(J/Jp)/(m_kPar*dt));
	    		
	    		// Algorithmic volumetric tangent, (dp_n+1_alg)/(dC)
			//tens4ds CixCi = dyad1s(Ci);
			//tens4ds CioCi = dyad4s(Ci);
	    		tens4dmm CpAl = (J*svo/(m_kPar*dt))*((1.0+log(J/Jp))*CixCi-2.0*log(J/Jp)*CioCi);
	    		cpAl = CpAl.pp(F)/J;
		}
		if (m_MBE_opt == 1) p = svx - svo;	// (MBE) // p = 1.0/3.0/J*Sx.dotdot(C) - svo;
	    	
		
	    	// total PK-2 stress, S_n+1
	    	S = Sx - J*p*Ci;   
	    	//S = Se; // test
	    	//mat3ds S = Sx; // test
	 
	    	
	    	// Extra material tangent (for isotropic growth ONLY)
	    	tens4dmm Cx11 = pow(J,1.0/3.0)*dyad1mm(Se,Ci)*(1.0/3.0);
	    	tens4dmm Cx12 = pow(J,4.0/3.0)*dyad1mm(dSedJ,Ci);
	    	tens4dmm Cx2  = pow(J,-1.0/3.0)*Ce;
	    	//
		//tens4dmm Cx = pow(J,1.0/3.0)*(1.0/3.0)*dyad1mm(Se,Ci) + pow(J,-1.0/3.0)*Ce;
		tens4dmm Cx = Cx11 + Cx12 + Cx2;	
	    	
	    	// spatial extra tangent (pushforward)
		tens4dmm cx = Cx.pp(F)/J;	

	    	// extra Cauchy stress
		mat3ds sx = 1.0/J*(F*(Sx*F.transpose())).sym();

		// volumetric spatial tangent (continuous)
		tens4dmm cp = 1.0/3.0*(2.0*sx.tr()*IoI-2.0*dyad1mm(I,sx)-ddot(IxI,cx))+ svo*(IxI-2.0*IoI);
			
		// total spatial tangent, c_n+1
		c = cx + cp + cpAl + calph;	// [XXX] added contrib from var. of alpha
		//c = ce;	// to test
		//c = ce + cp;	// to test
		//c = cx;	// to test
	}
	    	
	//feLog("mu/cc = %f\t%f\n",mu,cc);

	// prints stuff
	//feLog("alpha = %f\n", alpha);
	if (m_writeRes) {
	    	char *msg = 
	    		//"X = %f\t%f\t%f\n"
	    		//"svo = %f\n"
	    		"F = %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"
			//"Fo = %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"
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
		//X.x,X.y,X.z, svo,
    		F(0,0),F(1,1),F(2,2),F(0,1),F(1,0),F(0,2),F(2,0),F(1,2),F(2,1),
		//Fo(0,0),Fo(1,1),Fo(2,2),Fo(0,1),Fo(1,0),Fo(0,2),Fo(2,0),Fo(1,2),Fo(2,1), 
    		S.xx(),S.yy(),S.zz(),S.xy(),S.xz(),S.yz(), //); //,
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

	//double dt = GetTimeInfo().timeIncrement;
	
	//double endtime = GetFEModel()->GetCurrentStep()->m_tend;
	//double endtime = GetFEModel()->GetCurrentStep()->m_tend;

	//double endtime = 11.0;				// 11.0 | 31.0-32.0 (TEVG)
	double partialtime = m_sgrEnd;			// partialtime <= endtime | 10.0 | 10.4 (for TI calculation)
	double sgr = std::min(t,partialtime);		// min(t,partialtime) | min(t,9.0)
	
	// get the KinematicGrowth point data
	FEKinematicMaterialPointZ& pt = *mp.ExtractData<FEKinematicMaterialPointZ>();
    	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
    	mat3d &F = et.m_F; double J = et.m_J;

	mat3dd I(1.0);				// identity (2nd order)
	
	// Get growth deformation gradient, Fg_n+1
	mat3d  Fg = pow(J,1.0/3.0)*I;		// isotropic-volumetric growth ONLY
	mat3d  Fgi = Fg.inverse();		// Fg_n+1 inverse
	double Jgi = Fgi.det();		// 1/J_n
	
	// elastic deformation gradient, Fe_n+1 (for elastic stress calculation)
	mat3d Fe; // = et.m_F*Fgi;		// Note: px.m_Fe is Fe_n
	//double Je = et.m_J*Jgi;
	
	mat3d Fo = pt.m_Fo;
	
	Fe = F*Fgi;		// 
	if (m_lo_opt && t > m_sgrSta + eps){
		Fe = Fo;	// elastic strain does not evolve during G&R of aneurysms?
	}

	// substitute elastic deformation in material point
	//et.m_F = Fe;
	//et.m_J = Je;

	// MATERIAL PARAMETERS ----------------------------
	//double mu = m_mu; //100.0;	// NEO-HOOKE / ST.VENANT-K (UNIAXIAL EXAMPLE)

	// mass/volume fractions
	///*
	double phieo = m_phieo;
	double phico = m_phico;
	double phimo = 1.0 - phieo - phico;
	if (phieo + phico > 1.0){
		phieo = 1.0; //0.34;
		phico = 1.0; //0.33;
		phimo = 1.0; //
	}
	//*/
	//double phieo = 1.0; //0.34;
	//double phico = 1.0; //0.33;
	//double phimo = 1.0; //
	//
	
	double mu  = m_mu; 
	
	double cm = m_cm;
	double dm = m_dm;
	double cc = m_cc;
	double dc = m_dc;
	
	//double ccz = m_ccz;
	//double dcz = m_dcz;
	//double cct = m_cct;
	//double dct = m_dct;
	
	double Gm = m_Gm;
	double Gc = m_Gc;
	double Get = m_Get;
	double Gez = m_Gez;
	
	double betat = m_betat;		// 0.045
	double betaz = m_betaz;		// 0.062
	
	double betad = 0.5*(1.0 - betat - betaz);
	if (betat + betaz > 1.0){
		betat = 1.0; 
		betaz = 1.0; 
		betad = 0.5; 
	}
	
	double alpha_d = m_alpha;	
	
	// ======================================
	/*
	// KNOCKOUTS
	//double phieo = 0.34;
	//double phico = 0.33;
	//double phimo = 1.0 - phieo - phico;
	
	double mu = 89.71; 
	double Gm = m_Gm; //1.20; // dep stretch smc
	double Gc = m_Gc; //1.25; // dep stretch col
	//double Get = 1.90; //m_Get;		
	//double Gez = 1.62; //m_Gez;
	
	// ---------------
	double cm = 261.4;	// kPa
	double dm = 0.24;
	double cc = 234.9;	// kPa
	double dc = 4.08;
	double alpha_d = 29.9;	// deg
	//
	double betat = 0.056;	// orientation fractions for collagen
	double betaz = 0.067;
	*/ 
	//======================================
	/*
	// TORTUOSITY
	//double phieo = 0.331;
	//double phimo = 0.265;
	//double phico = 0.404;

	double mu = 88.46;									// 88.46

	double alpha_d = 38.1; //0.665;							// 0.615 (110) / 0.665 (120)

	double cm = 419.1;									// 419.1
	double dm = 0.105;									// 0.105
	double Gm = m_Gm; //1.198;									// 1.198
	double cc = 1278.;									// 1278.
	double dc = 11.10;									// 11.10
	double Gc = m_Gc; //1.092;									// 1.092
	
	double Get = 2.153;			// 1.993 (110) / 2.153 (120)
	double Gez = 1.572;			// 1.572

	double betat = 0.045;		// 0.045
	double betaz = 0.062;		// 0.062
	*/
	// ======================================
	//
	/*	
	// ELASTIN/COLL. REDUCTION FACTORS
	double m_muin; double m_ccin;
	if (m_asym == 1) {	// A-symmetric aneurysms
		m_muin = 0.50; 
		m_ccin = 0.40;	
	}
	if (m_asym == 0) {	// AXI-symmetric aneurysms
		m_muin = 0.70; 
		m_ccin = 0.70;	
	
	}
	*/
	///*
	//double m_muout = 1.00;		// 1.00 (alone) | 1.00-0.30 = 0.70 (aging) | 1.00 (hypertension) | 1.00 (alone & contractility (50))
	//double m_muin  = 0.40;		// 0.38 (alone) | 0.68-0.30 = 0.38 (aging) | 0.63 (hypertension) | 0.30 (alone & contractility (50))
	
	//double m_ccout = 1.00;
	//double m_ccin  = 0.60;		// 0.77 (alone) | 0.91 (aging) | 0.91 (hypertension)
	//*/
	// ------------------------------------------------
	double aexp = 1.0;									// 1.0 (KNOCKOUTS) | 0.0 (CMAME and TORTUOSITY)
	double alpha = alpha_d*M_PI/180;	// rad
	
	// retrieve material position
	vec3d X = mp.m_r0;
	vec3d Xcl = {0.0, 0.0, X.z};				// center line
	vec3d NX = {X.x-Xcl.x,X.y-Xcl.y,X.z-Xcl.z};		// radial vector
	double ro = sqrt(NX*NX);
	NX /= ro;

	// retrieve local element basis directions
	//vec3d N[3];

	// pointwise, consistent with mesh generated with Matlab script
	//N[2] = {0.0, 0.0, 1.0}; N[2] = N[2]/sqrt(N[2]*N[2]);	// axial = d(Xcl)/d(z)
	//N[1] = {-NX.y, NX.x, NX.z};				// circumferential
	//N[0] = N[2]^N[1];					// radial
	//N[2] = {1.0, 0.0, 0.0}; N[2] = N[2]/sqrt(N[2]*N[2]);	// axial = d(Xcl)/d(z)
	//N[1] = {0.0, 1.0, 0.0}; N[1] = N[1]/sqrt(N[1]*N[1]);	// circumferential
	
	
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// Unified fiber directions 
	vec3d Nz(Q(0,2),Q(1,2),Q(2,2)); 	// axial
	//vec3d Nt(Q(0,1),Q(1,1),Q(2,1)); 	// circ
	vec3d Nt = vec3d(-NX.y, NX.x, 0.0);  // circunferencial (tangente al radio)
	vec3d Nr = Nz ^ Nt; 			// radial
	double lenz = Nz.norm(); double lent = Nt.norm(); double lenr = Nr.norm();
	if (lenz > 1e-12) Nz /= lenz;
	if (lent > 1e-12) Nt /= lent;
	if (lenr > 1e-12) Nr /= lenr;
	vec3d N[3] = {Nr, Nt, Nz};
	

	/*
	if (m_degr==0) {
    		// Copy the local element basis directions to n
		//vec3d n[2];
    		//n[0].x = Q[0][0]; n[0].y = Q[1][0]; n[0].z = Q[2][0];	// axial? 	
    		//n[1].x = Q[0][1]; n[1].y = Q[1][1]; n[1].z = Q[2][1];	// circ.?
    		N[2].x = Q[0][0]; N[2].y = Q[1][0]; N[2].z = Q[2][0];	// axial? 	
    		N[1].x = Q[0][1]; N[1].y = Q[1][1]; N[1].z = Q[2][1];	// circ.?
    		Np = N[2]*cos(alpha) + N[1]*sin(alpha); 
		Nn = N[2]*cos(alpha) - N[1]*sin(alpha); 
    	}
    	*/
    	if (t > m_sgrSta + eps && m_alpm==1) {
    	
    		// update alpha 
    		// [TODO] UPDATE TANGENT TO ACCOUNT FOR VARIATION OF ALPHA
		double lt = (F * N[1]).norm();
		double lz = (F * N[2]).norm();
		alpha = atan(tan(alpha) * pow(lt / lz, aexp));
    	   	
		// pointwise, consistent with mesh generated with Matlab script
		/*
		N[2] = {0.0, 0.0, 1.0}; N[2] = N[2]/sqrt(N[2]*N[2]);	// axial = d(Xcl)/d(z)
		N[1] = {-NX.y, NX.x, NX.z};				// circu.
		Np = N[1]*sin(alpha)+N[2]*cos(alpha);	// original diagonal fiber dir.
		Nn = N[1]*sin(alpha)-N[2]*cos(alpha);	// idem for symmetric
		*/
    	}
    	
    	vec3d Np = Nt * sin(alpha) + Nz * cos(alpha);
    	vec3d Nn = Nt * sin(alpha) - Nz * cos(alpha);
    	
    	// Evaluate the structural direction in the reference configuration
	//vec3d  Np = N[1]*sin(alpha)+N[2]*cos(alpha);	// original diagonal fiber dir.
	//vec3d  Nn = N[1]*sin(alpha)-N[2]*cos(alpha);	// idem for symmetric
	//vec3d  Np = N[2]*cos(alpha) + N[1]*sin(alpha); 
	//vec3d  Nn = N[2]*cos(alpha) - N[1]*sin(alpha); 
	
	double azimuth = acos(-NX.y);			// azimuth wrt axis -Y
	
	
	mat3ds Ge = I;	// Ge (unit matrix)
	///*
	//if (m_degr==1) {
	// Ge from spectral decomposition
	N[0] = N[2]^N[1];
	Ge = 1.0/Get/Gez*dyad(N[0]) + Get*dyad(N[1]) + Gez*dyad(N[2]);
	//} 
	//*/
	
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
		ccin   = m_ccin*ccout;		// 0.77 (alone) | 0.91 (aging) | 0.91 (hypertension)
		//double cctout = m_ccout*cct;
		//double cctin  = m_ccin*cct;		
		//double cczout = m_ccout*ccz;
		//double cczin  = m_ccin*ccz;		// 
		
		if (m_asym == 0) {
			cc = ccout+(ccin-ccout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz));	// axisymmetric aneurysm
			//cct = cctout+(ccin-cctout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),5));	// axisymmetric aneurysm
			//ccz = cczout+(ccin-cczout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),5));	// axisymmetric aneurysm
		}
		if (m_asym == 1) {
			cc = ccout+(ccin-ccout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),m_nuz))*exp(-pow(std::abs((azimuth-M_PI)/(M_PI/m_pif)),m_nut)); 	// asymmetric aneurysm	
			//cct = cctout+(ccin-cctout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),5))*exp(-pow(std::abs((azimuth-M_PI)/(M_PI/3.0)),5)); 
			//ccz = cczout+(ccin-cczout)*(sgr-m_sgrSta)/(m_sgrEnd-m_sgrSta)*exp(-pow(std::abs((X.z-m_z_center)/m_z_spread),5))*exp(-pow(std::abs((azimuth-M_PI)/(M_PI/3.0)),5)); 
		}
	}
	
	pt.m_elas_dmg = mu/m_mu; 
	pt.m_coll_dmg = cc/m_cc;
	
	// Calculate the new elastic strain energy density
	if (m_hypSel == 1){
		mat3ds Ce = (Fe.transpose()*Fe).sym();
		pt.m_sed = (Ce.dotdot(Ce)-3.0)*mu;	// Neo-hooke
	}
	if (m_hypSel == 2){				// St.Venant-Kirchhoff
		mat3ds Ae = ((Fe.transpose()*Fe-I).sym())*0.5; // Green-Lagrange elastic strain
		pt.m_sed = Ae.dotdot(Ae)*mu;	
	}
	if (m_hypSel == 3){
		
		double lto = (Fe*N[1]).norm();
		double lzo = (Fe*N[2]).norm();
		double lpo = (Fe*Np).norm();
		double lno = (Fe*Nn).norm();
		
		//int lo_opt = 1;
		/*
		if (m_lo_opt && t <= m_sgrSta + eps){
			lto = (F*N[1]).norm();
			lzo = (F*N[2]).norm();
			lpo = (F*Np).norm();
			lno = (F*Nn).norm();
			// hacer G = 1?? REV; MARCOS NO LO HACE
			// G > 1 uniaxial leads to SDB!?
			//Gm = 1.0;
			//Gc = 1.0;
		}
		if (m_lo_opt && t > m_sgrSta + eps){
			//mat3d Fo = pt.m_Fo;
			lto = (Fo*N[1]).norm();
			lzo = (Fo*N[2]).norm();
			lpo = (Fo*Np).norm();
			lno = (Fo*Nn).norm();
		}
		*/
		
		double lmt2 = (Gm*lto)*(Gm*lto);
		double lct2 = (Gc*lto)*(Gc*lto);
		double lcz2 = (Gc*lzo)*(Gc*lzo);
		double lcp2 = (Gc*lpo)*(Gc*lpo);						// deposition stretch calculation (computational purposes)
		double lcn2 = (Gc*lno)*(Gc*lno);						// idem for symmetric
	
		double sede = ( mu/2.0*((Fe*Ge).dotdot(Fe*Ge)-3.0) );

		double sedm = (cm/(4.0*dm)*(exp(dm*(lmt2-1.0)*(lmt2-1.0))-1.0) );
		double sedc =  (
				//cct/(4.0*dct)*(exp(dct*(lct2-1.0)*(lct2-1.0))-1.0)*betat +
				//ccz/(4.0*dcz)*(exp(dcz*(lcz2-1.0)*(lcz2-1.0))-1.0)*betaz +
				cc/(4.0*dc)*(exp(dc*(lct2-1.0)*(lct2-1.0))-1.0)*betat +
				cc/(4.0*dc)*(exp(dc*(lcz2-1.0)*(lcz2-1.0))-1.0)*betaz +
				cc/(4.0*dc)*(exp(dc*(lcp2-1.0)*(lcp2-1.0))-1.0)*betad +
				cc/(4.0*dc)*(exp(dc*(lcn2-1.0)*(lcn2-1.0))-1.0)*betad 
				);

		pt.m_sed = J*(phieo * sede + phimo * sedm + phico * sedc);
	}
	
	
    	double sed = pt.m_sed;
    
	double sedt = sed;
	
	//pt.m_Wt = sed/J;	// strain energy density per unit CURRENT volume
    
    	// restore safe copy
    	//et.m_F = F; et.m_J = J;

	// return the total strain energy density
	return sedt;
}
//-----------------------------------------------------------------------------
// Dummy symmetric tangent to satisfy C++ base class requirements
tens4ds FEKinematicGrowthZ::Tangent(FEMaterialPoint& mp)
{
	return tens4ds(0.0);
}
