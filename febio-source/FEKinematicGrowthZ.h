/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material point data for kinematic growth materials
class FEKinematicMaterialPointZ : public FEMaterialPointData
{
public:
	enum { MAX_TERMS = 6 };

public:
	//! constructor
	FEKinematicMaterialPointZ(FEMaterialPointData* mp = nullptr);

	//! copy material point data
	FEMaterialPointData* Copy();

	//! Initialize material point data
	void Init();

	//! Update material point data
	void Update(const FETimeInfo& timeInfo);

	//! Serialize data to archive
	void Serialize(DumpStream& ar);

public:
	mat3ds	m_Se;	//!< elastic PK2 stress
	mat3ds	m_Sep;	//!< elastic 2nd PK stress at previous time

	double	m_sed;	//!< elastic strain energy density
	double	m_sedp;	//!< elastic strain energy density at previous time
	
	double	m_Jp;	//!< J at previous time
	double	m_J;	//!< J at current time
	
	mat3d	m_Fgp;	//!< Fg at previous time
	mat3d	m_Fg;	//!< Fg at current time
	
	double	m_svo;			//!< volumetric stress at o
	
	mat3d m_Fo;		//!< F at the end of Stage I (original homeostatic state)
	
	double  m_elas_dmg;	//!< degraded c_e / input c_e
	double  m_coll_dmg;	//!< degraded c_c / input c_c
	
	double m_Wt;		//!< m_sed / J
};

//-----------------------------------------------------------------------------
//! Kinematic growth material: equilibrated (MBE) or rate-dependent isotropic growth with degradation
class FEKinematicGrowthZ :	public FEElasticMaterial
{
public:
	enum { MAX_TERMS = FEKinematicMaterialPointZ::MAX_TERMS };

public:
	//! default constructor
	FEKinematicGrowthZ(FEModel* pfem);

	//! get the elastic base material
	FEElasticMaterial* GetBaseMaterial();

	//! Set the base material
	void SetBaseMaterial(FEElasticMaterial* pbase);

public:
	//! stress function
	mat3ds Stress(FEMaterialPoint& pt) override;

	//! tangent function
	tens4ds Tangent(FEMaterialPoint& pt) override;

	//! tangent function (only minor symm)
	tens4dmm TangentMm(FEMaterialPoint& pt) override;

	//! strain energy density
	double StrainEnergyDensity(FEMaterialPoint& pt) override;

    // returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;
	
	// In order to use the secant tangent approximation, override the UseSecantTangent() and have it return true:
	bool UseSecantTangent() override { return m_secant_tangent; }

public: 
	// material parameters
	double  m_kPar;		// rate parameter K_sigma k
	double  m_sgrSta;		// end of Stage I, start of G&R
	double  m_sgrEnd;		// end of the degradation ramp
	double  m_kJ;		// volumetric penalty of Stage I
	bool    m_MBE_opt;		// option for MBE
	bool    m_secant_tangent;  	// flag for using secant tangent calculation
	bool	m_asym;		// type of aneurysm
	int     m_hypSel;		// select hyperelastic model
	bool	m_degr;		// elastin & collagen degradation
	bool    m_alpm;		// update alpha
	double  m_muin;
	double  m_ccin;
	double  m_muout;
	double  m_ccout;
	
	double  m_Gm;
	double  m_Gc;
	double  m_Get;
	double  m_Gez;
	bool    m_lo_opt;		// option for stretch calc (1: Fe = Fo)
	bool    m_writeRes;		// option for printing
	double m_mu;			// shear modulus
	double m_cm;
	double m_dm;
	double m_cc;
	double m_dc;
	
	double m_alpha;
	
	double m_betat;			
	double m_betaz;
	
	double m_phieo; 
	double m_phico; 
	
	double m_z_center;	// axial centre of the degradation profile
	double m_z_spread;	// axial spread of the degradation profile
	double m_nuz;		// axial exponent of the degradation profile
	double m_nut;		// circumferential exponent of the degradation profile
	double m_pif;		// circumferential spread pi/t_pif

private:
	FEElasticMaterial*	m_Base;	//!< pointer to elastic solid material

public:
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
