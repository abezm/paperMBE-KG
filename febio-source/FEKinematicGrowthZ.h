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

	//mat3ds	m_H[MAX_TERMS];		//!< internal variables
	//mat3ds	m_Hp[MAX_TERMS];	//!< internal variables at previous timestep

    //double  m_alpha[MAX_TERMS];     //!< exponent of right-stretch tensor in series spring
    //double  m_alphap[MAX_TERMS];    //!< alpha at previous time step
    
	double	m_sed;	//!< elastic strain energy density
	double	m_sedp;	//!< elastic strain energy density at previous time
	
	double	m_Jp;	//!< J at previous time
	double	m_J;	//!< J at current time
	
	mat3d	m_Fgp;	//!< Fg at previous time
	mat3d	m_Fg;	//!< Fg at current time
	
	double	m_svo;			//!< volumetric stress at o
	
	mat3d m_Fo;
	
	double  m_elas_dmg;
	double  m_coll_dmg;
	
	double m_Wt;
};


//-----------------------------------------------------------------------------
//! This class implements a large deformation visco-elastic material
//
class FEKinematicGrowthZ :	public FEElasticMaterial
{
public:
	// NOTE: make sure that this parameter is the 
	//       same as the MAX_TERMS in the FEKinematicMaterialPointZ class
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
    
    //! calculate exponent of right-stretch tensor in series spring
    //bool SeriesStretchExponent(FEMaterialPoint& pt);
    
    // returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;
	
	// In order to use the secant tangent approximation, override the UseSecantTangent() and have it return true:
	bool UseSecantTangent() override { return m_secant_tangent; }

public: 
	// material parameters
	//double	m_g0;			//!< intitial visco-elastic coefficient
	//double	m_g[MAX_TERMS];	//!< visco-elastic coefficients
	//double	m_t[MAX_TERMS];	//!< relaxation times
	double  m_kPar;		// GAIN*RATE parameter
	double  m_sgrSta;		// initial G&R time
	double  m_sgrEnd;		// final G&R time
	double  m_kJ;
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
	//double m_Get;
	//double m_Gez;
	//
	double m_mu;			// shear modulus
	double m_cm;
	double m_dm;
	double m_cc;
	double m_dc;
	//double m_cct;
	//double m_dct;
	//double m_ccz;
	//double m_dcz;
	
	double m_alpha;
	
	double m_betat;			
	double m_betaz;
	
	double m_phieo; 
	double m_phico; 
	
	double m_z_center;	// cylinder half length
	double m_z_spread;
	double m_nuz;		// axisymmetric deviation parameter
	double m_nut;		// asymmetric deviation parameter
	double m_pif;		// asymmetric angle (factor of PI)
	
	// Active stress
	//double m_Tmax;
	//double m_lamM;
	//double m_lam0;
	//double m_cb;
	

private:
	FEElasticMaterial*	m_Base;	//!< pointer to elastic solid material

public:
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
