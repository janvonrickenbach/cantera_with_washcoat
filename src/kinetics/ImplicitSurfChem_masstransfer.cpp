/**
 *  @file ImplicitSurfChem.cpp
 * Definitions for the implicit integration of surface site density equations
 *  (see \ref  kineticsmgr and class
 *  \link Cantera::ImplicitSurfChem ImplicitSurfChem\endlink).
 */
// Copyright 2001  California Institute of Technology

#include "cantera/kinetics/ImplicitSurfChem_masstransfer.h"
#include "cantera/numerics/Integrator.h"
#include "cantera/kinetics/solveSP.h"

using namespace std;

namespace Cantera
{

// Constructor
ImplicitSurfChem_masstransfer::ImplicitSurfChem_masstransfer(vector<InterfaceKinetics*> k) :
    FuncEval(),
    m_nsurf(0),
    m_nv(0),
    m_numBulkPhases(0),
    m_numTotalBulkSpecies(0),
    m_numTotalSpecies(0),
    m_integ(0),
    m_atol(1.e-25),
    m_rtol(1.e-7),
    m_maxstep(0.0),
    m_mediumSpeciesStart(-1),
    m_bulkSpeciesStart(-1),
    m_surfSpeciesStart(-1),
    m_surfSolver(0),
    m_commonTempPressForPhases(true),
    m_ioFlag(0)
{
    m_nsurf = k.size();
    size_t ns, nsp;
    size_t  ntmax = 0;
    size_t kinSpIndex = 0;
    // Loop over the number of surface kinetics objects
    for (size_t n = 0; n < m_nsurf; n++) {
        InterfaceKinetics* kinPtr = k[n];
        m_vecKinPtrs.push_back(kinPtr);
        ns = k[n]->surfacePhaseIndex();
        if (ns == npos)
            throw CanteraError("ImplicitSurfChem",
                               "kinetics manager contains no surface phase");
        m_surfindex.push_back(ns);
        m_surf.push_back((SurfPhase*)&k[n]->thermo(ns));
        nsp = m_surf.back()->nSpecies();
        m_nsp_surf.push_back(nsp);
        m_nsp_tot.push_back(k[n]->nTotalSpecies());
        m_nv += m_nsp_tot.back();
        if (m_nsp_tot.back() > ntmax) {
            ntmax = m_nsp_tot.back();
        }
        m_specStartIndex.push_back(kinSpIndex);
        kinSpIndex += m_nsp_tot.back();

        size_t nPhases = kinPtr->nPhases();
        vector_int pLocTmp(nPhases);
        size_t imatch = npos;
        for (size_t ip = 0; ip < nPhases; ip++) {
            if (ip != ns) {
                ThermoPhase* thPtr = & kinPtr->thermo(ip);
                 
                m_bulk_massfraction.push_back(new double[thPtr->nSpecies()]);
                m_diff_coeffs.push_back(new double[thPtr->nSpecies()]);
                thPtr->getMassFractions(m_bulk_massfraction.back());
                if ((imatch = checkMatch(m_bulkPhases, thPtr)) == npos) {
                    m_bulkPhases.push_back(thPtr);
                    m_numBulkPhases++;
                    nsp = thPtr->nSpecies();
                    m_nspBulkPhases.push_back(nsp);
                    m_numTotalBulkSpecies += nsp;
                    imatch = m_bulkPhases.size() - 1;
                }
                pLocTmp[ip] = int(imatch);
            } else {
                pLocTmp[ip] = -int(n);
            }
        }
        pLocVec.push_back(pLocTmp);

    }
    m_numTotalSpecies = m_nv + m_numTotalBulkSpecies;
    m_concSpecies.resize(m_numTotalSpecies, 0.0);
    m_concSpeciesSave.resize(m_numTotalSpecies, 0.0);

    m_integ = newIntegrator("CVODE");


    // use backward differencing, with a full Jacobian computed
    // numerically, and use a Newton linear iterator

    m_integ->setMethod(BDF_Method);
    m_integ->setProblemType(DENSE + NOJAC);
    m_integ->setIterator(Newton_Iter);
    m_work.resize(ntmax);
}

int ImplicitSurfChem_masstransfer::checkMatch(std::vector<ThermoPhase*> m_vec, ThermoPhase* thPtr)
{
    int retn = -1;
    for (int i = 0; i < (int) m_vec.size(); i++) {
        ThermoPhase* th = m_vec[i];
        if (th == thPtr) {
            return i;
        }
    }
    return retn;
}

/*
 * Destructor. Deletes the integrator.
 */
ImplicitSurfChem_masstransfer::~ImplicitSurfChem_masstransfer()
{
    delete m_integ;
    delete m_surfSolver;
}

// overloaded method of FuncEval. Called by the integrator to
// get the initial conditions.
void ImplicitSurfChem_masstransfer::getInitialConditions(doublereal t0, size_t lenc,
        doublereal* c)
{
    size_t loc = 0;
    for (size_t n = 0; n < m_nsurf; n++) {
        m_surf[n]->getCoverages(c + loc);
        m_bulkPhases[n]->getMassFractions(c + loc + m_nsp_surf[n]);
        loc += m_nsp_tot[n];
    }
}


/*
 *  Must be called before calling method 'advance'
 */
void ImplicitSurfChem_masstransfer::initialize(doublereal t0)
{
    m_integ->setTolerances(m_rtol, m_atol);
    m_integ->initialize(t0, *this);
}

/*
 *  Set masstransfer coefficient
 */
void ImplicitSurfChem_masstransfer::set_masstransfer_coefficient(doublereal h,doublereal wc_geo_area)
{
    m_masstransfer_coefficient = h;
    m_wc_geo_area = wc_geo_area;
}

void ImplicitSurfChem_masstransfer::set_transport(Transport* t)
{
    m_transport = t;
}

// Integrate from t0 to t1. The integrator is reinitialized first.
/*
 *   This routine does a time accurate solve from t = t0 to t = t1.
 *   of the surface problem.
 *
 *  @param t0  Initial Time -> this is an input
 *  @param t1  Final Time -> This is an input
 */
void ImplicitSurfChem_masstransfer::integrate(doublereal t0, doublereal t1)
{
    m_integ->initialize(t0, *this);
    m_integ->setMaxStepSize(t1 - t0);
    m_integ->integrate(t1,t0);
    updateState(m_integ->solution());
}

// Integrate from t0 to t1 without reinitializing the integrator.
/*
 *  Use when the coverages have not changed from
 *  their values on return from the last call to integrate or
 *  integrate0.
 *
 *  @param t0  Initial Time -> this is an input
 *  @param t1  Final Time -> This is an input
 */
void ImplicitSurfChem_masstransfer::integrate0(doublereal t0, doublereal t1)

{
    m_integ->integrate(t1,t0);
    updateState(m_integ->solution());
}

void ImplicitSurfChem_masstransfer::updateState(doublereal* vecConcSpecies)
{
    size_t kstart;
    for (size_t ip = 0; ip < m_nsurf; ip++) {
        kstart = m_specStartIndex[ip];
        m_surf[ip]->setCoverages(vecConcSpecies + kstart);

        ThermoPhase* TP_ptr_bulk = m_bulkPhases[ip];
        kstart = kstart+m_nsp_surf[ip];
        TP_ptr_bulk->setMassFractions(vecConcSpecies + kstart);
    }

}

/*
 * Called by the integrator to evaluate ydot given y at time 'time'.
 */
void ImplicitSurfChem_masstransfer::eval(doublereal time, doublereal* y,
                            doublereal* ydot, doublereal* p)
{
    updateState(y);   // synchronize the surface state(s) with y
    doublereal rs0, sum;
    size_t loc, kstart;
    doublereal site_density_save = 0.0;
    double Sc;
    for (size_t n = 0; n < m_nsurf; n++) {
        site_density_save = m_surf[n]->siteDensity();
        m_surf[n]->setSiteDensity(site_density_save*m_wc_geo_area);
        rs0 = 1.0/(m_surf[n]->siteDensity());
        m_vecKinPtrs[n]->getNetProductionRates(DATA_PTR(m_work));
        kstart = m_vecKinPtrs[n]->kineticsSpeciesIndex(0,m_surfindex[n]);
        sum = 0.0;
        loc = 0;
        for (size_t k = 1; k < m_nsp_surf[n]; k++) {
            ydot[k + loc] = m_work[kstart + k] * rs0 * m_surf[n]->size(k);
            sum -= ydot[k];
        }
        ydot[loc] = sum;

        m_transport->getMixDiffCoeffs(m_diff_coeffs.at(n));
        for (size_t k = 0; k < (m_nsp_tot[n]-m_nsp_surf[n]); k++) {
            Sc = m_transport->viscosity() /m_bulkPhases[n]->density() / m_diff_coeffs.at(n)[k];
            ydot[k + m_nsp_surf[n] +  loc] = m_work[k] * m_bulkPhases[n]->molecularWeight(k)
                          / m_bulkPhases[n]->density() 
                         - m_masstransfer_coefficient * m_diff_coeffs.at(n)[k]
                                                      * (y[k+m_nsp_surf[n]] - 
                                                         m_bulk_massfraction[n][k]) ;
        }
        loc += m_nsp_tot[n];
        m_surf[n]->setSiteDensity(site_density_save);
    }
}

/*
 * getConcSpecies():
 *
 * Fills the local concentration vector, m_concSpecies for all of the
 * species in all of the phases that are unknowns in the surface
 * problem.
 *
 * m_concSpecies[]
 */
void ImplicitSurfChem_masstransfer::getConcSpecies(doublereal* const vecConcSpecies) const
{
    size_t kstart;
    for (size_t ip = 0; ip < m_nsurf; ip++) {
        ThermoPhase* TP_ptr = m_surf[ip];
        kstart = m_specStartIndex[ip];
        TP_ptr->getConcentrations(vecConcSpecies + kstart);
    }
    for (size_t ip = 0; ip <  m_numBulkPhases; ip++) {
        ThermoPhase* TP_ptr = m_bulkPhases[ip];
        TP_ptr->getConcentrations(vecConcSpecies + kstart);
        kstart += TP_ptr->nSpecies();
    }
}

/*
 * setConcSpecies():
 *
 * Fills the local concentration vector, m_concSpecies for all of the
 * species in all of the phases that are unknowns in the surface
 * problem.
 *
 * m_concSpecies[]
 */
void ImplicitSurfChem_masstransfer::setConcSpecies(const doublereal* const vecConcSpecies)
{
    size_t kstart;
    for (size_t ip = 0; ip < m_nsurf; ip++) {
        ThermoPhase* TP_ptr = m_surf[ip];
        kstart = m_specStartIndex[ip];
        TP_ptr->setConcentrations(vecConcSpecies + kstart);
    }
    kstart = m_nv;
    for (size_t ip = 0; ip <  m_numBulkPhases; ip++) {
        ThermoPhase* TP_ptr = m_bulkPhases[ip];
        TP_ptr->setConcentrations(vecConcSpecies + kstart);
        kstart += TP_ptr->nSpecies();
    }
}

/*
 * setCommonState_TP():
 *
 *  Sets a common temperature and pressure amongst the
 *  thermodynamic objects in the interfacial kinetics object.
 *
 *  Units Temperature = Kelvin
 *        Pressure    = Pascal
 */
void ImplicitSurfChem_masstransfer::
setCommonState_TP(doublereal TKelvin, doublereal PresPa)
{
    for (size_t ip = 0; ip < m_nsurf; ip++) {
        ThermoPhase* TP_ptr = m_surf[ip];
        TP_ptr->setState_TP(TKelvin, PresPa);
    }
    for (size_t ip = 0; ip < m_numBulkPhases; ip++) {
        ThermoPhase* TP_ptr = m_bulkPhases[ip];
        TP_ptr->setState_TP(TKelvin, PresPa);
    }
}

}
