/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2014 University of North Carolina at Chapel Hill
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *  USA
 *
 *  For up-to-date contact information, please visit the repository homepage,
 *  https://github.com/somarhub.
 ******************************************************************************/
#include "SOChannelBCUtil.H"
#include "EllipticBCUtils.H"
#include "Constants.H"
#include "LevelGeometry.H"


// -----------------------------------------------------------------------------
// Default constructor
// -----------------------------------------------------------------------------
SOChannelBCUtil::SOChannelBCUtil ()
{
    // Nothing is needed here. Setup is done in PhysBCUtil constructor.
}


// -----------------------------------------------------------------------------
// Default destructor
// -----------------------------------------------------------------------------
SOChannelBCUtil::~SOChannelBCUtil ()
{;}


// -----------------------------------------------------------------------------
// This object is its own factory
// -----------------------------------------------------------------------------
PhysBCUtil* SOChannelBCUtil::newPhysBCUtil () const
{

    PhysBCUtil* newBCPtr = new SOChannelBCUtil();
    return newBCPtr;
}


// -----------------------------------------------------------------------------
// Fills a FAB with the initial scalars
// -----------------------------------------------------------------------------
void SOChannelBCUtil::setScalarIC (FArrayBox&           a_scalarFAB,
                                   const int            a_scalarComp,
                                   const LevelGeometry& a_levGeo,
                                   const DataIndex&     a_di) const
{
    a_scalarFAB.setVal(0.0, a_scalarComp);
}


// -----------------------------------------------------------------------------
// Holds Diri BCs for the scalar field at the top of the domain.
// This is just a linear profile in b_{total}.
// -----------------------------------------------------------------------------
void SOChannelBCUtil::topBCValueFunc (Real*           a_pos,
                                      int*            a_dir,
                                      Side::LoHiSide* a_side,
                                      Real*           a_value,
                                      const RealVect& a_dx,
                                      Real            a_time)
{
    // We had better be at the top!
    CH_assert(a_dir[0] == SpaceDim-1);
    CH_assert(a_side[0] == Side::Hi);

    // These params adjust the shape of the hot/cold interface.
    static const RealVect L = LevelGeometry::getDomainLength();
    static const Real minB = 0.0;
    static const Real deltaB = 1.0;

    Real fracX = a_pos[1] / L[1];
    a_value[0] = minB + deltaB*fracX;
}


// -----------------------------------------------------------------------------
// diffusiveSolveFuncBC (used in scalar TGA solves)
// This is where the adiabatic/hot/cold BCs are set.
// -----------------------------------------------------------------------------
BCMethodHolder SOChannelBCUtil::diffusiveSolveFuncBC () const
{
    BCMethodHolder holder;

    const IntVect vmask = BASISV(SpaceDim-1);
    const IntVect hmask = IntVect::Unit - vmask;

    // Adiabatic ghost BC
    RefCountedPtr<BCGhostClass> adGhostBCPtr (
        new EllipticConstNeumBCGhostClass(RealVect::Zero,
                                          RealVect::Zero,
                                          IntVect::Unit,
                                          hmask)
    );
    holder.addBCMethod(adGhostBCPtr);

    // Adiabatic flux BC
    RefCountedPtr<BCFluxClass> adFluxBCPtr (
        new EllipticConstNeumBCFluxClass(RealVect::Zero,
                                         RealVect::Zero,
                                         IntVect::Unit,
                                         hmask)
    );
    holder.addBCMethod(adFluxBCPtr);

    // Top BC
    RefCountedPtr<BCGhostClass> topGhostBCPtr (
        new EllipticDiriBCGhostClass(&topBCValueFunc,
                                     IntVect::Zero,
                                     vmask)
    );
    holder.addBCMethod(topGhostBCPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// diffusiveSourceFuncBC (Used to calculate the diffusive term nu.L[scalar])
// Just use the diffusive solver's BCs.
// -----------------------------------------------------------------------------
BCMethodHolder SOChannelBCUtil::diffusiveSourceFuncBC () const
{
    return SOChannelBCUtil::diffusiveSolveFuncBC();
}


BCMethodHolder SOChannelBCUtil::basicVelFuncBC (int a_veldir, bool a_isViscous) const
{
   BCMethodHolder holder;
   const int extrapOrder = 0;
   IntVect ymask(BASISV(1));
   IntVect nonymask = IntVect::Unit - ymask;
   
   //Tried also looking at BeamGenerationBCUtil.cpp -- seems to have decomposition
   //into velocity components in a similar manner
   
   if (a_veldir == 1) {
        if (a_isViscous) {
           // zero order extrap on high y side: dv/dy=0
            RefCountedPtr<BCGhostClass> yextrapBCPtr(
                new EllipticExtrapBCGhostClass(extrapOrder,
                                               IntVect::Zero,
                                               ymask)
            );
            holder.addBCMethod(yextrapBCPtr);
            
           // no flux BCs on high y side: dv/dy = 0
            RefCountedPtr<BCFluxClass> yfluxBCPtr (
                new EllipticConstNeumBCFluxClass(RealVect::Zero,
                                                 RealVect::Zero,
                                                 IntVect::Zero,
                                                 ymask)
            );
            holder.addBCMethod(yfluxBCPtr);
            
           // no slip, solid wall on all other sides: v=0
            RefCountedPtr<BCGhostClass> vBCPtr = RefCountedPtr<BCGhostClass>(
                new BasicVelocityBCGhostClass(0.0,             // inflowVel
                                              -1,              // inflowDir
                                              Side::Lo,        // inflowSide
                                              -1,              // outflowDir
                                              Side::Hi,        // outflowSide
                                              a_veldir,
                                              a_isViscous,
                                              IntVect::Unit,
                                              nonymask)
            );
            holder.addBCMethod(vBCPtr);       
            
        } else {
           // zero order extrap on high y side
            RefCountedPtr<BCGhostClass> yextrapBCPtr(
                new EllipticExtrapBCGhostClass(extrapOrder,
                                               IntVect::Zero,
                                               ymask)
            );
            holder.addBCMethod(yextrapBCPtr);
           // no flux BCs on high y side
            RefCountedPtr<BCFluxClass> yfluxBCPtr (
                new EllipticConstNeumBCFluxClass(RealVect::Zero,
                                                 RealVect::Zero,
                                                 IntVect::Zero,
                                                 ymask)
            );
            holder.addBCMethod(yfluxBCPtr);
            
           // free slip, solid wall all other sides
            RefCountedPtr<BCGhostClass> vslipBCPtr = RefCountedPtr<BCGhostClass>(
                new BasicVelocityBCGhostClass(0.0,             // inflowVel
                                              -1,              // inflowDir
                                              Side::Lo,        // inflowSide
                                              -1,              // outflowDir
                                              Side::Hi,        // outflowSide
                                              a_veldir,
                                              false,
                                              IntVect::Unit,
                                              nonymask)
            );
            holder.addBCMethod(vslipBCPtr);
        }
   } else {
        return PhysBCUtil::basicVelFuncBC (a_veldir, a_isViscous);
   }

   return holder;
}





// Code below this point has not been refactored for SOChannel problem yet.
#if 0

// -----------------------------------------------------------------------------
// viscousSolveFuncBC (Used in single-component velocity TGA solves)
// -----------------------------------------------------------------------------
BCMethodHolder SOChannelBCUtil::viscousSolveFuncBC (int a_dir) const
{
    BCMethodHolder holder;

    const IntVect normVect = BASISV(a_dir);
    const IntVect transVect = IntVect::Unit - normVect;

    if (a_dir == 1)
    {
        // Almost! a_dir is the velocity component, so the code inside this block
        // only sets BCs on v. I guess I shoulc have called a_dir something like
        // a_velComp instead for clarity. -ES
        RefCountedPtr<BCGhostClass> ytransGhostBCPtr (
            new EllipticConstNeumBCGhostClass(RealVect::Zero,
                                              RealVect::Zero,
                                              IntVect(1,0,1), //on low side want dw/dy, du/dy=0
                                              IntVect(1,1,0)) //on high side want du/dy,dv/dy=0
        );
        holder.addBCMethod(ytransGhostBCPtr);


        RefCountedPtr<BCFluxClass> ytransFluxBCPtr (
            new EllipticConstNeumBCFluxClass(RealVect::Zero,
                                             RealVect::Zero,
                                             IntVect(1,0,1), //on low side want dw/dy, du/dy=0
                                             IntVect(1,1,0)) //on high side want du/dy,dv/dy=0
        );
        holder.addBCMethod(ytransFluxBCPtr);
        RefCountedPtr<BCGhostClass> ynormGhostBCPtr (
            new EllipticConstDiriBCGhostClass(RealVect::Zero,
                                              RealVect::Zero,
                                              IntVect(0,1,0), //on low side want v=0
                                              IntVect(0,0,1)) //on high side want w=0
        );
        holder.addBCMethod(ynormGhostBCPtr);

    }
    else
    {
        // Transverse, no-slip BCs
        RefCountedPtr<BCGhostClass> transGhostBCPtr (
            new EllipticConstNeumBCGhostClass(RealVect::Zero,
                                              RealVect::Zero,
                                              transVect,
                                              transVect)
        );
        holder.addBCMethod(transGhostBCPtr);

        // Transverse, no-slip BCs (sets fluxes)
        RefCountedPtr<BCFluxClass> transFluxBCPtr (
            new EllipticConstNeumBCFluxClass(RealVect::Zero,
                                             RealVect::Zero,
                                             transVect,
                                             transVect)
        );
        holder.addBCMethod(transFluxBCPtr);

        // Normal, no flux BCs
        RefCountedPtr<BCGhostClass> normGhostBCPtr (
            new EllipticConstDiriBCGhostClass(RealVect::Zero,
                                              RealVect::Zero,
                                              normVect,
                                              normVect)
        );
        holder.addBCMethod(normGhostBCPtr);
    }

    return holder;
}

#endif // 0

