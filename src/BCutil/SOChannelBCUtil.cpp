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

//-------------------------------------------------------------------------------
//basicVelFuncBC
//-------------------------------------------------------------------------------
BCMethodHolder SOChannelBCUtil::basicVelFuncBC (int a_veldir, bool a_isViscous) const
{
    // // TEMPORARY!!
    // return PhysBCUtil::basicVelFuncBC(a_veldir, a_isViscous);

    BCMethodHolder holder;
    const int extrapOrder = 0;

    //Tried also looking at BeamGenerationBCUtil.cpp -- seems to have decomposition
    //into velocity components in a similar manner

    if (a_veldir == 0) {
        // x-dir
        {
            // Periodic. Do nothing
        }

        // y-dir
        {
            // Lo side
            if (a_isViscous) {
                // No-slip
                RefCountedPtr<BCGhostClass> ghostBCPtr(
                    new EllipticConstDiriBCGhostClass(RealVect::Zero,
                                                      RealVect::Zero,
                                                      BASISV(1),
                                                      IntVect::Zero)
                );
                holder.addBCMethod(ghostBCPtr);

            } else {
                // Free-slip
                RefCountedPtr<BCGhostClass> ghostBCPtr(
                    new EllipticExtrapBCGhostClass(extrapOrder,
                                                   BASISV(1),
                                                   IntVect::Zero)
                );
                holder.addBCMethod(ghostBCPtr);

                RefCountedPtr<BCFluxClass> fluxBCPtr(
                    new EllipticConstNeumBCFluxClass(RealVect::Zero,
                                                     RealVect::Zero,
                                                     BASISV(1),
                                                     IntVect::Zero)
                );
                holder.addBCMethod(fluxBCPtr);
            }

            // Hi side
            {
                // Free-slip
                RefCountedPtr<BCGhostClass> ghostBCPtr(
                    new EllipticExtrapBCGhostClass(extrapOrder,
                                                   IntVect::Zero,
                                                   BASISV(1))
                );
                holder.addBCMethod(ghostBCPtr);

                RefCountedPtr<BCFluxClass> fluxBCPtr(
                    new EllipticConstNeumBCFluxClass(RealVect::Zero,
                                                     RealVect::Zero,
                                                     IntVect::Zero,
                                                     BASISV(1))
                );
                holder.addBCMethod(fluxBCPtr);
            }
        } // end y-dir

        // z-dir
        if (a_isViscous) {
            //  No-slip
            RefCountedPtr<BCGhostClass> ghostBCPtr(
                new EllipticConstDiriBCGhostClass(RealVect::Zero,
                                                  RealVect::Zero,
                                                  BASISV(2),
                                                  BASISV(2))
            );
            holder.addBCMethod(ghostBCPtr);

        } else {
            // Free-slip
            RefCountedPtr<BCGhostClass> ghostBCPtr(
                new EllipticExtrapBCGhostClass(extrapOrder,
                                               BASISV(2),
                                               BASISV(2))
            );
            holder.addBCMethod(ghostBCPtr);

            RefCountedPtr<BCFluxClass> fluxBCPtr(
                new EllipticConstNeumBCFluxClass(RealVect::Zero,
                                                 RealVect::Zero,
                                                 BASISV(2),
                                                 BASISV(2))
            );
            holder.addBCMethod(fluxBCPtr);
        } // end z-dir

    } else if (a_veldir == 1) {
        // TEMPORARY
        // We will set standard no-slip/free-slip BCs all around,
        // then overwrite them one-by-one. This way, we can comment
        // out some BCs and have something to fall back on.
        holder = PhysBCUtil::basicVelFuncBC(a_veldir, a_isViscous);

        // x-dir
        {
            // Periodic. Do nothing
        }

        // y-dir
        {
            // Lo side
            {
                // No momentum flux through wall
                RefCountedPtr<BCGhostClass> ghostBCPtr(
                    new EllipticConstDiriBCGhostClass(RealVect::Zero,
                                                      RealVect::Zero,
                                                      BASISV(1),
                                                      IntVect::Zero)
                );
                holder.addBCMethod(ghostBCPtr);
            }

            // NOTE: This leads to an instability.
            // I think we need zero Diri Pressure BCs here.
            // But that may conflict with the other velocity comp BCs.
            // This needs some thought.
            // // Hi side
            // {
            //     // Zero strain rate normal to boundary.
            //     RefCountedPtr<BCGhostClass> ghostBCPtr(
            //         new EllipticExtrapBCGhostClass(extrapOrder,
            //                                        IntVect::Zero,
            //                                        BASISV(1))
            //     );
            //     holder.addBCMethod(ghostBCPtr);

            //     RefCountedPtr<BCFluxClass> fluxBCPtr(
            //         new EllipticConstNeumBCFluxClass(RealVect::Zero,
            //                                          RealVect::Zero,
            //                                          IntVect::Zero,
            //                                          BASISV(1))
            //     );
            //     holder.addBCMethod(fluxBCPtr);
            // }
        } // end y-dir

        // z-dir
        if (a_isViscous) {
            //  No-slip
            RefCountedPtr<BCGhostClass> ghostBCPtr(
                new EllipticConstDiriBCGhostClass(RealVect::Zero,
                                                  RealVect::Zero,
                                                  BASISV(2),
                                                  BASISV(2))
            );
            holder.addBCMethod(ghostBCPtr);

        } else {
            // Free-slip
            RefCountedPtr<BCGhostClass> ghostBCPtr(
                new EllipticExtrapBCGhostClass(extrapOrder,
                                               BASISV(2),
                                               BASISV(2))
            );
            holder.addBCMethod(ghostBCPtr);

            RefCountedPtr<BCFluxClass> fluxBCPtr(
                new EllipticConstNeumBCFluxClass(RealVect::Zero,
                                                 RealVect::Zero,
                                                 BASISV(2),
                                                 BASISV(2))
            );
            holder.addBCMethod(fluxBCPtr);
        } // end z-dir

    } else {
        // No slip / free slip on all dirs and sides.
        return PhysBCUtil::basicVelFuncBC(a_veldir, a_isViscous);
    }

    return holder;
}
