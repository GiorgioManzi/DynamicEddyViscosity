/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::RASModels::kEpsDynamic

Description
    Dynamic k-epsilon turbulence model based on the Abe, Kondoh and Nagano's low-Reynolds number formulation.
	
	References:

    References: 
    \verbatim 
        Cimarelli, A., Marras, N., Niceno, B., & Tessier Urrecha, Y. (2025).
        Unifying formalism and closures for coarse-grained approaches to turbulence.
        arXiv preprint arXiv:2512.19274.

        Tessier Urrecha, Y., Manzi, G., Cancellotti, F., & Cimarelli, A. (2026).
        Development of dynamic RANS closures for coarse-grained approaches to 
        turbulence. Under review.
    
        Abe K., Kondoh T., Nagano Y., 1994. A new turbulence model for predicting
	    fluid flow and heat transfer in separating and reattaching flows-I. Flow
	    field calculations. Int. J. Heat and Mass Transfer 37(1), pp 138-151.

    The default model coefficients are given as:
    \verbatim
        AbeKondohNaganoKECoeffs
        {
            Cmu         0.09;
            C1          1.5;
            C2          1.9;
            sigmaEps    1.4;
            sigmaK      1.4;
        }
    \endverbatim

SourceFiles
    kEpsDynamic.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsDynamic_H
#define kEpsDynamic_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class kEpsDynamic Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kEpsDynamic
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel>>
{
protected:

    // Protected data

        // Model coefficients

        //dimensionedScalar Cmu_;
        dimensionedScalar C1_;
        dimensionedScalar C2_;
        dimensionedScalar C3_;
        dimensionedScalar sigmak_;
        dimensionedScalar sigmaEps_;

        // dynamicCmu clipping limits
        dimensionedScalar dynamicCmuMin_;
        dimensionedScalar dynamicCmuMax_;

        label nStart_;

		const int window_;
		IOdictionary physicalProperties_;
			
		dimensionedScalar low_epsilon_;
		dimensionedScalar low_k_;
		dimensionedScalar low_M_;

        // Fields
			
        volScalarField k_;
        volScalarField epsilon_;

		volSymmTensorField S;
		volTensorField gU;
        volScalarField dynamicCmu_;

        // Low-Re damping function
		const volScalarField& y_;	
		volScalarField f2Field_;
		
        //Dynamic procedure buffer
		PtrList<volVectorField> UPrev_;
		PtrList<volScalarField> kPrev_;
	    PtrList<volScalarField> epsilonPrev_;
        PtrList<volSymmTensorField> SPrev_;
		PtrList<volVectorField> Uprint_;
		PtrList<volScalarField> kprint_;
	    PtrList<volScalarField> epsilonprint_;
		PtrList<volTensorField> gUprint_;
			
        label ringIndex_;

    // Private Member Functions

        //tmp<volScalarField> fMu() const;
        tmp<volScalarField> f2() const;

        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;

    //- Runtime type info
    TypeName("kEpsDynamic");

    // Constructors

        //- Construct from components
        kEpsDynamic
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = typeName
        );

        // Do not allow copy
        kEpsDynamic(const kEpsDynamic&) = delete;

    // Destructor
    virtual ~kEpsDynamic(){}

    // Member Functions

        // Read coeff. 
        virtual bool read();

        // Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New("DkEff",(this->nut_/sigmak_ + this->nu()));
        }

        // Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return volScalarField::New("DepsilonEff",(this->nut_/sigmaEps_ + this->nu()));
        }

        // Return the turbulent kinetic energy
        virtual tmp<volScalarField> k() const { return k_; }

        // Return the turbulence kinetic energy dissipation
        virtual tmp<volScalarField> epsilon() const{return epsilon_;}

        // Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return volScalarField::New
            (
                "omega",
                epsilon_/(dynamicCmu_*k_),
                epsilon_.boundaryField().types()
            );
        }

        // Solve equations and correct nut
        virtual void correct();

    // Member Operators
        // Do not allow assignment
        void operator=(const kEpsDynamic&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace RASModels
} // namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    # include "kEpsDynamic.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //