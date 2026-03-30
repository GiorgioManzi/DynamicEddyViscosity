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
    Foam::RASModels::kEpsDyn

Description
    Abe, Kondoh and Nagano's low-Reynolds number k-epsilon turbulence model
    for incompressible flows
	
    Wall boundary conditions are:

        k       = 1e-10
        epsilon = epsilonLowReWallFunction
	
	
	References:

    Abe K., Kondoh T., Nagano Y., 1994. A new turbulence model for predicting
	fluid flow and heat transfer in separating and reattaching flows-I. Flow
	field calculations. Int. J. Heat and Mass Transfer 37(1), pp 138-151

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
    kEpsDyn.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsDyn_H
#define kEpsDyn_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class kEpsDyn Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kEpsDyn
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
			const int finestra_;
			IOdictionary physicalProperties_;
			//dimensionedScalar nu_;
			dimensionedScalar low_epsilon_;
			dimensionedScalar low_k_;
			dimensionedScalar low_M_;


        // Fields
			
            volScalarField k_;
			volScalarField dynamicCmu_;
			volSymmTensorField S;
			//volTensorField Sns;
			volTensorField gU;
            volScalarField epsilon_;
			
			
			//const volScalarField::Internal& y_;
			const volScalarField& y_;
			
			volScalarField f2Field_;
			
			PtrList<volVectorField> UPrev_;
			PtrList<volScalarField> kPrev_;
	        PtrList<volScalarField> epsilonPrev_;
			PtrList<volVectorField> Uprint_;
			PtrList<volScalarField> kprint_;
	        PtrList<volScalarField> epsilonprint_;
			PtrList<volTensorField> gUprint_;
			
			


    // Private Member Functions

        //tmp<volScalarField> fMu() const;
        tmp<volScalarField> f2() const;

        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("kEpsDyn");


    // Constructors

        //- Construct from components
        kEpsDyn
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        kEpsDyn(const kEpsDyn&) = delete;


    //- Destructor
    virtual ~kEpsDyn()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                (this->nut_/sigmak_ + this->nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return volScalarField::New
            (
                "DepsilonEff",
                (this->nut_/sigmaEps_ + this->nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return volScalarField::New
            (
                "omega",
                epsilon_/(dynamicCmu_*k_),
                epsilon_.boundaryField().types()
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const kEpsDyn&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kEpsDyn.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
