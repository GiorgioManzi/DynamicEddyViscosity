/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "kEpsDynamic.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsDynamic<BasicMomentumTransportModel>::f2() const
{
	const volScalarField yStar(pow(this->nu()*epsilon_,scalar(0.25))*y_/this->nu()); 
	
    volScalarField eps(epsilon_);
    tmp<volScalarField> Rt = sqr(k_)/(this->nu()*bound(eps, this->epsilonMin_));	
	
    return
        min((scalar(1)-0.3*exp(-sqr(Rt/6.5)))*sqr(scalar(1)-exp(-yStar/3.1)),scalar(1.0));
}


template<class BasicMomentumTransportModel>
void kEpsDynamic<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = dynamicCmu_*sqr(k_)/(epsilon_+low_epsilon_);
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>kEpsDynamic<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>(
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>kEpsDynamic<BasicMomentumTransportModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>(
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()/dimTime
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kEpsDynamic<BasicMomentumTransportModel>::kEpsDynamic
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.5
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.9
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.4
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.4
        )
    ),
	window_
    (
        this->coeffDict().lookupOrDefault("nWindow", 2)
    ),

	physicalProperties_
    (
        IOobject
        (
            "physicalProperties",         
            this->runTime_.constant(),    
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    nStart_
    (
        this->coeffDict().template lookupOrDefault<label>("nStart", 10)
    ),
	
	low_epsilon_
	(
	    "low_epsilon",
		dimensionSet(0,2,-3,0,0),
		SMALL
	),

	low_k_
	(
	    "low_k",
		dimensionSet(0, 2, -2, 0, 0),
		SMALL
	),

	low_M_
	(
	    "low_M",
		dimensionSet(0, 4, -4, 0, 0),
		SMALL
	),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    dynamicCmu_
    (
        IOobject
        (
            "dynamicCmu",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
	S
    (
        IOobject
        (
            "S",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
		dimensionedSymmTensor(dimensionSet(0, 0, -1, 0, 0), Zero)
    ),

	gU
    (
        IOobject
        (
            "gU",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
		dimensionedTensor(dimensionSet(0, 0, -1, 0, 0), Zero)
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
		
	y_(wallDist::New(this->mesh_).y()),
	
	f2Field_
    (
        IOobject
        (
            "f2",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        f2()()
    ),

    ringIndex_(0)

{
    Info << "kOmegaDynamic constructor called!" << nl;

////////////////////////////////    PtrLists    ////////////////////////////////////

	kPrev_.setSize(window_ - 1);
	epsilonPrev_.setSize(window_ - 1);
	UPrev_.setSize(window_ - 1);
    SPrev_.setSize(window_ - 1);
	kprint_.setSize(window_ - 1);
	epsilonprint_.setSize(window_ - 1);
	Uprint_.setSize(window_ - 1);
	gUprint_.setSize(window_ - 1);

	for (int i = 0; i < window_ - 1; ++i)
    {
		kPrev_.set
		(
		    i,
			new volScalarField
			(
			    IOobject
                (
				    "kPrev_" + name(i),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar(dimensionSet(0,2,-2,0,0), 0)
			)
		);

		epsilonPrev_.set
		(
		    i,
			new volScalarField
			(
			    IOobject
                (
				    "epsilonPrev_" + name(i),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar(dimensionSet(0,2,-3,0,0), SMALL)
			)
		);

		UPrev_.set
		(
		    i,
			new volVectorField
			(
			    IOobject
                (
				    "UPrev_" + name(i),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedVector(dimensionSet(0,1,-1,0,0), Zero)
			)
		);

        SPrev_.set
        (
            i,
            new volSymmTensorField
            (
                IOobject
                (
                    "SPrev_" + name(i),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedSymmTensor(dimensionSet(0,0,-1,0,0), Zero)
            )
        );

		kprint_.set
		(
		    i,
			new volScalarField
			(
			    IOobject
                (
				    "kprint_" + name(i+1),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar(dimensionSet(0,2,-2,0,0), 0)
			)
		);

		epsilonprint_.set
		(
		    i,
			new volScalarField
			(
			    IOobject
                (
				    "epsilonprint_" + name(i+1),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar(dimensionSet(0,2,-3,0,0), 1e-12)
			)
		);

		Uprint_.set
		(
		    i,
			new volVectorField
			(
			    IOobject
                (
				    "Uprint_" + name(i+1),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedVector(dimensionSet(0,1,-1,0,0), Zero)
			)
		);

		gUprint_.set
		(
		    i,
			new volTensorField
			(
			    IOobject
                (
				    "gUprint_" + name(i+1),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedTensor(dimensionSet(0,0,-1,0,0), Zero)
			)
		);
    }

    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool kEpsDynamic<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicMomentumTransportModel>
void kEpsDynamic<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local ref.
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;

    volScalarField& nut = this->nut_;

    const Foam::fvModels& fvModels
    (
        Foam::fvModels::New(this->mesh_)
    );

    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    volScalarField divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );

    tmp<volTensorField> tgradU = fvc::grad(U);

    volScalarField G
    (
        this->GName(), 
        nut*(tgradU() && dev(twoSymm(tgradU())))
    );

        // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha*rho*G*epsilon_/(k_+low_k_)
      - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha*rho*divU, epsilon_)
      - fvm::Sp(C2_*f2()*alpha*rho*epsilon_/(k_+low_k_), epsilon_)
      + epsilonSource()
      + fvModels.source(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvConstraints.constrain(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G 
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(alpha*rho*epsilon_/(k_+low_k_), k_)
      + kSource()
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);
	
//////////////////////////////////  Dynamic procedure  ////////////////////////////////////

	const volTensorField& gradU = tgradU();

	S = symm(gradU);

	gU = gradU;
	
	volSymmTensorField sum_uu = symm(U * U);  
	volVectorField sum_u = U;
	
	volSymmTensorField kES_mean = 2 * sqr(k_) / (epsilon_ + low_epsilon_) * S;    
	volSymmTensorField sum_S = S;
	
	volScalarField k_mean = k_;
	volScalarField epsilon_mean = epsilon_;
	
	volScalarField squaredS_mean = 
        sqr(S.component(symmTensor::XX))
      + sqr(S.component(symmTensor::YY))
      + sqr(S.component(symmTensor::ZZ));

	tgradU.clear();
	
	for (int i = 0; i < window_ - 1; ++i)
    {
        sum_uu += symm(UPrev_[i] * UPrev_[i]);    
		sum_u += UPrev_[i];

        const volSymmTensorField& S_i = SPrev_[i];
		
		kES_mean += 2 * sqr(kPrev_[i]) / (epsilonPrev_[i] + low_epsilon_) * S_i;
		sum_S += S_i;
		
		k_mean += kPrev_[i];
		epsilon_mean += epsilonPrev_[i];
		
        squaredS_mean +=
            sqr(S_i.component(symmTensor::XX))
          + sqr(S_i.component(symmTensor::YY))
          + sqr(S_i.component(symmTensor::ZZ));
    }
	
	sum_u = sum_u / window_;
	volSymmTensorField uu_average = sum_uu / window_;
	volSymmTensorField u_average_u_average = symm((sum_u) * (sum_u));
		
	volSymmTensorField L = uu_average - u_average_u_average;
	
	kES_mean = kES_mean / window_;
	volSymmTensorField S_average = sum_S / window_;

	k_mean = k_mean / window_;
	k_mean = k_mean + 0.5 * (tr(uu_average) - (sum_u & sum_u));
	
	epsilon_mean = epsilon_mean / window_;
	
	squaredS_mean = squaredS_mean / window_;
	
	volScalarField SSquared_mean = 
        sqr(S_average.component(symmTensor::XX))
	  + sqr(S_average.component(symmTensor::YY))
      + sqr(S_average.component(symmTensor::ZZ));

	epsilon_mean = epsilon_mean + this->nu() * (squaredS_mean - SSquared_mean);	
	
	volSymmTensorField M = kES_mean - 2*sqr(k_mean)/(epsilon_mean + low_epsilon_)*S_average;
	
	dynamicCmu_ = (M && L) / ((M && M)+low_M_);       

    // Starting procedure 
	if (this->runTime_.timeIndex() < nStart_)
    {
        dynamicCmu_ = dimensionedScalar("cmuStart", dimless, 0.09);
    }

    // Clipping
    dynamicCmu_ = max(dynamicCmuMin_, min(dynamicCmu_, dynamicCmuMax_));

    // Update stored fields
    Uprint_[ringIndex_] = UPrev_[ringIndex_];
    kprint_[ringIndex_] = kPrev_[ringIndex_];
    epsilonprint_[ringIndex_] = epsilonPrev_[ringIndex_];
    gUprint_[ringIndex_] = fvc::grad(UPrev_[ringIndex_]);

    // Insert new values
    UPrev_[ringIndex_] = this->U_;
    kPrev_[ringIndex_] = this->k_;
    epsilonPrev_[ringIndex_] = this->epsilon_;
    SPrev_[ringIndex_] = S;

    // Advance index
    ringIndex_ = (ringIndex_ + 1) % (window_ - 1);
	
    correctNut();
	
    Info << "Time = " << this->runTime_.timeName() << nl << endl;
}

////////////////////////////////////////////////////////////////////////////////

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //